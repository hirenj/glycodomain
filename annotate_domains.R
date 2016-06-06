if ( ! exists('interpro_raw_relations') ) {
  download.file('ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt','interpro_relations.txt')
  interpro_raw_relations = gsub('::.*','',readLines('interpro_relations.txt'))
}

if ( ! exists('interpro_raw_go')) {
  download.file('http://www.geneontology.org/external2go/interpro2go','interpro_go.txt')
  interpro_raw_go = gsub('^InterPro:','',Filter(function(x) grepl("^InterPro",x),readLines('interpro_go.txt')))
  interpro_raw_go = as.data.frame(matrix( unlist(lapply(interpro_raw_go,function(x) { lapply(strsplit(x,' '), function(vals) { vals[c(1,length(vals))] }) })), ncol=2,byrow=T))
  interpro_raw_go = setNames(interpro_raw_go,c('interpro','go'))
}

go_manual_mappings = read.delim('go_manual.tsv')
interpro_manual_mappings = read.delim('interpro_manual.tsv')

if ( ! exists('all_sites.human')) {
  devtools::source_gist('b79cb9f63bf6d7c6cb68')
  gator.flatpack <- as.data.frame(gator.flatpack)
  msdatas = do.call(rbind,lapply(grep("msdata",ls(envir = Rgator:::getDataEnvironment()),value=T), function(set) {
    datas = get(set)
    datas = datas[datas$multiple.protein == FALSE & grepl("[STY]",datas$peptide),]
    hexnac_ratios = Rgator:::hexnac.ratios(datas)
    if (! is.null(hexnac_ratios) && nrow(hexnac_ratios) > 0) {
      hexnac_ratios = lapply(split(hexnac_ratios$ratio,hexnac_ratios$spectra),function(x) median(x,na.rm=T) )
      datas = datas[datas$spectra %in% names(Filter(function(median) { median <= 1.1 },hexnac_ratios)),]
    }
    if (nrow(datas) < 1) {
      return(datas)
    }
    datas$biosrc = gsub(".msdata.json","",set,fixed=T)
    datas$source.free = datas$source
    datas$compartment[grepl("tcl",datas$source.free,ignore.case=T)] <- 'tcl'
    datas$compartment[grepl("sec",datas$source.free,ignore.case=T)] <- 'sec'
    datas[,c('uniprot','site','peptide.start','peptide.end','biosrc','source.free','compartment')]
  }))

  gator.flatpack = gator.flatpack[! gator.flatpack$biosrc %in% c('fibroblast_hcmv','fibroblast_hcmv2','fibroblast_vzv','fibroblast_vzv2','hacat_hsv2'),]
  gator.flatpack = rbind(gator.flatpack,setNames(msdatas,names(gator.flatpack)))

  gator.flatpack$site = ifelse(is.na(gator.flatpack$sites.site),floor((as.numeric(gator.flatpack$peptide.start)+as.numeric(gator.flatpack$peptide.end))/2),gator.flatpack$sites.site )
  all_sites = unique(gator.flatpack[!is.na(gator.flatpack$site),c('uniprot','site')])
  all_sites.human = all_sites[all_sites$uniprot %in% tolower(Rgator::getUniprotIds(9606)),]
}

if ( is.null(get('transmembrane.uniprot.9606')) ) {
  Rgator:::downloadTransmembrane(9606)
}
if ( is.null(get('interpro.9606')) ) {
  Rgator::downloadInterproDomains(9606)
}

human.domains = rbind( interpro.9606.annotated[interpro.9606.annotated$class == 'Domain',c('uniprot','dom','start','end')], transmembrane.uniprot.9606 )


get_classes = function(interpro) {
  go_mapped = merge(interpro_raw_go[interpro_raw_go$interpro %in% interpro,], go_manual_mappings,by='go')[,c('interpro','Class')]
  manual_mapped = interpro_manual_mappings[interpro_manual_mappings$interpro %in% interpro,c('interpro','Class')]
  all_mapped = rbind(go_mapped,manual_mapped)
  if (length(unique(all_mapped$interpro)) == length(unique(interpro))) {
    return (all_mapped)
  }
  return (rbind(all_mapped, data.frame(interpro=interpro[!interpro %in% all_mapped$interpro],Class=NA)))
}


group_interpros = function(interpros,depth=4) {
	grouped = split(interpros,cumsum(seq_along(interpros) %in% grep("^I",interpros)))
	group_names = sapply(grouped,function(group) { group[1] })
	grouped = lapply(grouped, function(x) {  if (length(x) > 1) { gsub('^--','',x[2:length(x)]) } else { NA }  })
	names(grouped) = group_names
	if ( depth > 0 ) {
		lapply(grouped, function(group) { if (is.na(group[1])) { group } else { group_interpros(group,depth - 1) } })
	} else {
		grouped
	}
}

parent_terms = function(term,groups,distance=1) {
  val = unique(unlist(lapply(
    strsplit(
      gsub(paste(term,".*$",sep=""),'',grep(term,names(unlist(groups)),value=T)),
      '.',
      fixed=T),
      function(x) { distval = (length(x)-distance+1); x[ifelse(distval < 0, 0,distval):length(x)] }
    )))
  if (is.null(val)) {
    NA
  } else {
    val
  }
}

child_terms = function(term,groups,distance=1) {
  val = Filter(
    function(x) !is.na(x) & x != '',
    unique(unlist(lapply(
    strsplit(
      gsub(paste("^.*",term,sep=""),'',grep(term,names(unlist(groups)),value=T)),
      '.',
      fixed=T
      ), function(x) x[2:(distance+1)]
    )))
  )
  if (is.null(val)) {
    NA
  } else {
    val
  }
}

add_child_counts = function(term_counts,groups) {
	unlist(lapply(names(term_counts), function(term) {  term_counts[term] + sum( term_counts[child_terms(term,groups,10)], na.rm=T ) }))
}

overlapping_domains = function(query_interpro,target_interpro,all_domains) {
  domains = all_domains[all_domains$dom %in% c(query_interpro,target_interpro),]
  plyr::ddply(domains,c('uniprot'),function(doms) {
    if ( length(intersect(doms$dom,query_interpro)) < 1 ) {
      return ()
    }
    target_doms = doms[doms$dom %in% target_interpro,]
    target_ranges = IRanges::IRanges(as.numeric(target_doms$start), as.numeric(target_doms$end),names=target_doms$dom)
    query_doms = doms[doms$dom %in% query_interpro,]
    allranges =  IRanges::IRanges(as.numeric(query_doms$start), as.numeric(query_doms$end),names=query_doms$dom)
    over_func = get('%over%',getNamespace('IRanges'))
    overlaps = IRanges::mergeByOverlaps(allranges,target_ranges)
    if (nrow(overlaps) < 1) {
      return ()
    }
    wanted = sapply(1:nrow(overlaps),function(idx){
      if (overlaps[idx,1]@NAMES == overlaps[idx,2]@NAMES) {
        return (FALSE)
      }
      coverage = IRanges::intersect(overlaps[idx,1],overlaps[idx,2])@width
      query_width = overlaps[idx,1]@width
      target_width = overlaps[idx,1]@width
      if (coverage / query_width < 0.5 || coverage / target_width < 0.5) {
        return (FALSE)
      }
      return (TRUE)
    })
    overlaps = overlaps[wanted,]
    if (nrow(overlaps) < 1) {
      return ()
    }
    data.frame(query = sapply(1:nrow(overlaps),function(idx) { overlaps[idx,1]@NAMES }),
               target = sapply(1:nrow(overlaps),function(idx) { overlaps[idx,2]@NAMES }) )
  },.progress='text')
}

expand_classes = function(classes,groups) {
  plyr::ddply(classes,'interpro',function(df) {
    children = child_terms(df$interpro[1],all_groups,10)
    rbind(df,expand.grid(interpro=children[!is.na(children)],Class=df$Class))
  })
}

apply_classes = function(domains,classes) {
  classlist = classes$Class
  names(classlist) = classes$interpro
  domains$dom = ifelse(domains$dom %in% names(classlist),classlist[domains$dom],domains$dom)
  return (domains)
}

all_groups = group_interpros(interpro_raw_relations)

if ( ! exists('domainsets')) {
  domainsets = Rgator::calculateDomainSets(all_sites.human,'site',human.domains,stem_distance=50,remove_tm_overlaps = TRUE)
}

raw_counts = table(unique(rbind(domainsets$inside[,c('uniprot','dom')],domainsets$between[,c('uniprot','dom')]))$dom)
child_summed_counts = sort(add_child_counts(raw_counts,all_groups))

high_count_terms = child_summed_counts[child_summed_counts >= 10]

low_count_terms = child_summed_counts[child_summed_counts < 10 & ! names(child_summed_counts) %in% unlist(sapply(names(high_count_terms),function(x) child_terms(x,all_groups,10))) ]

all_classes = get_classes(c(names(high_count_terms),names(low_count_terms[low_count_terms >= 2])))

overlaps = overlapping_domains(all_classes[is.na(all_classes$Class),'interpro'], all_classes[! is.na(all_classes$Class),'interpro'], human.domains[human.domains$uniprot %in% all_sites.human$uniprot,])
overlaps$key = paste(overlaps$query,overlaps$target)
overlap_counts = table(unique(overlaps)$key)

overlap_mapping = unique(overlaps[ ! grepl("TMhelix",overlaps$key) & overlaps$key %in% names(overlap_counts[overlap_counts > 2]),c('query','target')])

overlap_classes = setNames(merge(overlap_mapping,all_classes,by.x='target',by.y='interpro')[,c('query','Class')],c('interpro','Class'))

all_classes = get_classes(unique(domainsets$all$dom))

all_classes = rbind(all_classes[! is.na(all_classes$Class),],overlap_classes)

missing_doms = unique(c( domainsets$inside[! domainsets$inside$dom %in% all_classes$interpro,]$dom, domainsets$between[! domainsets$between$dom %in% all_classes$interpro,]$dom))

missing_doms = missing_doms[ ! missing_doms %in% unlist(sapply(names(high_count_terms),function(x) child_terms(x,all_groups,10))) ]

missing_low_frequency = sort(base::table(unique(human.domains[human.domains$uniprot %in% all_sites.human$uniprot & human.domains$dom %in% missing_doms,c('uniprot','dom')])$dom))

final_classes = expand_classes(all_classes)

human.glycodomains = apply_classes(human.domains,final_classes)

#cytosolic = read.delim('cytosolic.txt',header=F)$V1

domainsets.glycodomain = Rgator::calculateDomainSets(all_sites.human,'site',human.glycodomains,stem_distance=50,remove_tm_overlaps = TRUE)

#sort(table(unique(domainsets.glycodomain$inside[,c('uniprot','dom')])$dom))