if ( ! exists('interpro_raw_relations') ) {
  download.file('ftp://ftp.ebi.ac.uk/pub/databases/interpro/ParentChildTreeFile.txt','interpro_relations.txt')
  interpro_raw_relations = gsub('::.*','',readLines('interpro_relations.txt'))
}

if ( ! exists('interpro_go_mappings')) {
  download.file('http://www.geneontology.org/external2go/interpro2go','interpro_go.txt')
  interpro_raw_go = gsub('^InterPro:','',Filter(function(x) grepl("^InterPro",x),readLines('interpro_go.txt')))
  interpro_raw_go = as.data.frame(matrix( unlist(lapply(interpro_raw_go,function(x) { lapply(strsplit(x,' '), function(vals) { vals[c(1,length(vals))] }) })), ncol=2,byrow=T))
  interpro_raw_go = setNames(interpro_raw_go,c('interpro','go'))
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


all_groups = group_interpros(interpro_raw_relations)

domainsets = Rgator::calculateDomainSets(all_sites.human,'site',human.domains,stem_distance=50,remove_tm_overlaps = TRUE)
raw_counts = table(unique(rbind(domainsets$inside[,c('uniprot','dom')],domainsets$between[,c('uniprot','dom')]))$dom)
child_summed_counts = sort(add_child_counts(raw_counts,all_groups))

high_count_terms = child_summed_counts[child_summed_counts >= 10]

low_count_terms = child_summed_counts[child_summed_counts < 10 & ! names(child_summed_counts) %in% unlist(sapply(names(high_count_terms),function(x) child_terms(x,all_groups,10))) ]
