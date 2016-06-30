if ( ! exists('interpro_raw_relations') ) {
  download.file('ftp://ftp.ebi.ac.uk/pub/databases/interpro/current/ParentChildTreeFile.txt','interpro_relations.txt')
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

get_classes = function(interpro) {
  go_mapped = merge(interpro_raw_go[interpro_raw_go$interpro %in% interpro & interpro_raw_go$go %in% go_manual_mappings$go,], go_manual_mappings,by='go')[,c('interpro','Class')]
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

expand_classes = function(classes,groups) {
  plyr::ddply(classes,'interpro',function(df) {
    children = child_terms(df$interpro[1],all_groups,10)
    child_classes = get_classes(children[!is.na(children)])
    defined_classes = child_classes[!is.na(child_classes$Class),]
    undefined_classes = child_classes[is.na(child_classes$Class),]
    rbind(df,defined_classes,expand.grid(interpro=undefined_classes$interpro,Class=df$Class))
  },.progress='text')
}

all_groups = group_interpros(interpro_raw_relations)

# We should write the overlapping domains out somewhere...

overlap_mapping = read.delim('overlap_manual.tsv',header=T)

overlap_classes = setNames(merge(overlap_mapping,get_classes(overlap_mapping$target),by.x='target',by.y='interpro')[,c('source','Class')],c('interpro','Class'))

all_classes = get_classes(unique(c(interpro_raw_go$interpro,interpro_manual_mappings$interpro)))

all_classes = rbind(all_classes[! is.na(all_classes$Class),],overlap_classes)

final_classes = expand_classes(all_classes)

# Export final_classes to a file that we upload somewhere.

write.table(unique(final_classes[final_classes$Class != 'Skipped',]),file='Glycodomain-class.tsv',row.names=F,col.names=T,quote=F)