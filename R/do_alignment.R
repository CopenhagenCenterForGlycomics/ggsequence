#' @importFrom methods as
do_alignment = function(sequences) {
	if (is.null(names(sequences))) {
		names(sequences) = as.character(1:length(sequences))
	}
	seqs = lapply( 1:length(sequences), function(seq_id) { 
		views = Biostrings::Views(sequences[seq_id],start=1,end=nchar(sequences[seq_id]))
		names(views) = names(sequences)[seq_id]
		views
	})
	names(seqs) = names(sequences)
	aln = msa::msaClustalOmega( Reduce(c,Map(function(view) { methods::as(view,'AAStringSet')}, seqs)) )
	attributes(aln)$sequences = seqs
	aln
}


do_alignment.ranges = function(ranges=data.frame(uniprot=c('O00533','P12345'),start=c(1,20),end=c(23,25))) {
	seqranges = ranges %>%
	dplyr::mutate(uniprot = tolower(uniprot)) %>%	
	merge( Rgator::getUniprotSequences(.$uniprot)) %>%
	dplyr::mutate(uniprot = toupper(uniprot)) %>%
			dplyr::group_by(uniprot) %>%
			dplyr::group_map( ~ {
				views = Biostrings::Views(.$sequence[1],start=.$start,end=.$end)
				names(views) = paste(.$uniprot,.$start,.$end,sep='.')
				views
			},.keep=T)

	aln = msa::msaClustalOmega( Reduce(c,Map(function(view) { methods::as(view,'AAStringSet')}, seqranges)) )
	names(seqranges) = toupper(unique(ranges$uniprot))
	attributes(aln)$sequences = seqranges
}