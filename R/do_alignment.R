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