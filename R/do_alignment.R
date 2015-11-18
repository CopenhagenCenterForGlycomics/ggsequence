do_alignment = function(sequences) {
	seqs = lapply( 1:length(sequences), function(seq_id) { 
		views = Biostrings::Views(sequences[seq_id],start=1,end=nchar(sequences[seq_id]))
		names(views) = paste(seq_id,1,sapply(sequences[seq_id],nchar),sep='#')
		views
	})
	names(seqs) = as.character(1:length(sequences))
	aln = msa::msaClustalOmega( Reduce(c,Map(function(view) { as(view,'AAStringSet')}, seqs)) )
	attributes(aln)$sequences = seqs
	aln
}