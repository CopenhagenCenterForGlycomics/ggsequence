plot_sequence <- function(...) {
	sequence=list(...)
	if (length(sequence) > 1) {
		stop("More than 1 sequence needs to be aligned")
	}
	seqnames = names(sequence)
	seq.ids = names(sequence)
	sequence = sequence[[1]]
	starts = BiocGenerics::start(sequence)[1]
	ends = BiocGenerics::end(sequence)[1]
	seqdata = data.frame(seq.id=seq.ids,seqname=seqnames,start=starts,end=ends)
	aas.frame = as.data.frame(matrix(unlist(strsplit(as.character(sequence),NULL)),nrow=1))
	colnames(aas.frame) <- 1:ncol(aas.frame)
	aas.frame$seqname = seqnames
	aas.melted = reshape2::melt( merge(aas.frame,seqdata,by='seqname'), c('seqname','seq.id','start','end'),variable.name='pos',value.name='aa')
	aas.melted$pos = as.numeric(aas.melted$pos)
	aas.melted$seqname = factor(aas.melted$seqname,levels=rev(unique(seqnames)))
	out.plot = 	ggplot2::ggplot(aas.melted,ggplot2::aes(seqname=seqname,x=pos,y=seqname,pos=pos,start=start,end=end,seq.id=seq.id,aa=aa))+
				ggplot2::scale_x_discrete(name="Amino acid position",limit=1:max(aas.melted$pos),breaks=seq(from=0,to=max(aas.melted$pos),by=10))+
				ggplot2::theme_minimal()+
				ggplot2::theme(panel.grid=ggplot2::element_blank())
	out.plot = out.plot + ggplot2::scale_y_discrete(name="Sequence")
	out.plot
}

#' Constructor for a plot from an XString object
#' Give the sequence a name when calling ggplot(my.sequence=AAString('STAT'))
#' @importFrom ggplot2 ggplot
#' @export
ggplot.XString = function(...) {
	plot_sequence(...)
}