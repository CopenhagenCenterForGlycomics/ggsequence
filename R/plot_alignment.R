
plot_alignment <- function(alignment) {
	starts = unlist(as.vector(sapply(attributes(alignment)$sequences,BiocGenerics::start)))
	ends = unlist(as.vector(sapply(attributes(alignment)$sequences,BiocGenerics::end)))
	seqnames = unlist(as.vector(sapply(attributes(alignment)$sequences,names)))
	seq.ids = unlist(sapply(names(attributes(alignment)$sequences),function(seq.id) { rep(seq.id,length(attributes(alignment)$sequences[[seq.id]])) },simplify=F) )
	seqdata = data.frame(seq.id=seq.ids,seqname=seqnames,start=starts,end=ends)
	aas.frame = as.data.frame(Biostrings::as.matrix(alignment))
	colnames(aas.frame) <- 1:ncol(aas.frame)
	aas.frame$seqname = rownames(aas.frame)
	aas.melted = reshape2::melt( merge(aas.frame,seqdata,by='seqname'), c('seqname','seq.id','start','end'),variable.name='pos',value.name='aa')
	aas.melted$pos = as.numeric(aas.melted$pos)
	aas.melted$seqname = factor(aas.melted$seqname,levels=rev(unique(seqnames)))
	out.plot = 	ggplot(aas.melted,aes(seqname=seqname,x=pos,y=seqname,pos=pos,start=start,end=end,seq.id=seq.id,aa=aa))+
				scale_x_discrete(name="Amino acid position",limit=1:max(aas.melted$pos),breaks=seq(from=0,to=max(aas.melted$pos),by=10))+
				theme_minimal()+
				theme(panel.grid=ggplot2::element_blank())
	out.plot = out.plot + scale_y_discrete(name="Sequence")
	class(out.plot) = c('aligned.plot',class(out.plot))
	out.plot
}

#' Constructor for a plot from an alignment
#' @export
ggplot.MsaAAMultipleAlignment = function(alignment) {
	plot_alignment(alignment)
}