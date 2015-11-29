
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
				scale_x_discrete(name="Amino acid position",limit=1:max(aas.melted$pos),breaks=seq(from=0,to=max(aas.melted$pos),by=50))+
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

#' Create a grob from an alignment plot
#' @export
ggplotGrob.aligned.plot = function(x) {
	g = ggplot2::ggplotGrob(x)

	# We should really check whether there is a geom_conservation layer
	# here, and if the overlay is true or false, and then decide
	# what we do with the clipping

	g$layout$clip[g$layout$name == "panel"] = "off"
	g
}

#' Default method for printing an aligned plot
#' @export
print.aligned.plot = function(x) {
  grid::grid.newpage()
  g = ggplotGrob.aligned.plot(x)
  grid::grid.draw(g)
}

get_plot = function() {

	site_data = data.frame(seq.ids=c('bar','foo'),site=c(4,5),class=c('galnac','gal(b1-3)galnac'))

	sigpep_data = data.frame(seq.ids=c('bar','foo'),start=c(1,1),end=c(3,4))

	plot = 	ggplot(do_alignment(c(foo='MNTTTMMMNPPPP',bar='NNSMMMPP')))+
			geom_barcode(overlay=F)+
			geom_segment(aes(x=..seqstart..,xend=..seqend..),stat="gappedSequence",size=2,colour="black")+
			geom_sugar(aes(x=..site..),stat="alignedSite",annotations=site_data,id.column='seq.ids',columns=c('site'),offset=0)+
			geom_segment(aes(x=..start..,xend=..end..),stat="alignedSite",colour="red",size=4,alpha=0.5,annotations=sigpep_data,id.column='seq.ids',columns=c('start','end'))
	plot
}
