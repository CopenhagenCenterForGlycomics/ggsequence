get_plot = function() {

	site_data = data.frame(seq.ids=c('bar','foo'),site=c(4,5),sugar=c('galnac','gal(b1-3)[glcnac(b1-6)]galnac'))

	sigpep_data = data.frame(seq.ids=c('bar','foo'),start=c(1,1),end=c(3,4))

	plot = 	ggplot(ggsequence:::do_alignment(c(foo='MNTTTMMMNPPPPMNTTTMMMNPPPPMNTTTMMMNPPPPMNNSMMMPPPP',bar='NNSMMMPPNNSMMMPPNNSMMMPPNNSMMMPP')))+
		    geom_text(aes(label=stat(aa),color=stat(conservation)),stat="conservation")+
			geom_segment(aes(x=stat(seqstart),xend=stat(seqend)),stat="gappedSequence",size=2,colour="black",alpha=0.5)
	
	  
	if (requireNamespace('ggsugar',quietly = TRUE)) {
	  plot = plot +
			ggsugar::geom_sugar(aes(x=stat(site)),stat="alignedSite",annotations=site_data,id.column='seq.ids',columns=c('site'),size=4,offset=2)
	}
	plot = plot +
			geom_bracket(aes(x=stat(start),xend=stat(end)),offset=1,size=0.5,stat="alignedSite",annotations=sigpep_data,id.column='seq.ids',columns=c('start','end'))+
			geom_segment(aes(x=stat(start),xend=stat(end)),stat="alignedSite",colour="red",size=4,alpha=0.5,annotations=sigpep_data,id.column='seq.ids',columns=c('start','end'))+
			coord_conservation()
	plot
}
