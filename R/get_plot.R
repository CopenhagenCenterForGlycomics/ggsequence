get_plot = function() {

	site_data = data.frame(seq.ids=c('bar','foo'),site=c(4,5),sugar=c('galnac','gal(b1-3)galnac'))

	sigpep_data = data.frame(seq.ids=c('bar','foo'),start=c(1,1),end=c(3,4))

	plot = 	ggplot(do_alignment(c(foo='MNTTTMMMNPPPPMNTTTMMMNPPPPMNTTTMMMNPPPPMNTTTMMMNPPPP',bar='NNSMMMPPNNSMMMPPNNSMMMPPNNSMMMPP')))+
			geom_barcode(overlay=F)+
	    geom_text(aes(label=..aa..,color=..conservation..),stat="conservation")+
			geom_segment(aes(x=..seqstart..,xend=..seqend..),stat="gappedSequence",size=2,colour="black")+
			ggsugar::geom_sugar(aes(x=..site..),stat="alignedSite",annotations=site_data,id.column='seq.ids',columns=c('site'),offset=0)+
			geom_bracket(aes(x=..start..,xend=..end..),offset=1,size=0.5,stat="alignedSite",annotations=sigpep_data,id.column='seq.ids',columns=c('start','end'))+
			geom_segment(aes(x=..start..,xend=..end..),stat="alignedSite",colour="red",size=4,alpha=0.5,annotations=sigpep_data,id.column='seq.ids',columns=c('start','end'))
	plot
}