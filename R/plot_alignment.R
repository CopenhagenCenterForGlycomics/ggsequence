
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
	# aas.melted = rbind(aas.melted,data.frame(seqname="conservation",seq.id="conservation",start=-1,end=-1,pos=-1,aa='X'))
	aas.melted$seqname = factor(aas.melted$seqname,levels=unique( unique(aas.melted$seqname)))
	out.plot = ggplot(aas.melted)+geom_text(aes(x=pos,y=seqname,label=aa),size=3)+coord_fixed(2)+scale_x_discrete(limit=1:max(aas.melted$pos),breaks=seq(from=0,to=max(aas.melted$pos),by=50))+theme_minimal()
	out.plot + scale_y_discrete(labels=function(lab) { unlist(Map(function(lbls) { paste(lbls[1],lbls[2],sep='@')} ,strsplit(lab,'#'))) })
}

get_plot = function() {
	site_data = data.frame(seqname=c('2#1#6','1#1#13'),site=c(4,5),class=c('galnac','gal(b1-3)galnac'))
	plot = 	plot_alignment(do_alignment(c('MNTTTMMMNPPPP','NNSMMM')))+
			geom_barcode(overlay=F)+
			stat_gapped_sequence(aes(x=..seqstart..),geom="segment",size=2,colour="black")+
			geom_sugar(aes(x=..site..),annotations=site_data,columns=c('site'))
	# browser()
	# ggplot2::ggplot_build(plot)
	plot
}
