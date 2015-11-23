
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
	aas.melted$seqname = gdata::reorder.factor(aas.melted$seqname,new.order=rev(unique(seqnames)))
	# aas.melted = rbind(aas.melted,data.frame(seqname="conservation",seq.id="conservation",start=-1,end=-1,pos=-1,aa='X'))
	# +coord_fixed(2)+
	#aas.melted$seqname = factor(aas.melted$seqname,levels=unique( unique(aas.melted$seqname)))
	out.plot = ggplot(aas.melted,aes(seqname=seqname,x=pos,y=seqname,pos=pos,start=start,end=end,seq.id=seq.id,aa=aa))+scale_x_discrete(limit=1:max(aas.melted$pos),breaks=seq(from=0,to=max(aas.melted$pos),by=50))+theme_minimal()+theme(panel.grid=ggplot2::element_blank())
	out.plot + scale_y_discrete()
}

get_alignment = function(df) {
	df$identifier = factor(df$identifier,levels=c('HSV-1','HSV-2','VZV','HCMV','EBV'))
	df = df[order(df$identifier),]
	seqs = Rgator::getUniprotSequences(df$uniprot_id)
	seqs$uniprot = toupper(seqs$uniprot)
	seqs = merge(df,seqs,by.x='uniprot_id',by.y='uniprot')$sequence
	seqnames = paste(df$protein,df$identifier,sep='.')
	seqs = lapply( 1:length(seqs), function(seq_id) { 
		views = Biostrings::Views(seqs[seq_id],start=1,end=nchar(seqs[seq_id]))
		names(views) = seqnames[seq_id]
		views
	})
	names(seqs) = seqnames
	aln = msa::msaClustalOmega( Reduce(c,Map(function(view) { as(view,'AAStringSet')}, seqs)) )
	attributes(aln)$sequences = seqs
	aln
}

plot_ieva = function() {
	annotations = read.csv('ieva annotations 150821v2.csv',header=T)
	annotations = annotations[order(annotations$protein),]
	proteins = unique(annotations$protein)
	annotations$seqname = paste(annotations$protein,annotations$identifier,sep='.')
	sites = annotations[annotations$annotation=='UNAMBIGUOUS',c('seqname','annotation','start')]
	sites$class = 'galnac'
	sigpeps = annotations[annotations$annotation=='SP',c('seqname','annotation','start','end')]
	tms = annotations[annotations$annotation=='TM',c('seqname','annotation','start','end')]
	ambiguous = unique(annotations[annotations$annotation=='AMBIGUOUS',c('seqname','annotation','start','end','number_of_ambiguous')])

	plots = lapply(proteins,function(prot) {
		plot_alignment(get_alignment(unique(annotations[annotations$protein == prot,c('protein','uniprot_id','identifier')])))+
		geom_barcode(overlay=F)+
		geom_segment(aes(x=..seqstart..,xend=..seqend..),stat="gappedSequence",size=2,colour="black")+
		geom_sugar(aes(x=..start..),stat="alignedSite",annotations=sites,columns=c('start'),offset=0)+
		geom_segment(aes(x=..start..,xend=..end..),stat="alignedSite",colour="red",size=4,alpha=0.5,annotations=sigpeps,columns=c('start','end'))+
		geom_segment(aes(x=..start..,xend=..end..),stat="alignedSite",colour="blue",size=4,alpha=0.5,annotations=tms,columns=c('start','end'))+
		geom_segment(aes(x=..start..,xend=..end..),stat="alignedSite",colour="yellow",size=4,alpha=0.5,annotations=ambiguous,columns=c('start','end'))+
		geom_text(aes(x=..start..,label=..number_of_ambiguous..),stat="alignedSite",colour="black",nudge_y=-0.25,size=5,alpha=1,annotations=ambiguous,columns=c('start','end'))
	})
	do.call(gridExtra::arrangeGrob,c( plots,ncol=1))
}

get_plot = function() {

	site_data = data.frame(seqname=c('2','1'),site=c(4,5),class=c('galnac','gal(b1-3)galnac'))

	sigpep_data = data.frame(seqname=c('2','1'),start=c(1,1),end=c(3,4))

	plot = 	plot_alignment(do_alignment(c('MNTTTMMMNPPPP','NNSMMMPP')))+
			geom_barcode(overlay=F)+
			geom_segment(aes(x=..seqstart..,xend=..seqend..),stat="gappedSequence",size=2,colour="black")+
			geom_sugar(aes(x=..site..),stat="alignedSite",annotations=site_data,columns=c('site'),offset=0)+
			geom_segment(aes(x=..start..,xend=..end..),stat="alignedSite",colour="red",size=4,alpha=0.5,annotations=sigpep_data,columns=c('start','end'))


#			geom_text(aes(x=pos,y=seqname,label=aa),size=3)

	plot
}
