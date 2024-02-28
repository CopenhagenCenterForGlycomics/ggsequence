#' @importFrom methods is
plot_alignment <- function(data,mapping=ggplot2::aes(),...) {
  stopifnot('sequences' %in% names(attributes(data)))
  starts = unlist(as.vector(sapply(attributes(data)$sequences,function(x) BiocGenerics::start(x)  )))
  ends = unlist(as.vector(sapply(attributes(data)$sequences,function(x) BiocGenerics::end(x) )))
  seqnames = unlist(as.vector(sapply(names(attributes(data)$sequences),function(seq.id) {
    if ( methods::is( attributes(data)$sequences[[seq.id]], 'XString' )) {
      seq.id
    } else {
      names( attributes(data)$sequences[[seq.id]] )
    }
  })))
  seq.ids = unlist(sapply(names(attributes(data)$sequences),function(seq.id) {
    if ( methods::is( attributes(data)$sequences[[seq.id]], 'XStringViews' ) |
       methods::is( attributes(data)$sequences[[seq.id]], 'ViewSet' )
       ) {
      rep(seq.id,length(attributes(data)$sequences[[seq.id]]))
    } else {
      seq.id
    }
  },simplify=F) )
  seqdata = data.frame(seq.id=seq.ids,seqname=seqnames,start=starts,end=ends)
  aas.frame = as.data.frame(Biostrings::as.matrix(data))
  colnames(aas.frame) <- 1:ncol(aas.frame)
  aas.melted = reshape2::melt( cbind(aas.frame,seqdata), c('seqname','seq.id','start','end'),variable.name='pos',value.name='aa')
  aas.melted$pos = as.numeric(aas.melted$pos)
  aas.melted$seqname = factor(aas.melted$seqname,levels=rev(unique(seqnames)))

  desired_aes = ggplot2::aes_string(seqname="seqname",x="pos",y="seqname",pos="pos",start="start",end="end",seq.id="seq.id",aa="aa")
  desired_aes = c(desired_aes,mapping)
  class(desired_aes) <- "uneval"

  out.plot =  ggplot2::ggplot(aas.melted,desired_aes,...)+
        ggplot2::scale_x_discrete(name="Amino acid position",limit=1:max(aas.melted$pos),breaks=seq(from=0,to=max(aas.melted$pos),by=10))+
        ggplot2::theme_minimal()+
        ggplot2::theme(panel.grid=ggplot2::element_blank())
  out.plot = out.plot + ggplot2::scale_y_discrete(name="Sequence")
  class(out.plot) = c('aligned.plot',class(out.plot))
  out.plot
}

#' Constructor for a plot from an alignment
#' @importFrom ggplot2 ggplot
#' @export
ggplot.AAMultipleAlignment = function(data,mapping=ggplot2::aes(),...,environmemt=parent.frame()) {
  plot_alignment(data,mapping,...)
}