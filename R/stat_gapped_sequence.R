alignment_default_aes <- function(mapping) {
  default_mapping = ggplot2::aes(seqname=seqname,x=pos,y=seqname,pos=pos,start=start,end=end,seq.id=seq.id,aa=aa )
  if (is.null(mapping)) {
    mapping = default_mapping
  }
  default_mapping = default_mapping[ ! names(default_mapping) %in% names(mapping), ]
  mapping = do.call( ggplot2::aes, c( mapping, default_mapping ) )
}


#' Stat to re-align columns in a data frame using the alignment information
#' @export
stat_gapped_sequence <- function(mapping = NULL, data = NULL, geom = "rect",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T,...) {
  ggplot2::layer(
    data = data,
    mapping = alignment_default_aes(mapping),
    stat = StatGappedSequence,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm=na.rm,
      ...
    )
  )
}

#' @export
StatGappedSequence <- ggplot2::ggproto("StatGappedSequence", ggplot2::Stat,
                        default_aes = ggplot2::aes(x=..seqstart..,xend=..seqend..,y=..y..,yend=..y..),
                        compute_panel = function(data,scales,annotations=NULL,columns=c()) {
                          alignment =  apply( t(sapply(1:max(data$y), function(y) { data[data$y == y,'aa']})), 1, function(row) paste(row,collapse='') )
                          trans = plyr::ldply(1:max(data$y),function(yval) {
                            aln = alignment[yval]
                            transitions = sort( c( unlist(gregexpr("\\w-",aln,perl=T)), unlist(gregexpr("-\\w",aln,perl=T))+1))
                            transitions = transitions[ transitions > 0 ]
                            if (substr(aln,1,1) != "-") {
                              transitions = c(1,transitions)
                            }
                            if (substr(aln,nchar(aln),nchar(aln)) != "-") {
                              transitions = c(transitions,nchar(aln))
                            }
                            groups = matrix(transitions,ncol=2)
                            groups = as.data.frame(groups)
                            names(groups) = c('seqstart','seqend')
                            groups$y = yval
                            groups                         
                          })
                          browser()
                          return(trans)
                        }
)