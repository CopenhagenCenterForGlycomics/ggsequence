#' Stat to calculate the starts and ends of sequence in an alignment.
#' Populate the ..seqstart.. and ..seqend.. so that they contain the
#' start and ends of each segment of sequence.
#' @export
stat_gapped_sequence <- function(mapping = NULL, data = NULL, geom = "rect",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T,...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
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
                            groups = matrix(transitions,ncol=2,byrow=T)
                            groups = as.data.frame(groups)
                            names(groups) = c('seqstart','seqend')
                            groups$y = yval
                            groups
                          })
                          trans$seqstart = trans$seqstart - 0.5
                          trans$seqend = trans$seqend + 0.5
                          trans[trans$seqstart < 0,] = 0
                          trans$seqend = ifelse(trans$seqend > max(data$pos),trans$seqend-0.5,trans$seqend)
                          return(trans)
                        }
)
