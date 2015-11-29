#' Stat to re-align columns in a data frame using the alignment information
#' @export
stat_aligned_site <- function(mapping = NULL, data = NULL, geom = "point",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T,annotations=NULL,columns=c(),...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatAlignedSite,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm=na.rm,
      annotations=annotations,
      columns=columns,
      ...
    )
  )
}

#' @export
StatAlignedSite <- ggplot2::ggproto("StatAlignedSite", ggplot2::Stat,
                        default_aes = ggplot2::aes(y=..y..,yend=..y..),
                        compute_panel = function(self,data,scales,id.column=NULL,annotations=NULL,columns=c()) {
                          data$pos = as.numeric(data$pos)
                          return(unique(compute_sites(data,annotations,id.column,columns)))
                        }
)


compute_sites <- function(alignment.data,sites,id.column,site.columns) {
  alignment =  apply( t(sapply(1:max(alignment.data$y), function(y) { alignment.data[alignment.data$y == y,'aa']})), 1, function(row) paste(row,collapse='') )
  merged = merge(setNames(alignment.data[,c('start','end','seq.id','y')], c('aln.start','aln.end','seq.id','y')),sites,by.x='seq.id',by.y=id.column)
  wanted.rows = Reduce(function(old,col) { old & ( ((merged[[col]] - merged$aln.start) * (merged$aln.end - merged[[col]])) >= 0 ) },site.columns,rep(TRUE,nrow(merged)))
  merged = unique(merged[which(wanted.rows),])
  # merged$seqid.draw = sapply(merged$seqid,function(id) { paste(id,'.')})
  for (col in site.columns) {
    merged[[col]] = merged[[col]] - merged$aln.start
    merged[[col]] = apply(merged,1,function(row) {
      rescale_site(alignment[as.numeric(row['y'])],as.numeric(row[col])+1)
    })    
  }
  return (merged[,c(1,4:ncol(merged))])
}

rescale_site <- function(seq,site) {
  seq_parts = Map( function(part) {
    if (grepl('[A-Za-z]',part)) {
      part
    } else {
      nchar(part)
    }
  }, as.list(unlist(strsplit(seq,'\\b',perl=T))))
  pos = 1
  for (part in 1:length(seq_parts)) {
    if (is.numeric(seq_parts[[part]])) {
      pos = pos + seq_parts[[part]]
    } else {
      if (site == 1) {
        return (pos)
      } else {
        site = site - 1
        pos = pos + 1
      }
    }
  }
}