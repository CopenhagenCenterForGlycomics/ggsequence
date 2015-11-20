get_conservation <- function(data) {

  strong_res = c( "[STA]",
                  "[NEQK]",
                  "[NHQK]",
                  "[NDEQ]",
                  "[QHRK]",
                  "[MILV]",
                  "[MILF]",
                  "[HY]",
                  "[FYW]")
 
  weak_res = c( "[CSA]",
                "[ATV]",
                "[SAG]",
                "[STNK]",
                "[STPA]",
                "[SGND]",
                "[SNDEQK]",
                "[NDEQHK]",
                "[NEQHRK]",
                "[FVLIM]",
                "[HFY]")

  consvals = plyr::ddply(data,'pos',function(aas) {
    if ( length(unique(aas$aa)) == 1 && aas$aa[1] != '-' ) {
      aas$cons <- 1
    } else if ( max(sapply(strong_res,function(re) { length(grep(re,aas$aa)) })) == length(aas$aa) ) {
      aas$cons <- 2/3
    } else if ( max(sapply(weak_res,function(re) { length(grep(re,aas$aa)) })) == length(aas$aa) ) {
      aas$cons <- 1/3
    } else {
      aas$cons <- 0
    }
    aas
  })$cons

  consvals
}


conservation_default_aes <- function(mapping) {
  default_mapping = ggplot2::aes(seqname=seqname,x=pos,y=seqname,pos=pos,start=start,end=end,seq.id=seq.id,aa=aa,alpha=..conservation.. )
  if (is.null(mapping)) {
    mapping = default_mapping
  }
  default_mapping = default_mapping[ ! names(default_mapping) %in% names(mapping), ]
  mapping = do.call( ggplot2::aes, c( mapping, default_mapping ) )
}

#pBrackets

overlay_conservation <- function(alignment.plot) {
  alignment.plot + geom_tile(aes(x=pos,y=seqid),fill='black',alpha=get_conservation(alignment.plot$data))
}

#' Wrap the venn diagram into a geom
#' @export
geom_barcode <- function(mapping = NULL, data = NULL, stat = "conservation",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T,overlay=F,...) {
  ggplot2::layer(
    data = data,
    mapping = conservation_default_aes(mapping),
    stat = StatConservation,
    geom = GeomBarcode,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm=na.rm,
      overlay=overlay,
      ...
    )
  )
}

#' @export
GeomBarcode <- ggplot2::ggproto("GeomBarcode", ggplot2::GeomTile,
                        draw_panel = function(self,data, panel_scales, coord,overlay=F) {
                          if ( ! overlay ) {
                            data$ymin = 0.95
                            data$ymax = 1
                            panel_scales$y.range = c( panel_scales$y.range[1], 1)
                          }
                          data = unique(data[,c('xmin','xmax','ymin','ymax','alpha','fill','size','colour','linetype')])
                          return (self$super$draw_panel(unique(data[,c('xmin','xmax','ymin','ymax','alpha','fill','size','colour','linetype')]),panel_scales,coord))
                        }
)


#' Stat to calculate a conservation value using the clustalw scoring
#' @export
stat_conservation <- function(mapping = NULL, data = NULL, geom = "barcode",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T,alone=T,...) {
  ggplot2::layer(
    data = data,
    mapping = conservation_default_aes(mapping),
    stat = StatConservation,
    geom = "barcode",
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm=na.rm,
      fill="black",
      alone=alone,
      ...
    )
  )
}

#' @export
StatConservation <- ggplot2::ggproto("StatConservation", ggplot2::Stat,
                        default_aes = ggplot2::aes(),
                        setup_data = function(data,params) {
                          data[data$seqname != 'conservation',]
                        },
                        compute_panel = function(data,scales,alone=T) {
                          data$conservation = get_conservation(data)
                          data$pos = as.numeric(data$pos)
                          if (alone) {
                            data$seqname = data$seqname #factor(rep("conservation",nrow(data)),levels=levels(data$seqname) )
                            data$y = as.numeric(data$seqname)
                          }
                          return(unique(data))
                        }
)