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

  consframe = plyr::dlply(data,'pos',function(aas) {
    if ( length(unique(aas$aa)) == 1 && aas$aa[1] != '-' ) {
      1
    } else if ( max(sapply(strong_res,function(re) { length(grep(re,aas$aa)) })) == length(aas$aa) ) {
      2/3
    } else if ( max(sapply(weak_res,function(re) { length(grep(re,aas$aa)) })) == length(aas$aa) ) {
      1/3
    } else {
      0
    }
  })
  unlist(consframe)
}


conservation_default_aes <- function(mapping) {
  default_mapping = ggplot2::aes_string(alpha="stat(conservation)")
  if (is.null(mapping)) {
    mapping = default_mapping
  }
  missing_names = ! names(default_mapping) %in% names(mapping)
  if (any(missing_names)) {
    default_mapping = default_mapping[missing_names]
  }
  mapping = do.call( ggplot2::aes, c( mapping, default_mapping ) )
}

#' Wrap the venn diagram into a geom
#' @export
geom_barcode <- function(mapping = NULL, data = NULL, stat = "conservation",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T,height=2,...) {
  layer = ggplot2::layer(
    data = data,
    mapping = conservation_default_aes(mapping),
    stat = StatConservation,
    geom = GeomBarcode,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm=na.rm,
      height=height,
      ...
    )
  )

  layer
}

is.zero <- function(x) is.null(x) || inherits(x, "zeroGrob")

#' @export
GeomBarcode <- ggplot2::ggproto("GeomBarcode", ggplot2::GeomTile,
                        draw_panel = function(self,data, panel_scales, coord,height=10) {
                          data = unique(data[,c('xmin','xmax','ymin','ymax','alpha','fill','size','colour','linetype')])
                          data$ymin = min(data$ymin)
                          data$ymax = max(data$ymax)
                          return(self$super()$draw_panel(unique(data), panel_scales, coord))
                        }
)

#' Stat to calculate a conservation value using the clustalw scoring
#' @export
stat_conservation <- function(mapping = NULL, data = NULL, geom = "barcode",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T,alone=T,...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = StatConservation,
    geom = geom,
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
                          data$conservation = unlist(sapply(get_conservation(data), function(x) rep(x,max(data$y)-min(data$y)+1),simplify=F))
                          data$pos = as.numeric(data$pos)
                          if (alone) {
                            data$y = as.numeric(data$seqname)
                          }
                          return(unique(data))
                        }
)