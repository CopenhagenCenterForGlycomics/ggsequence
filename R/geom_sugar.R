alignment_default_aes <- function(mapping) {
  default_mapping = ggplot2::aes(seqname=seqname,x=pos,y=seqname,pos=pos,start=start,end=end,seq.id=seq.id,aa=aa )
  if (is.null(mapping)) {
    mapping = default_mapping
  }
  default_mapping = default_mapping[ ! names(default_mapping) %in% names(mapping), ]
  mapping = do.call( ggplot2::aes, c( mapping, default_mapping ) )
}

#' Draw sugars
#' @export
geom_sugar <- function(mapping = NULL, data = NULL, stat = "alignedSite",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T,...) {
  ggplot2::layer(
    data = data,
    mapping = alignment_default_aes(mapping),
    stat = stat,
    geom = GeomSugar,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm=na.rm,
      ...
    )
  )
}

draw_sugar = function(x,y,sugar) {
  if (sugar == "galnac") {
    return (draw_galnac(x,y))
  }
  if (sugar == "gal(b1-3)galnac") {
    return (draw_galb13galnac(x,y))
  }
}

draw_galnac = function(x,y) {
  grid::rectGrob(
    x, grid::unit(y,"native") + grid::unit(2 * .pt,"mm"),
    width = grid::unit(1 * .pt ,"mm"),
    height = grid::unit(1 * .pt,"mm"),
    default.units = "native",
    just = c("centre", "centre"),
    gp = grid::gpar(
      col = ggplot2::alpha("black", 1),
      fill = ggplot2::alpha("yellow", 1),
      lwd = 0.5 * .pt,
      lty = 1,
      lineend = "butt"
    )
  )
}

draw_galb13galnac = function(x,y) {
  galnac = draw_galnac(x,y)
  gal = grid::circleGrob(
    x, grid::unit(y,"native") + grid::unit(3 * .pt,"mm"),
    r = grid::unit(0.5 * .pt ,"mm"),
    default.units = "native",
    gp = grid::gpar(
      col = ggplot2::alpha("black", 1),
      fill = ggplot2::alpha("yellow", 1),
      lwd = 0.5 * .pt,
      lty = 1,
      lineend = "butt"
    )
  )
  grid::gList(gal,galnac)
}

.pt <- 72.27 / 25.4

#' @export
GeomSugar <- ggplot2::ggproto("GeomSugar", ggplot2::Geom,
                        required_aes='x',
                        draw_panel = function(data, panel_scales, coord) {
                          coords <- coord$transform(data, panel_scales)
                          draw_sugar_vec = Vectorize(draw_sugar)
                          results = draw_sugar_vec(coords$x,coords$y,coords$class)
                          do.call(grid::gList,results)
                        }
)
