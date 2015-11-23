#' Draw sugars
#' @export
geom_sugar <- function(mapping = NULL, data = NULL, stat = "identity",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T,offset=2,...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSugar,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm=na.rm,
      offset=offset,
      ...
    )
  )
}

draw_sugar = function(x,y,sugar,offset) {
  if (sugar == "galnac") {
    return (draw_galnac(x,y,offset))
  }
  if (sugar == "gal(b1-3)galnac") {
    return (draw_galb13galnac(x,y,offset))
  }
}

draw_galnac = function(x,y,offset=0) {
  grid::rectGrob(
    x, grid::unit(y,"native") + grid::unit(offset * .pt,"mm"),
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

draw_galb13galnac = function(x,y,offset=0) {
  galnac = draw_galnac(x,y,offset)
  gal = grid::circleGrob(
    x, grid::unit(y,"native") + grid::unit((offset+1) * .pt,"mm"),
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
                        draw_panel = function(data, panel_scales, coord,offset=2) {
                          coords <- coord$transform(data, panel_scales)
                          draw_sugar_vec = Vectorize(draw_sugar,SIMPLIFY=F)
                          results = draw_sugar_vec(coords$x,coords$y,coords$class,offset)
                          do.call(grid::gList,results)
                        }
)
