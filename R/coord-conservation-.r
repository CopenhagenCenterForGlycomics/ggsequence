require(grid)

coord_conservation <- function(xlim = NULL, ylim = NULL, expand = TRUE,
                            default = FALSE, clip = "on") {
  ggproto(NULL, CoordConservation,
    limits = list(x = xlim, y = ylim),
    expand = expand,
    default = default,
    clip = clip
  )
}

pad_conservation <- function(conservation,range) {
  padding=0
  conservation = c(rep(-10,padding),0, conservation , rep(-10,padding))
  indices = seq(-1*padding,length(conservation)-padding-1)
  conservation[indices[which(indices >= ceiling(range[1]) & indices <= floor(range[2]))]+padding+1]
}

select_conservation <- function(conservation,range) {
  conservation = c(rep(-10,padding),0, conservation , rep(-10,padding))
  indices = seq(-1*padding,length(conservation)-padding-1)
  conservation[indices[which(indices >= ceiling(range[1]) & indices <= floor(range[2]))]+padding+1]
}


mytemp_grob <- function(conservation,range,position,theme) {
  # gTree(children=gList(
  # grid::rectGrob(
  #   0, 0,
  #   width = 1,
  #   height = grid::unit(2*.pt,"mm"),
  #   default.units = "native",
  #   just = c("left", "top"),
  #   gp = grid::gpar(
  #     #fill = alpha("black", 0.1+seq(0,0.9,by=0.3)),
  #     fill = alpha("black",0),
  #     col='red',
  #     lwd = 1,
  #     lineend = "butt"
  #   )
  # ), // OTHER GROB HERE)
  
  grid::rectGrob(
    conservation$x - 0.5*conservation$width, 0,
    width = conservation$width,
    height = grid::unit(2*.pt,"mm"),
    default.units = "native",
    just = c("left", "top"),
    gp = grid::gpar(
      fill = alpha("black",0.1+conservation$value*0.9),
      col='black',
      lwd = 0.3,
      lineend = "butt"
    )
  )
}

get_conservation_from_scale <- function(scale_obj,range) {
  cons=attributes(scale_obj$scale)$conservation
  limits=c( min(scale_obj$limits), max(scale_obj$limits))
  rescaler=scale_obj$rescale
  if (limits[1] < 1) {
    limits[1] = 1
  }
  if (limits[2] > length(cons)) {
    limits[2] = length(cons)
  }
  wanted_cons = cons[limits[1]:limits[2]]
  x = sapply( limits[1]:limits[2], rescaler )
  return( data.frame( x = x, value = wanted_cons, width=rescaler(2)-rescaler(1) ))
}


#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
CoordConservation <- ggproto("CoordConservation", CoordCartesian,

  is_linear = function() TRUE,
  is_free = function() TRUE,

  render_axis_h = function(panel_params, theme) {
    bottom_grob = panel_guides_grob(panel_params$guides, position = "bottom", theme = theme)
    list(
      top = mytemp_grob(get_conservation_from_scale(panel_params$x), position = "bottom", theme = theme),#ggplot2::zeroGrob(),#panel_guides_grob(panel_params$guides, position = "top", theme = theme),
      bottom = bottom_grob
    )
  },

  setup_panel_params = function(self,scale_x,scale_y,params=list()) {
    if ('conservation' %in% names(params)) {
      attributes(scale_x)$conservation <- params$conservation
    }
    parent <- ggproto_parent(CoordCartesian, self)
    panel_params <- parent$setup_panel_params(scale_x, scale_y, params)
    panel_params
  },
  
  setup_params = function(data) {
    alignment_data = data[[1]]
    rescaled_idxs = lapply( unique(alignment_data$seq.id), function(seqid) {
      max_idx = max(alignment_data[alignment_data$seq.id == seqid, 'end'])
      min_idx = max(alignment_data[alignment_data$seq.id == seqid, 'start'])
      sapply( min_idx:max_idx, function(site) {
        rescale_site( alignment_data[alignment_data$seq.id == seqid, 'aa'], site )
      })
    })
    list(conservation=get_conservation(alignment_data),aligned_indexes=rescaled_idxs)
  }

)

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

"%|W|%" <- function(a, b) {
  if (!is.waive(a)) a else b
}

panel_guides_grob <- function(guides, position, theme) {
  guide <- guide_for_position(guides, position) %||% guide_none()
  guide_gengrob(guide, theme)
}

guide_for_position <- function(guides, position) {
  has_position <- vapply(
    guides,
    function(guide) identical(guide$position, position),
    logical(1)
  )

  guides <- guides[has_position]
  guides_order <- vapply(guides, function(guide) as.numeric(guide$order)[1], numeric(1))
  Reduce(guide_merge, guides[order(guides_order)])
}