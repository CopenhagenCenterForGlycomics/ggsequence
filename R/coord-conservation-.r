require(grid)

#' Generate barcodes and amino acid indexes for aligned proteins
#' @param alignment_axis Show the amino acids positions for the individual aligned sequences
#' @param conservation_axis Show the conservation barcode for the aligned sequences
#' @export
coord_conservation <- function(alignment_axis=TRUE,conservation_axis=TRUE) {
  ggplot2::ggproto(NULL, CoordConservation,
    limits = list(x = NULL, y = NULL),
    expand = TRUE,
    default = FALSE,
    clip = "on",
    alignment_axis = alignment_axis,
    conservation_axis = conservation_axis
  )
}

draw_conservation_grobs <- function(conservation) {
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
  conservation = conservation[conservation$x >= 0,]
  grid::rectGrob(
    conservation$x - 0.5*conservation$width, 0,
    width = conservation$width,
    height = grid::unit(2*.pt,"mm"),
    default.units = "native",
    just = c("left", "top"),
    gp = grid::gpar(
      fill = ggplot2::alpha("black",0.1+conservation$value*0.9),
      col='black',
      lwd = 0.3,
      lineend = "butt"
    )
  )
}

get_conservation_from_scale <- function(scale_obj) {
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

get_aa_indexes_from_scale <- function(scale_obj) {
  indexes = attributes(scale_obj$scale)$aligned_indexes
  limits=c( min(scale_obj$limits), max(scale_obj$limits))
  rescaler=scale_obj$rescale
  breaks = scale_obj$get_breaks()
  breaks = breaks[! is.na(breaks) ]
  break_size = 10
  if (length(breaks) > 1) {
    break_size = breaks[2] - breaks[1]
  }
  indexes_breaks = lapply(indexes, function(seq_idxs) {
    seq_idxs = stats::setNames(seq_idxs,1:length(seq_idxs))
    seq_idxs = seq_idxs[seq_idxs >= limits[1] & seq_idxs <= limits[2]]
    if (length(seq_idxs) < 1) {
      return (data.frame())
    }
    wanted_indexes = c(1,seq(0,max(as.numeric(names(seq_idxs))),by=break_size))
    wanted_positions = rescaler(seq_idxs[wanted_indexes[wanted_indexes > 0]])
    wanted_positions = wanted_positions[!is.na(wanted_positions)]
    data.frame(x=wanted_positions,label=names(wanted_positions))
  })
  indexes_breaks
}

draw_axis_labels <- function(indexes,element) {
  lapply(1:length(indexes), function(offset) {
    df = indexes[[offset]]
    df = df[df$x >= 0,]
    if (nrow(df) < 1) {
      return(ggplot2::zeroGrob())
    }
    rlang::exec(ggplot2::element_grob,element,
      df$x,
      TRUE,
      label = df$label,
      check.overlap = FALSE)
  })
}

#' @export
CoordConservation <- ggproto("CoordConservation", CoordCartesian,

  is_linear = function() TRUE,
  is_free = function() TRUE,
  render_axis_h = function(panel_params, theme) {

    bottom_axis_grob = panel_guides_grob(panel_params$guides, position = "bottom", theme = theme)

    margin = unit(0.5, "line")
    if (panel_params$alignment_axis) {
      new_axis = rev(draw_axis_labels(get_aa_indexes_from_scale(panel_params$x),calc_element('axis.text.x',theme)))
      new_heights = do.call(grid::unit.c, lapply(new_axis,grid::grobHeight))
      new_heights = new_heights
      all_heights = grid::unit.c( grid::grobHeight(bottom_axis_grob), new_heights)+2*margin
      bottom_axis_grob <- gridExtra::arrangeGrob(grobs = c(list(bottom_axis_grob), new_axis ),
        heights=all_heights,ncol=1
      )
    }
    top_axis_grob = zeroGrob()
    if (panel_params$conservation_axis) {
      top_axis_grob <- draw_conservation_grobs(get_conservation_from_scale(panel_params$x))
    }
    list(
      top = top_axis_grob,
      bottom = bottom_axis_grob
    )
  },

  setup_panel_params = function(self,scale_x,scale_y,params=list()) {
    if ('conservation' %in% names(params)) {
      attributes(scale_x)$conservation <- params$conservation
    }
    if ('aligned_indexes' %in% names(params)) {
      attributes(scale_x)$aligned_indexes <- params$aligned_indexes      
    }
    parent <- ggproto_parent(CoordCartesian, self)
    panel_params <- parent$setup_panel_params(scale_x, scale_y, params)
    panel_params$conservation_axis = self$conservation_axis
    panel_params$alignment_axis = self$alignment_axis
    panel_params
  },
  
  setup_params = function(data) {
    alignment_data = data[[1]]
    rescaled_idxs = lapply( unique(alignment_data$seq.id), function(seqid) {
      max_idx = max(alignment_data[alignment_data$seq.id == seqid, 'end'])
      min_idx = max(alignment_data[alignment_data$seq.id == seqid, 'start'])
      fullseq = paste(alignment_data[alignment_data$seq.id == seqid, 'aa'],collapse='')
      sapply( min_idx:max_idx, function(site) {
        rescale_site( fullseq , site )
      })
    })
    list(conservation=get_conservation(alignment_data),aligned_indexes=rescaled_idxs)
  }

)

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

panel_guides_grob <- function(guides, position, theme) {
  guide <- guide_for_position(guides, position) %||% ggplot2::guide_none()
  ggplot2::guide_gengrob(guide, theme)
}

guide_for_position <- function(guides, position) {
  has_position <- vapply(
    guides,
    function(guide) identical(guide$position, position),
    logical(1)
  )

  guides <- guides[has_position]
  guides_order <- vapply(guides, function(guide) as.numeric(guide$order)[1], numeric(1))
  Reduce(ggplot2::guide_merge, guides[order(guides_order)])
}