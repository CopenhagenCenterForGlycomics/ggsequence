#' Draw brackets
#' @export
geom_bracket <- function(mapping = NULL, data = NULL, stat = "identity",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T,offset=0,...) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomBracket,
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

# Adapted from the pBrackets library

a_cb_brackets <-
function(phi=1, ticks=0.5, type=1){
  if(phi<0.01) phi <-0.01
  n      <- 1000
  d_start <- 1
  d_end   <- 1
  add_s   <- 0
  add_e   <- 0
  if(!is.null(ticks)){
    ticks  <- unique(ticks)
    tsigns <- sign(ticks)
    o      <- order(abs(ticks))
    ticks  <- abs(ticks)[o]
    tsigns <- tsigns[o]
    if(ticks[1]==0){
      add_s <- 2
      d_start<- -1
      ticks<- ticks[-1]
      tsigns <- tsigns[-1]
    }
    if(length(ticks)>0){
      if(ticks[length(ticks)]==1){
        add_e <- 2
        d_end<- -1
        ticks  <- ticks[-length(ticks)]
        tsigns <- tsigns[-length(tsigns)]
      }
    }
  }
  nt     <- length(ticks)
  np     <- 2+2*nt
  rx     <- (1/np)*phi
  md     <- min(diff(c(0, ticks, 1)))
  if(md<(rx*2)) rx<- md/2
  if(type==1) p <- -rev(exp(seq(0,5,length.out=round(n*rx))))
  if(type==2) p <-  sqrt(seq(0,5,length.out=round(n*rx)))
  if(type==3) p <- seq(0,1,length.out=round(n*rx))
  if(type==4) p <- c(0,rep(1, round(n*rx)-1))
  if(type==5) p <- -rev((seq(0,5,length.out=round(n*rx)))^2)
  p      <- a_st(p)
  pb     <- length(p)
  sy <- c(p*d_start+add_s)
  location <- pb
  if(nt>0){
    for(i in 1:nt){
      on <- round(n*ticks[i])
      add<- 2
      if(tsigns[i]==-1) add<- 0
      sy <- c(sy, rep(1, on-location-pb), tsigns[i]*(-rev(p)+add), tsigns[i]*(-p+add))
      location <- on+pb
    }
  }
  sy <- c(sy, rep(1, n-location-pb), rev(p)*d_end+add_e)
  sy <- a_st(sy)
  sx <- seq(0, 1, length.out=length(sy))
  rbind(sx, sy)
}


a_st <-
function(x, na.rm=F)
{
  x<-x-min(x, na.rm=na.rm)
  x<-x/max(x, na.rm=na.rm)
  x
}


grid.brackets <-
function(x1, y1, x2, y2, h=NULL, ticks=0.5, curvature=0.5, type=1, col=1, lwd=1, lty='solid')
{

  is.unit = grid::is.unit
  unit = grid::unit
  grid.lines = grid::grid.lines
  gpar = grid::gpar

  if(is.null(h)) h <- 0.05
  if(!is.unit(x1) | !is.unit(y1) | !is.unit(x2) | !is.unit(y2) | !is.unit(h)){
    x1 <- unit(x1, 'native')
    x2 <- unit(x2, 'native')
    y1 <- unit(y1, 'native')
    y2 <- unit(y2, 'native')
  }

  if(!is.numeric(curvature)) stop('curvature must be numeric')
  if(!is.numeric(type))      stop('type must be a integer, 1 to 5')
  if(length(ticks)==1)       if(is.na(ticks)) ticks<- NULL
  if(!is.numeric(ticks) & !is.null(ticks))     stop('ticks must be numeric or NULL')
  if(length(ticks)>1){
  if(any(duplicated(abs(ticks)))) stop('duplicated ticks')
  }
  if(curvature<0) curvature<- 0
  if(curvature>1) curvature<- 1

  xd <- (x2-x1)
  yd <- (y2-y1)

  brackets<- a_cb_brackets(phi=curvature, ticks=ticks, type=type)
  grid.lines(x1+xd*brackets[1,], y1+h*brackets[2,], gp=gpar(col=col, lwd=lwd, lty=lty))
}



draw_bracket = function(x,xend,y,size=1,offset=0) {
  if ( ! grid::is.unit(offset) ) {
    offset = grid::unit(offset*.pt,"mm")
  }
  grid.brackets(grid::unit(x,"native"),grid::unit(y,"native")+offset,grid::unit(xend,"native"),grid::unit(y,"native")+offset,grid::unit(size*.pt,"mm"),ticks=NULL,curvature=0.2)
}

.pt <- 72.27 / 25.4

#' @export
GeomBracket <- ggplot2::ggproto("GeomBracket", ggplot2::Geom,
                        required_aes=c('x','xend'),
                        draw_panel = function(data, panel_scales, coord,size=1,offset=0) {
                          coords <- coord$transform(data, panel_scales)
                          draw_bracket_vec = Vectorize(draw_bracket,SIMPLIFY=F)
                          results = draw_bracket_vec(coords$x,coords$xend,coords$y,size,rep(offset,nrow(coords)))
                          do.call(grid::gList,results)
                        }
)
