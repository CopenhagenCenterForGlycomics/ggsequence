get_conservation <- function(data) {

 strong_res = c("[STA]",
 "[NEQK]",
 "[NHQK]",
 "[NDEQ]",
 "[QHRK]",
 "[MILV]",
 "[MILF]",
 "[HY]",
 "[FYW]")
 
 weak_res = c("[CSA]",
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
			aas$cons <- 0.75
		} else if ( max(sapply(strong_res,function(re) { length(grep(re,aas$aa)) })) == length(aas$aa) ) {
			aas$cons <- 0.5
		} else if ( max(sapply(weak_res,function(re) { length(grep(re,aas$aa)) })) == length(aas$aa) ) {
			aas$cons <- 0.1
		} else {
			aas$cons <- 0
		}
		aas
	})$cons
	consvals
}

overlay_conservation <- function(alignment.plot) {
	alignment.plot + geom_tile(aes(x=pos,y=seqid),fill='black',alpha=get_conservation(alignment.plot$data))
}


#' Wrap the venn diagram into a geom
#' @export
stat_conservation <- function(mapping = NULL, data = NULL, stat = "identity",
                          position = "identity",
                          show.legend = NA, inherit.aes = TRUE,na.rm=T,alone=T,...) {
  ggplot2::layer(
    data = data,
    mapping = ggplot2::aes(seqname=seqname,x=pos,y=seqname,pos=pos,start=start,end=end,seq.id=seq.id,aa=aa),
    stat = StatConservation,
    geom = "tile",
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
                        default_aes = ggplot2::aes(alpha=..conservation..),
                        setup_data = function(data,params) {
                        	data[data$seqname != 'conservation',]
                        },
                        compute_panel = function(data,scales,alone=T) {
                        	data$conservation = get_conservation(data)
                        	data$pos = as.numeric(data$pos)
                        	if (alone) {
                        		data$seqname = factor(rep("conservation",nrow(data)),levels=levels(data$seqname) )
                        		data$y = as.numeric(data$seqname)
                        	}
							return(unique(data))
                        }
)