
#' read_vplot
#'
#' @param file .VMat file 
#' @param xlims limits for position axis
#' @return returns matrix with vplot informatiom
#' @seealso \code{\link{plotV}} 
#' @export
read_vplot<-function(file, xlims= NA){
  mat = read.delim(file,skip = 7, header=F)
  i_lower = scan(file,skip = 3,nmax=1,quiet=T)
  i_upper = scan(file,skip = 5,nmax=1,quiet=T)-1
  rownames(mat) = i_lower:i_upper
  if (is.na(xlims)){
    w = round(dim(mat)[2]/2)
    colnames(mat) = -w:w
  }
  else{
    colnames(mat) = xlims[1]:xlims[2]
  }
  return(mat)
}

#' vplot_theme
#'
#' @param base_size size for text
#' @param base_family font family
#' @return returns theme for ggplot
#' @seealso \code{\link{plotV}} 
#' @import ggplot2
#' @import grid
#' @export
vplot_theme <-function(base_size = 7, base_family="Helvetica"){
  theme(line = element_line(colour = "black", size = 0.5, linetype = 1, 
                            lineend = "butt"), 
        rect = element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1), 
        text = element_text(family = base_family, face = "plain", colour = "black", size = base_size, hjust = 0.5, 
                            vjust = 0.5, angle = 0, lineheight = 0.9), 
        strip.text = element_text(size = rel(0.8)), 
        axis.line = element_blank(),#element_line(colour = "black", size = 0.5), 
        axis.text = element_text(size = rel(0.8), colour = "black"),
        axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(hjust = 1), 
        axis.ticks = element_line(colour = "black", size=0.25), 
        axis.title.x = element_text(), 
        axis.title.y = element_text(angle = 90), 
        axis.ticks.length = unit(0.10, "cm"), 
        axis.ticks.margin = unit(0.1, "cm"), 
        legend.background = element_rect(colour = NA), 
        legend.margin = unit(0.2, "cm"), 
        legend.key = element_blank(), 
        legend.key.size = unit(1.2, "lines"), 
        legend.key.width = unit(0.05, "inches"), 
        legend.key.height = unit(0.2, "inches"),
        legend.text = element_text(size = rel(0.8)), 
        legend.text.align = NULL, 
        legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0), 
        legend.title.align = NULL, 
        legend.position = "right", 
        legend.direction = NULL, 
        legend.justification = "center", 
        legend.box = NULL, 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA,color = "black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.margin = unit(0.25, "lines"), 
        panel.margin.x = NULL, 
        panel.margin.y = NULL, 
        strip.background = element_blank(), 
        strip.text.x = element_text(), 
        strip.text.y = element_text(angle = -90), 
        plot.background = element_blank(), 
        plot.title = element_text(size = rel(1.2)), 
        plot.margin = unit(c(0.5, 0.5, 0.25, 0.25), "lines"), complete = TRUE)
}


#' plotV
#'
#' @param X vplot matrix
#' @param xlabel label for x axis
#' @param ylabel label for y axis
#' @param guide label for legend
#' @param name title
#' @param palette color palette
#' @param limits limits for color
#' @param xbreaks where to include breaks for x axis
#' @param ylabel label for y axis
#' @return returns ggplot object
#' @seealso \code{\link{read_vplot}} \code{\link{vplot_theme}} 
#' @import ggplot2
#' @export
plotV <- function(X, xlabel = "Center position relative to dyad (bp)", ylabel = "Fragment size (bp)", guide = "Density", name = NA, 
                   palette = "BuPu", limits = NA, xbreaks = NA, ybreaks = NA){
  df = cbind(data.frame("y" = factor(rownames(X), levels = rownames(X), ordered=T)),X)
  mdf = melt(df, id = "y")
  p = ggplot(mdf, aes(x=variable, y=y, col = value)) + 
    geom_raster(aes(fill=value), interpolate=T) + 
    coord_fixed()+
    xlab(xlabel) + ylab(ylabel)+
    vplot_theme(7)
  if (is.na(limits[1])){
    limits = c(0, max(mdf$value))
  }
  p = p + scale_fill_gradientn(colours = c("white",colorRampPalette(brewer.pal(9,palette))(9)), name = guide, 
                               limits=limits, breaks = limits, labels = pretty_scale_format, expand=c(0,0), guide = guide_colorbar(ticks = F))+
    scale_colour_gradientn(colours = c("white",colorRampPalette(brewer.pal(9,palette))(9)), name = guide, 
                           limits=limits,breaks = limits, labels = pretty_scale_format,expand=c(0,0), guide = F)
  if (is.na(xbreaks[1])){
    xbreaks = seq(as.numeric(colnames(X)[1]),as.numeric(colnames(X)[ncol(X)]),10)
  }
  if (is.na(ybreaks[1])){
    ybreaks = seq(as.numeric(rownames(X)[1]),as.numeric(rownames(X)[nrow(X)]),10)
  }  
  p = p + scale_x_discrete(breaks = xbreaks) + 
    scale_y_discrete(breaks = ybreaks)
  if (!is.na(name)){
    p = p + ggtitle(name) + theme( plot.title = element_text(size = 6, colour = "black"))
  }
  return(p)
}
