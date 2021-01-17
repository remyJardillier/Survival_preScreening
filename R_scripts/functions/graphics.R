pal <- wes_palette("Zissou1", 100, type = "continuous")

# Parameters for the the graphics -----------------------------------------

# @description: theme for ggplot
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(family = 'Helvetica'), 
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(fill = NA, colour = "black", size = rel(1)),
            axis.title = element_text(face = "bold",size = rel(1.2)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = 15), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="lightgrey", linetype = "dashed", size = 0.3),
            panel.grid.minor = element_line(colour="#f0f0f0", linetype = "dashed", size = 0.3),
            legend.key = element_rect(colour = NA),
            legend.position = "none",
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold"),
            plot.subtitle = element_text(hjust = 0.5)
    ))
  
}

# @description: theme for ggplot with legend on the right
theme_Publication_legend_right <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(family = 'Helvetica'), 
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(fill = NA, colour = "black", size = rel(1)),
            axis.title = element_text(face = "bold",size = rel(1.2)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = 15), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="lightgrey", linetype = "dashed", size = 0.3),
            panel.grid.minor = element_line(colour="#f0f0f0", linetype = "dashed", size = 0.3),
            legend.key = element_rect(colour = NA),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold"), 
            plot.subtitle = element_text(hjust = 0.5),
            legend.title = element_text(size = 18),
            legend.text = element_text(size = 18)
    ))
  
}


# @description: parameters for "fill" in ggplot
scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#a6cee3","#ef3b2c","#662506","#fb9a99","#984ea3","#ffff33")), ...)
  
}

# @description: parameters for "fill" in ggplot for grids (figure 2A)
scale_fill_Publication_grid <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#FF0000FF", "#FF8000FF", "#FFFF00FF", "#FFFF80FF", "blue")), ...)
  
}

# @description: parameters for "fill" in ggplot only for selection methods
scale_fill_Publication_slc <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#fdb462","#7fc97f","#a6cee3","#ef3b2c","#662506","#fb9a99","#984ea3","#ffff33")), ...)
  
}

# @description: parameters for "colour" in ggplot
scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}


# @description: parameters for 'classic' (non-ggplot) graphics
par(font.axis = 1.5, # font for x and y axis
    font.lab=1, # font for x and y label
    cex.lab = 1.5, # magnification for x and y label
    cex.axis = 1.5, # magnification for x and y axis
    mar=c(5.1, 4.5, 3.1, 2.1) # margins (bottom, left, top, right)
)



# stars (p-values) and colors ---------------------------------------------

# 
# @description: Stars according to p-value
#
# @parameters:
#   - p_val: a p-value
#  
# @return:
#   - string containing stars or 'n.s' according to the significancer level
my_stars.pval <- function(p_val){
  
  if(is.na(p_val)){
    return(NA)
  }else if(p_val <= 0.001){ 
    return("***")
  }else if(p_val <= 0.01 & p_val > 0.001){
    #return("\U1F7B6\U1F7B6")
    return("**")
  }else if(p_val <= 0.05 & p_val > 0.01){
    return("*")
  }else if(p_val <= 0.1 & p_val > 0.05){
    return("+")
  }else{
    return("n.s.")
  }
}

my_stars.pval_vect <- function(x){
  y <- sapply(x, my_stars.pval)
  return(y)
}

my.col <- function(p_val){
  
  if(is.na(p_val)){
    return(NA)
  }else if(p_val <= 0.001){ 
    return("darkred")
  }else if(p_val <= 0.01 & p_val > 0.001){
    return("red")
  }else if(p_val <= 0.05 & p_val > 0.01){
    return("darkorange")
  }else if(p_val <= 0.1 & p_val > 0.05){
    return("#FFCC33")
  }else{
    return("grey")
  }
}

