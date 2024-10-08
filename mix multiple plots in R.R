# Save ggplots with semi-transparent colors
#Use cairo-based postscript graphics devices 
#adapted from
#https://www.r-bloggers.com/2017/07/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
p <- plot(iris)
plot(iris)
#Using build in tools in Rstudio
#plot -> Save as /Export
##################################################################
# JPEG device
jpeg("my_plot.jpeg", quality = 100)
# iris
plot(iris)
# Close device
dev.off()
##################################################################
#Saving Image usingggsave()
#ggsave(
#  filename,
#  plot = last_plot(),
#  device = NULL,
#  path = NULL,
#  scale = 1,
#  width = NA,
#  height = NA,
#  units = c("in", "cm", "mm"),
#  dpi = 300,
#  limitsize = TRUE,
#  ...
#  )
library(ggplot2)
p <- ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width)) + geom_point()
p
ggsave("ggplot2save.jpg")
ggsave("ggplot2save.png")
ggsave("ggplot2save.tiff")
ggsave("ggplot2save.pdf")
#################################################################
#Changing Image size and DPI
ggsave("ggplot2save.jpg", width = 10, height = 8, units = c("cm"))
ggsave("ggplot2save.jpg", width = 10, height = 8, units = c("cm"), dpi = 300)
##plotting mulitple graph onto the same image using ggrid()
##################################################################
#base R
jpeg("hist_gpa_sat.jpg")
par(mfrow=c(2,1))
hist(gpa)
hist(sat)
dev.off()

#using ggpubr
library(ggpubr)
p <- ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width)) + geom_point()
q <- ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width)) + geom_line()
r <- ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width)) + geom_jitter()
s <- ggplot(iris,aes(x=Petal.Length,y=Petal.Width)) + geom_jitter()
ggarrange(p, q, r,s ,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
ggsave("ggpubrsave.jpg", width = 20, height = 16, units = c("cm"), dpi = 300)
###################################################################
#Using corplot
library("cowplot")
ggdraw() +
  draw_plot(p, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(q, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(r, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))
ggsave("cowplot.ggdraw.jpg", width = 20, height = 16, units = c("cm"), dpi = 300)
###################################################################
##Something a bit more Complicated
# Scatter plot colored by groups ("Species")
sp <- ggscatter(iris, x = "Sepal.Length", y = "Sepal.Width",
                color = "Species", palette = "jco",
                size = 3, alpha = 0.6)+
  border()
sp
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(iris, "Sepal.Length", fill = "Species",
                   palette = "jco")
yplot <- ggdensity(iris, "Sepal.Width", fill = "Species",
                   palette = "jco")+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme()
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, sp, yplot,
          ncol = 2, nrow = 2,  align = "hv",
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)
##plotting mulitple graph onto the same image using ggrid()
##################################################################
#base R
jpeg("Base.R.Multiple.jpg")
par(mfrow=c(1,2))# one row 2 columns
hist(rnorm(200)) # this is the first plot in to first place in table
hist(rnorm(200)) # this is the second plot in to second place in table
dev.off()
getwd()
#adapted from
#https://www.r-bloggers.com/2017/07/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
plot(iris)
#Using build in tools in Rstudio
#plot -> Save as /Export
##################################################################
# JPEG device
jpeg("my_plot.jpeg", quality = 75)
# iris
plot(iris)
# Close device
dev.off()
##################################################################
#Saving Image usingggsave()
#ggsave(
#  filename,
#  plot = last_plot(),
#  device = NULL,
#  path = NULL,
#  scale = 1,
#  width = NA,
#  height = NA,
#  units = c("in", "cm", "mm"),
#  dpi = 300,
#  limitsize = TRUE,
#  ...
#  )
library(ggplot2)
p <- ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width)) + geom_point()
p
ggsave("ggplot2save.jpg")
ggsave("ggplot2save.png")
ggsave("ggplot2save.tiff")
ggsave("ggplot2save.pdf")
#################################################################
#Changing Image size and DPI
ggsave("ggplot2save.jpg", width = 10, height = 8, units = c("cm"))
ggsave("ggplot2save.jpg", width = 10, height = 8, units = c("cm"), dpi = 300)
##plotting mulitple graph onto the same image using ggrid()
##################################################################
#base R
jpeg("Base.R.Multiple.jpg")
par(mfrow=c(1,2))
hist(rnorm(200))
hist(rnorm(200))
dev.off()
#using ggpubr
library(ggpubr)
p <- ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width)) + geom_point()
q <- ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width)) + geom_line()
r <- ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width)) + geom_jitter()
s <- ggplot(iris,aes(x=Petal.Length,y=Petal.Width)) + geom_jitter()
ggarrange(p, q, r,s ,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
ggsave("ggpubrsave.jpg", width = 20, height = 16, units = c("cm"), dpi = 300)
###################################################################
#Using corplot
library("cowplot")
ggdraw() +
  draw_plot(p, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(q, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(r, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))
ggsave("cowplot.ggdraw.jpg", width = 20, height = 16, units = c("cm"), dpi = 300)
###################################################################
##Something a bit more Complicated
# Scatter plot colored by groups ("Species")
sp <- ggscatter(iris, x = "Sepal.Length", y = "Sepal.Width",
                color = "Species", palette = "jco",
                size = 3, alpha = 0.6)+
  border()
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(iris, "Sepal.Length", fill = "Species",
                   palette = "jco")
yplot <- ggdensity(iris, "Sepal.Width", fill = "Species",
                   palette = "jco")+
  rotate()
# Cleaning the plots
yplot <- yplot + clean_theme()
xplot <- xplot + clean_theme()
# Arranging the plot
ggarrange(xplot, NULL, sp, yplot,
          ncol = 2, nrow = 2,  align = "hv",
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)
plot(iris)
# JPEG device
jpeg("my_plot.jpeg", quality = 75)
# iris
plot(iris)
# Close device
dev.off()
getwd()
library(ggplot2)
p <- ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width)) + geom_point()
p
ggsave("ggplot2save.jpg")
ggsave("ggplot2save.png")
ggsave("ggplot2save.tiff")
ggsave("ggplot2save.pdf")
ggsave("ggplot2save.jpg", width = 10, height = 8, units = c("cm"), dpi = 300)
##plotting mulitple graph onto the same image using ggrid()
##################################################################
#base R
jpeg("Base.R.Multiple.jpg")
par(mfrow=c(1,2))
hist(rnorm(200))
hist(rnorm(200))
dev.off()
#using ggpubr
library(ggpubr)
p <- ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width)) + geom_point()
q <- ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width)) + geom_line()
r <- ggplot(iris,aes(x=Sepal.Length,y=Sepal.Width)) + geom_jitter()
s <- ggplot(iris,aes(x=Petal.Length,y=Petal.Width)) + geom_jitter()
ggarrange(p, q, r,s ,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
ggsave("ggpubrsave.jpg", width = 20, height = 16, units = c("cm"), dpi = 300)
###################################################################
#Using cowplot
library("cowplot")
###################################################################
#Using cowplot
library("cowplot")
ggdraw() +
  draw_plot(p, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(q, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(r, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))
ggsave("cowplot.ggdraw.jpg", width = 20, height = 16, units = c("cm"), dpi = 300)
# Scatter plot colored by groups ("Species")
sp <- ggscatter(iris, x = "Sepal.Length", y = "Sepal.Width",
                color = "Species", palette = "jco",
                size = 3, alpha = 0.6)+
  border()
sp
# Marginal density plot of x (top panel) and y (right panel)
xplot <- ggdensity(iris, "Sepal.Length", fill = "Species",
                   palette = "jco")
xplot
yplot()
yplot <- ggdensity(iris, "Sepal.Width", fill = "Species",
                   palette = "jco")+
  rotate()
yplot()
# Cleaning the plots
yplot <- yplot + clean_theme()
xplot <- xplot + clean_theme()
yplot
# Arranging the plot
ggarrange(xplot, NULL, sp, yplot,
          ncol = 2, nrow = 2,  align = "hv",
          widths = c(2, 1), heights = c(1, 2),
          common.legend = TRUE)