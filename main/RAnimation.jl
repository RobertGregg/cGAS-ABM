using RCall

include("Driver.jl")


timeLength = length(sol.t)

heatAll = DataFrame()
heatAll.x = repeat(repeat(1:N,N),timeLength)
heatAll.y = repeat(repeat(1:N,inner=N),timeLength)
heatAll.t = repeat(sol.t,inner=N^2)
heatAll.IFN = vec(sol[:,:,7,:])

@rput heatAll
R"""
library(ggplot2)
library(ggpubr)
library(gganimate)

commonFigureOptions <- list(scale_fill_distiller(palette = "Spectral",guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")),
                            scale_x_continuous(expand=c(0,0)),
                            scale_y_continuous(expand=c(0,0)),
                            theme_bw(base_size = 14),
                            ylab("Cell"),
                            xlab("Cell"),
                            labs(fill="IFN (nM)"),
                            theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1))

p1 <- ggplot(heatAll, aes(x, y, fill=IFN)) +
  geom_raster(aes(fill=IFN)) +
  labs(fill="IFN (nM)") +
  ggtitle("Drop (ISD): \n Time = 10 hours") +
  commonFigureOptions +
  labs(title = 'Time: {format(round(frame_time, 2), nsmall = 2)} hours') +
  transition_time(t)

anim_save("InterferonAnimation.gif",animation = p1)
"""
