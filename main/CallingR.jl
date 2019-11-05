using RCall

include("Driver.jl")


########################################################################
#Plotting States (Figure 2)
########################################################################
solFigure2 = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)

allStateData = DataFrame()
for (i,st) in enumerate(statesNames)
    allStateData[!,Symbol(st)] = vec(solFigure2[:,:,i,:])
end

timeLength = length(solFigure2.t)
allStateData.Cell = repeat(1:nCells,timeLength)
allStateData.Time = repeat(solFigure2.t,inner=nCells)
allStateData.Infected = repeat(vec(u0[:,:,2] .> 0.0),timeLength)

@rput allStateData
R"""
library(ggplot2)
library(ggpubr)

stateAve = aggregate(allStateData[,1:14], list(allStateData$Infected,allStateData$Time), mean)
colnames(stateAve)[1] <- "Cell.State"
colnames(stateAve)[2] <- "Time"

logic<- unlist(lapply(stateAve$Cell.State, function(x) x == TRUE))
stateAve$Cell.State[which(logic == TRUE)] <- "Primary"
stateAve$Cell.State[which(logic == FALSE)] <- "Secondary"
low = aggregate(allStateData[,1:14], list(allStateData$Infected,allStateData$Time), FUN = 'quantile',probs=0.05)
high = aggregate(allStateData[,1:14], list(allStateData$Infected,allStateData$Time), FUN = 'quantile',probs=0.95)

commonFigureOptions <- list( scale_x_continuous(breaks=seq(0, 48, 12)),
  theme_pubr(border=TRUE),
  ylab("nM"),
  theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1))

p1 <- ggplot(stateAve) + geom_line(aes(y=DNA, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=low$DNA, ymax=high$DNA, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("DNA") +
  commonFigureOptions

p2 <- ggplot(stateAve) + geom_line(aes(y=IFNb, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=low$IFNb, ymax=high$IFNb, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IFNb") +
  commonFigureOptions

p3 <- ggplot(stateAve) + geom_line(aes(y=IRF7, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=low$IRF7, ymax=high$IRF7, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("IRF7") +
  commonFigureOptions

p4 <- ggplot(stateAve) + geom_line(aes(y=TREX1, x=Time, group=Cell.State, color = Cell.State))+
  geom_ribbon(aes(ymin=low$TREX1, ymax=high$TREX1, x=Time,group=Cell.State,fill = Cell.State), alpha = 0.3) +
  ggtitle("TREX1") +
  commonFigureOptions

figure <- ggarrange(p1, p2, p3, p4,
                    labels = c("A", "B", "C", "D"),
                    common.legend = TRUE, legend = "right",
                    align = "hv",
                    ncol = 2, nrow = 2)

ggsave("../Figures/Figure2.pdf")
"""

########################################################################
#Plot Headtmaps (Figure 3)
########################################################################

solFigure2 = @time solve(prob,CVODE_BDF(linear_solver=:GMRES),saveat=0.1)

heat10hr = DataFrame()
heat10hr.x = repeat(1:N,N)
heat10hr.y = repeat(1:N,inner=N)
heat10hr.IFN = vec(solFigure2(10.0)[:,:,7])

heat10hrVirus = CSV.read("heat10hrVirus.csv")

@rput heat10hr
@rput heat10hrVirus
R"""
  library(ggplot2)
  library(ggpubr)

  commonFigureOptions <- list(scale_fill_distiller(palette = "Spectral",guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")),
    scale_x_continuous(expand=c(0,0)),
    scale_y_continuous(expand=c(0,0)),
    theme_bw(base_size = 14),
    ylab("Cell"),
    xlab("Cell"),
    labs(fill="IFN (nM)"),
    theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1))

  p1 <- ggplot(heat10hr, aes(x, y, fill=IFN)) +
    geom_raster(aes(fill=IFN)) +
    labs(fill="IFN (nM)") +
    ggtitle("Drop (ISD): \n Time = 10 hours") +
    commonFigureOptions

  p2 <- ggplot(heat10hrVirus, aes(x, y, fill=IFN)) +
    geom_raster(aes(fill=IFN)) +
    ggtitle("Wash (Virus): \n Time = 10 hours") +
    commonFigureOptions

  figure <- ggarrange(p1, p2,
                      labels = c("A", "B"),
                      common.legend = TRUE, legend = "right",
                      align = "h",
                      ncol = 2, nrow = 2)

  ggsave("../Figures/Figure3.pdf")
"""


########################################################################
#Plot Knockdown (Figure 4)
########################################################################


@rput KDData
R"""
library(ggplot2)
library(ggpubr)

stateAve = aggregate(KDData[,1:2], list(KDData$Percent,KDData$Time), mean)
colnames(stateAve)[1] <- "Percent"
colnames(stateAve)[2] <- "Time"


low = aggregate(KDData[,1:2], list(KDData$Percent,KDData$Time), FUN = 'quantile',probs=0.05)
high = aggregate(KDData[,1:2], list(KDData$Percent,KDData$Time), FUN = 'quantile',probs=0.95)

commonFigureOptions <- list(scale_x_continuous(breaks=seq(0, 48, 12)),
  theme_pubr(border=TRUE),
  ylab("nM"),
  theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 1))

p1 <- ggplot(stateAve) + geom_line(aes(y=TREXKD, x=Time, group=factor(Percent), color = factor(Percent))) +
  geom_ribbon(aes(ymin=low$TREXKD, ymax=high$TREXKD, x=Time,group=factor(Percent),fill = factor(Percent)), alpha = 0.2) +
  ggtitle("TREX1 Knockdown") +
  commonFigureOptions

p2 <- ggplot(stateAve) + geom_line(aes(y=IRF7KD, x=Time, group=factor(Percent), color = factor(Percent))) +
  geom_ribbon(aes(ymin=low$IRF7KD, ymax=high$IRF7KD, x=Time,group=factor(Percent),fill = factor(Percent)), alpha = 0.2) +
  ggtitle("IRF7 Knockdown") +
  commonFigureOptions

  figure <- ggarrange(p1, p2,
                      labels = c("A", "B"),
                      common.legend = TRUE, legend = "right",
                      align = "h",
                      ncol = 2, nrow = 2)

ggsave("../Figures/Figure4.pdf")
"""


#full = vcat(primary, secondary)
#CSV.write("MOI1e-2.csv",full)
lowInfect = CSV.read("MOI1e-2.csv")
medInfect = CSV.read("MOI1e-1.csv")
highInfect = CSV.read("MOI1e0.csv")

@rput lowInfect
@rput medInfect
@rput highInfect
R"""
library(ggplot2)
library(ggpubr)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

commonFigureOptions <- list(geom_split_violin(scale="width"),
  theme_pubr(border=TRUE),
  xlab("Initial Condition Variance"),
  ylab("IFN (nM)"),
  scale_x_discrete(labels= unlist(lapply(seq(0.0,1.0,0.1),toString))),
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = 0.25, position = position_dodge(width = .25)),
  theme(plot.title = element_text(hjust = 0.5),aspect.ratio = 0.3)
  )


p1 <- ggplot(lowInfect, aes(factor(percent), IFN, fill = CellState)) +
  ggtitle("1% Primary Cells Initially") +
  commonFigureOptions

p2 <- ggplot(medInfect, aes(factor(percent), IFN, fill = CellState)) +
  ggtitle("10% Primary Cells Initially") +
  commonFigureOptions

p3 <- ggplot(highInfect, aes(factor(percent), IFN, fill = CellState)) +
  ggtitle("63% Primary Cells Initially") +
  commonFigureOptions

figure <- ggarrange(p1, p2, p3,
                    labels = c("A","B","C"),
                    common.legend = TRUE, legend = "right",
                    align = "h",
                    ncol = 1, nrow = 3)

ggsave("../Figures/Figure5.pdf")
"""
