#!/usr/bin/Rscript
#setwd("~/Documents/Zack/PyCharm/Internship/Denovo/")

args = commandArgs(trailingOnly=TRUE)
mydata<- read.csv(args[1], sep=",")

#mydata<- read.csv("Mix-13_stats.csv", sep=",")
#setwd("~/Documents/Zack/")

#Subset to genome data
genome<- mydata[grep("GenomePlot", mydata$Plot), ]
#Subset to pie data
pie<- mydata[grep("PiePlot", mydata$Plot), ]
colnames(pie)[1] <- "Base"
#Subset to coverage data
coverage<- mydata[grep("CoveragePlot", mydata$Plot), ]
#Subset to contig data
contigs<- mydata[grep("ContigPlot",mydata$Plot), ]


#Set libraries
library(ggplot2)
library(gridExtra)
library("cowplot", lib.loc="~/RPackages")

geomchart<- ggplot() + geom_rect(data=genome, aes(xmin=genome$X.Min, xmax=genome$X.Max, ymin=genome$Y.Min, ymax=genome$Y.Max, fill=Figure, colour=Figure)) + theme(axis.line.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank()) + scale_y_continuous(limits = c(0,5))
#geomchart

#Make pie chart
bp<- ggplot(pie, aes(x="", y=pie$X.Min, fill=Base)) + geom_bar(width = 1, stat = "identity")
piechart<- bp + ylab("") + xlab("") +
  geom_text(aes(label = pie$X.Min), position = position_stack(vjust = 0.5), size = 3) +
  coord_polar(theta = "y")

#piechart

#Make coverageplot
covplot<- ggplot(coverage, aes(x=coverage$X.Min))
coverageplot<- covplot + geom_area(aes(y=coverage$X.Max, colour=Figure), fill="#F8766D") + ylab("") + xlab("Genome position")+ scale_y_reverse() + theme(axis.line.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank())
#coverageplot

#Make contig plot
#contigplot<- ggplot() + geom_rect(data=contigs, aes(xmin=contigs$X.Min, xmax=contigs$X.Max, ymin=contigs$Y.Min, ymax=contigs$Y.Max, fill=Figure, colour=Figure, alpha = 0.1)) + theme(axis.text.y = element_blank()) + xlab("Genome Location")
#contigplot 
contigplot<- ggplot(contigs, aes(xmin=contigs$X.Min,xmax=contigs$X.Max,ymin=contigs$Y.Min,ymax=contigs$Y.Max, colour=Figure)) + geom_rect(alpha=0.2,fill=alpha("grey",0.2)) + theme(axis.text.y = element_blank()) + xlab("Genome Location")

test<- plot_grid(geomchart, coverageplot, contigplot, align="v", nrow = 3, rel_heights = c(1/3,1/3,1/3))

#pdf("plot.pdf",width = 20, height = 7)
pdf(args[2],width = 20, height = 7)
plot_grid(piechart, test, ncol = 1, rel_heights = c(0.8,2))
dev.off()



