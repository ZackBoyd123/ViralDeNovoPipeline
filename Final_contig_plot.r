#!/usr/bin/Rscript
args = commandArgs(trailingOnly = TRUE)
coverage_and_genome<- read.csv(args[1], sep = ",")
contigs<- read.csv(args[2],sep = ",")
#Subset to genome data
genome<- coverage_and_genome[grep("GenomePlot", coverage_and_genome$Plot), ]
#Subset to coverage data
coverage<- coverage_and_genome[grep("CoveragePlot", coverage_and_genome$Plot), ]
library(ggplot2)
library("cowplot", lib.loc="~/RPackages")
#Genome chart
geomchart<- ggplot() + geom_rect(data=genome, aes(xmin=X.Min, xmax=X.Max, ymin=Y.Min, ymax=Y.Max, fill=Figure, colour=Figure)) + theme(axis.line.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), legend.text = element_text(size = 8), legend.title = element_blank()) + scale_y_continuous(limits = c(0,5))

#Coverage Plot
covplot<- ggplot(coverage, aes(x=coverage$X.Min))
coverageplot<- covplot + geom_area(aes(y=coverage$X.Max, colour=Figure), fill="#F8766D") + ylab("") + xlab("Genome position")+ scale_y_reverse() + theme(axis.line.x = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), legend.title = element_blank(), legend.text = element_text(size = 8))

#Plot of contigs
contigplot<- ggplot(contigs, aes(xmin=X.Min,xmax=X.Max,ymin=Y.Min,ymax=Y.Max, colour=Figure,fill=Figure)) + geom_rect(alpha=0.1,fill=alpha(0.1)) + xlab("Genome Location")+ theme(legend.title = element_blank(), legend.text = element_text(size = 8))
contigplot<- contigplot + facet_wrap(~Aligner, ncol= 1, scales = "free_y")    #, strip.position = "right")
#Final plot, all aligned
final_plot<- plot_grid(geomchart,coverageplot,contigplot,align = "v",nrow = 3, rel_heights = c(1/12,1/12,5/6))

pdf(args[3], width = 20, height = 15)
final_plot
dev.off()

