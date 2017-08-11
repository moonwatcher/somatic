#!/usr/local/bin/Rscript

library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
data_filename = args[1]
title = args[2]
diagram_filename = gsub('.csv', '.pdf', data_filename)

intensity = read.csv(data_filename, header=T)
intensity_plot <- ggplot(intensity, aes(index, log(intensity + 1, 2), fill=pivot)) +
geom_bar(stat="identity", position="dodge") + facet_grid(study ~ .) +
scale_alpha_manual(name="Sample", values=c(0.4, 1.0)) +
scale_fill_manual(name="Region", values=c('#0B345B', '#D97A04', '#581A4D')) +
scale_x_continuous(breaks = intensity$index[which(intensity$index%%2 == 0)], labels = intensity$start[which(intensity$index%%2 == 0)]) +
ggtitle(title) +
xlab("chromosome") + 
ylab("log2(expression + 1)") +  
theme (
    panel.background = element_rect(colour = "#FFFFFF", fill="#E4E4E4"), 
    axis.text.x = element_text(angle = 90, size = rel(0.8)),
    axis.title = element_text(size = rel(1), color="#444444"),
    panel.grid.minor = element_line(colour = "#F0F0F0", size = rel(0.25), linetype = "solid"),
    panel.grid.major = element_line(colour = "#F2F2F2", size = rel(0.4), linetype = "solid"),
    plot.title = element_text(vjust = 1)
)
ggsave(diagram_filename, plot = intensity_plot, width = 16, height = length(unique(intensity$study)) * 3, limitsize=FALSE)
