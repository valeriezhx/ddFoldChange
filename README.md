ddFoldChange  

R package

This is a R package to analyze RT-PCR results by using delta-delta CT (ΔΔCT) method.

Install the package from GitHub in R:

devtools::install_github("valeriezhx/ddFoldChange_0.0.0.9000.tar")

library(ddFoldChange)

There are two steps to analysis the data

1) run .csv input file prepared in the correct format and obtain a new .csv file containing gene fold changes 

process.data("xxx.csv", replicate=X, output="yyy.csv")


2) plot bar graph for the gene fold changes relative to the control group

df.process <- read.csv("yyy.output", header = TRUE, sep = ',')

plot_FC(df.process)
