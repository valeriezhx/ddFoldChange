# ddFoldChange

R package

This is a R package to analyze RT-PCR results by using delta-delta CT (ΔΔCT) method.

Install the package from GitHub in R:

  devtools::install_github("valeriezhx/ddFoldChange", build_vignettes = TRUE, force=TRUE)

library(ddFoldChange)

There are two steps to analysis the data

run .csv input file prepared in the correct format and obtain a new .csv file containing gene fold changes
ProcessData("xxx.csv", replicate=X, output="yyy.csv")

plot bar graph for the gene fold changes relative to the control group
PlotFC(df.process)
