##' Plot Fold Change for Target Genes
##'
##' We want to use ggplot to graph the fold change for each of the target genes.
##' @title Bar Graph for Gene Fold Change
##' @author Xia Zhao
##' @param df.process The output data created from the `ProcessData` function, including the fold change
##'
##'  \describe{
##'     \item{trt}{A factor column indicating treatment groups.}
##'     \item{Fold Change Columns}{Fold change values for each target gene (e.g., `FC.<gene>`).}
#'   }
##' @import ggplot2
##' @examples
##' # Example usage with internal CSV data bundled with the package
##' # Step 1: Load the data and process it to get df.process
##' file_path <- system.file("extdata", "ct_test.csv", package = "ddFoldChange")
##' output_file <- "processed_ct_test.csv"
##' df.process<- ProcessData(file.path = file_path, replicate = 3, output.file = output_file)
##' head(df.process)
##'
##' # Step 2: Plot the fold change for target genes using the processed data
##' PlotFC(df.process)
##' @export

PlotFC <- function(df.process) {


  targets = colnames(df.process)[5:((ncol(df.process) - 4)/4 + 4)]  #create a name list of gene interested


  for (target in targets) {
    fold.change.col <- paste0("FC.ddCt.dCt.", target)

    p <- ggplot2::ggplot(df.process, ggplot2::aes(x = trt, y = .data[[fold.change.col]], fill = trt)) + ggplot2::geom_bar(stat = "summary", fun = "mean", position = ggplot2::position_dodge(),
                                                                                                                          color = "black", width = 0.6) + ggplot2::labs(title = paste("Target", target, "gene expression"), x = "Group", y = paste("Fold Change", "(", target, ")")) +
      ggplot2::theme_minimal() + ggplot2::theme(plot.title = element_text(hjust = 0.5))  # Center the title

    print(p)
  }
}
