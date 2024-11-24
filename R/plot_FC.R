##' Plot Fold Change for Target Genes
##' 
##' We want to use ggplot to graph the fold change for each of the target genes.
##' @title Bar Graph for Gene Fold Change
##' @author Xia Zhao
##' @param df.process The output data created from the `process.data` function, including the fold change
##' 
##'  \describe{
##'     \item{trt}{A factor column indicating treatment groups.}
##'     \item{Fold Change Columns}{Fold change values for each target gene (e.g., `FC.<gene>`).}
#'   }
##' @import ggplot2
##' @importFrom utils read.csv
##' @export

plot_FC <- function(df.process) {


    targets = colnames(df.process)[5:((ncol(df.process) - 4)/4 + 4)]  #create a name list of gene interested


    for (target in targets) {
        fold.change.col <- paste0("FC.ddCt.dCt.", target)

        p <- ggplot2::ggplot(df.process, ggplot2::aes(x = trt, y = .data[[fold.change.col]], fill = trt)) + ggplot2::geom_bar(stat = "summary", fun = "mean", position = ggplot2::position_dodge(),
            color = "black", width = 0.6) + ggplot2::labs(title = paste("Target", target, "gene expression"), x = "Group", y = paste("Fold Change", "(", target, ")")) +
            ggplot2::theme_minimal() + ggplot2::theme(plot.title = element_text(hjust = 0.5))  # Center the title

        print(p)
    }
}
