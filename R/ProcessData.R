##' Process the input Data for Fold Change Calculations
##'
##' We want to perform normalization of Ct values for the target genes and calculate fold-change calculations relative to the control group
##' @title Gene fold change calculation
##' @author Xia Zhao
##' @param file.path Path to the input CSV file
##' @param replicate Number of replicates for each group
##' @param output.file Created Output file for fold change plotting in the next function
##' \describe{
##'   \item{means.internal}{Mean Ct values of the housekeeping gene for each group}
##'   \item{dCt.col.stat}{The column corresponding to the first target gene}
##'   \item{columns.to.normalize}{The columns corresponding to all target genes}
##'   \item{baseline.means}{The mean delta CT value for the control group for each target gene}
##'   \item{new.col.name}{Create new columns for delta CT, delta CT, and fold change in each loop}
##'   }
##' @return A data frame containing original Ct value, delta CT value, delta delta CT value and fold changes
##'   \describe{
##'     \item{Ct values}{Original Ct values from the input file.}
##'     \item{dCt}{Delta Ct values for each target gene.}
##'     \item{ddCt}{Delta delta Ct values for each target gene.}
##'     \item{FC}{Calculated fold changes for each target gene.}
##'   }
##' @importFrom utils read.csv write.csv
##' @export
##'
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
ProcessData <- function(file.path, replicate, output.file) {

  # Add a check to make sure the file exists
  if (!file.exists(file.path)) {
    stop("File not found: ", file.path)
  }

  # load the data
  df <- read.csv(file.path, header = TRUE, sep = ",")

  # Calculate the mean for every replicate rows in column 3 (the internal control) and add to new column 'means_internal'
  means.internal <- rep(sapply(seq(1, nrow(df), by = replicate), function(i) {
    mean(df[i:min(i + (replicate - 1), nrow(df)), 3])
  }), each = replicate)

  # combine the internal means and other column in the data, treating it as reference
  df <- cbind(df[, 1:2], means.internal, df[, 3:ncol(df)])
  reference.column <- "means.internal"

  # Skip columns 1-4, Calculate delta Ct value for each target column
  for (col in names(df)) {
    if (!col %in% colnames(df)[1:4]) {
      new.col.name <- paste("dCt", col, sep = ".")
      df[[new.col.name]] <- df[[col]] - df[[reference.column]]
    }
  }

  dCt.col.stat <- (ncol(df) - 4)/2 + 5
  columns.to.normalize <- names(df)[dCt.col.stat:ncol(df)]
  baseline.means <- colMeans(df[1:replicate, columns.to.normalize])

  for (col in columns.to.normalize) {
    new.col.name <- paste("ddCt", col, sep = ".")
    df[[new.col.name]] <- df[[col]] - baseline.means[col]
  }

  fold.col.stat <- (ncol(df) - 4)/3 * 2 + 5
  columns.for.fold <- names(df)[fold.col.stat:ncol(df)]

  for (col in columns.for.fold) {
    new.col.name <- paste("FC", col, sep = ".")
    df[[new.col.name]] <- 2^(-df[[col]])
  }

  # Save the processed data to CSV in inst/extdata/
  output_path <- file.path(system.file("extdata", package = "ddFoldChange"), output.file)
  write.csv(df, file = output_path, row.names = FALSE)

  df.process<-df
  return(df.process)
}
