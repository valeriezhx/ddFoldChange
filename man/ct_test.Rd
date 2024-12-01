
\name{ct_test}
\docType{data}
\alias{ct_test}
\title{Example Ct Data for qPCR Analysis}
\description{
  This dataset contains Ct values from qPCR experiments used to calculate fold changes for target genes.
  The dataset includes Ct values for both target genes and an internal control gene (housekeeping gene).
}
\format{
  A data frame with 9 rows and 5 columns:
  \describe{
    \item{ID}{The sample name or ID.}
    \item{trt}{The treatment group (e.g., control, treated).}
    \item{x}{Ct values for the housekeeping gene (internal control).}
    \item{y}{Ct values for the first target gene x.}
    \item{z}{Ct values for the second target gene y.}
  }
}
\examples{
  # Load the dataset
  data(ct_test)
  
  # View the first few rows of the dataset
  head(ct_test)
}

