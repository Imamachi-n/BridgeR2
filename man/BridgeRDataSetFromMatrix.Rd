% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc-Relative_RPKM.R
\name{BridgeRDataSetFromMatrix}
\alias{BridgeRDataSetFromMatrix}
\title{Calculate relative RPKM expression from data.table format.}
\usage{
BridgeRDataSetFromMatrix(inputFile, group = c("Control", "Knockdown"),
  hour = c(0, 1, 2, 4, 8, 12), cutoff = 0.1, cutoffBelow = 0.1,
  inforColumn = 4, save = T, outputPrefix = "BridgeR_1")
}
\arguments{
\item{inputFile}{The vector of tab-delimited matrix file.}

\item{group}{The vector of group names.}

\item{hour}{The vector of time course about BRIC-seq experiment.}

\item{cutoff}{Cutoff value of RPKM at 0hr.}

\item{cutoffBelow}{Cutoff value of RPKM at all time points.}

\item{inforColumn}{The number of information columns.}

\item{save}{Whether to save the output matrix file.}

\item{outputPrefix}{The prefix for the name of the output.}
}
\value{
data.table object about relative RPKM values.
}
\description{
\code{BridgeRDataSetFromMatrix} calculates the relative RPKM values compared with 0hr,
importing data.table format.
}
\examples{
library(data.table)
rpkm_matrix <- data.table(gr_id = c(8, 9, 14),
                          symbol = c("AAAS", "AACS", "AADAT"),
                          accession_id = c("NM_015665", "NM_023928", "NM_182662"),
                          locus = c("chr12", "chr12", "chr4"),
                          CTRL_1_0h = c(41, 5, 5),
                          CTRL_1_1h = c(48, 7, 6),
                          CTRL_1_2h = c(56, 10, 6),
                          CTRL_1_4h = c(87, 12, 10),
                          CTRL_1_8h = c(124, 20, 11),
                          CTRL_1_12h = c(185, 22, 15),
                          gr_id = c(8, 9, 14),
                          symbol = c("AAAS", "AACS", "AADAT"),
                          accession_id = c("NM_015665", "NM_023928", "NM_182662"),
                          locus = c("chr12", "chr12", "chr4"),
                          KD_1_0h = c(21, 10, 3),
                          KD_1_1h = c(33, 11, 3),
                          KD_1_2h = c(42, 15, 4),
                          KD_1_4h = c(60, 20, 5),
                          KD_1_8h = c(65, 37, 6),
                          KD_1_12h = c(70, 42, 6))
group <- c("Control", "Knockdown")
hour <- c(0, 1, 2, 4, 8, 12)
test_table <- BridgeRDataSetFromMatrix(inputFile = rpkm_matrix,
                                       group = group,
                                       hour = hour,
                                       cutoff = 0.1,
                                       inforColumn = 4,
                                       save = FALSE)

}

