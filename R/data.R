#' test BRIC-seq dataset for RNA half-life comparison
#'
#' A dataset containing the RPKM for each time point and information column
#'  about 200 genes. The variables are as follows:
#'
#' @format A data frame with 200 rows and 20 variables:
#' \describe{
#'   \item{gr_id}{Group id}
#'   \item{symbol}{Gene symbol}
#'   \item{accession_id}{Gene accession id (RefSeq)}
#'   \item{locus}{Genome locus}
#'   \item{CTRL_1_0h}{RPKM value at 0h in control condition}
#'   \item{CTRL_1_1h}{RPKM value at 1h in control condition}
#'   \item{CTRL_1_2h}{RPKM value at 2h in control condition}
#'   \item{CTRL_1_4h}{RPKM value at 4h in control condition}
#'   \item{CTRL_1_8h}{RPKM value at 8h in control condition}
#'   \item{CTRL_1_12h}{RPKM value at 12h in control condition}
#'   \item{gr_id}{Group id}
#'   \item{symbol}{Gene symbol}
#'   \item{accession_id}{Gene accession id (RefSeq)}
#'   \item{locus}{Genome locus}
#'   \item{KD_1_0h}{RPKM value at 0h in knockdown condition}
#'   \item{KD_1_1h}{RPKM value at 1h in knockdown condition}
#'   \item{KD_1_2h}{RPKM value at 2h in knockdown condition}
#'   \item{KD_1_4h}{RPKM value at 4h in knockdown condition}
#'   \item{KD_1_8h}{RPKM value at 8h in knockdown condition}
#'   \item{KD_1_12h}{RPKM value at 12h in knockdown condition}
#' }
"RNA_halflife_comparison"


#' test BRIC-seq dataset for p-value estimation using grubbs test
#'
#' A dataset containing the RPKM for each time point and information column
#'  about 200 genes. The variables are as follows:
#'
#' @format A data frame with 200 rows and 40 variables:
#' \describe{
#'   \item{gr_id}{Group id}
#'   \item{symbol}{Gene symbol}
#'   \item{accession_id}{Gene accession id (RefSeq)}
#'   \item{locus}{Genome locus}
#'   \item{CTRL_1_0h}{RPKM value at 0h in control condition}
#'   \item{CTRL_1_1h}{RPKM value at 1h in control condition}
#'   \item{CTRL_1_2h}{RPKM value at 2h in control condition}
#'   \item{CTRL_1_4h}{RPKM value at 4h in control condition}
#'   \item{CTRL_1_8h}{RPKM value at 8h in control condition}
#'   \item{CTRL_1_12h}{RPKM value at 12h in control condition}
#'   \item{gr_id}{Group id}
#'   \item{symbol}{Gene symbol}
#'   \item{accession_id}{Gene accession id (RefSeq)}
#'   \item{locus}{Genome locus}
#'   \item{CTRL_2_0h}{RPKM value at 0h in control condition}
#'   \item{CTRL_2_1h}{RPKM value at 1h in control condition}
#'   \item{CTRL_2_2h}{RPKM value at 2h in control condition}
#'   \item{CTRL_2_4h}{RPKM value at 4h in control condition}
#'   \item{CTRL_2_8h}{RPKM value at 8h in control condition}
#'   \item{CTRL_2_12h}{RPKM value at 12h in control condition}
#'   \item{gr_id}{Group id}
#'   \item{symbol}{Gene symbol}
#'   \item{accession_id}{Gene accession id (RefSeq)}
#'   \item{locus}{Genome locus}
#'   \item{CTRL_3_0h}{RPKM value at 0h in control condition}
#'   \item{CTRL_3_1h}{RPKM value at 1h in control condition}
#'   \item{CTRL_3_2h}{RPKM value at 2h in control condition}
#'   \item{CTRL_3_4h}{RPKM value at 4h in control condition}
#'   \item{CTRL_3_8h}{RPKM value at 8h in control condition}
#'   \item{CTRL_3_12h}{RPKM value at 12h in control condition}
#'   \item{gr_id}{Group id}
#'   \item{symbol}{Gene symbol}
#'   \item{accession_id}{Gene accession id (RefSeq)}
#'   \item{locus}{Genome locus}
#'   \item{KD_1_0h}{RPKM value at 0h in knockdown condition}
#'   \item{KD_1_1h}{RPKM value at 1h in knockdown condition}
#'   \item{KD_1_2h}{RPKM value at 2h in knockdown condition}
#'   \item{KD_1_4h}{RPKM value at 4h in knockdown condition}
#'   \item{KD_1_8h}{RPKM value at 8h in knockdown condition}
#'   \item{KD_1_12h}{RPKM value at 12h in knockdown condition}
#' }
"RNA_halflife_grubbs_test"
