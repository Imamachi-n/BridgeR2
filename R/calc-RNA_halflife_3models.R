#' Calculate RNA half-life for each gene.
#'
#' \code{BridgeRHalfLifeCalc3models} returns the dataframe of
#' RNA half-life for each gene.
#'
#' @param inputFile The vector of tab-delimited matrix file.
#'
#' @param group The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param inforColumn The number of information columns.
#'
#' @param CutoffTimePointNumber The number of minimum time points for calc.
#'
#' @param save Whether to save the output matrix file.
#'
#' @param outputPrefix The prefix for the name of the output.
#'
#' @export
#'
#' @import data.table
#' @importFrom stats optim
#' @importFrom stats cor

BridgeRHalfLifeCalc3models <- function(inputFile,
                                       group = c("Control","Knockdown"),
                                       hour = c(0, 1, 2, 4, 8, 12),
                                       inforColumn = 4,
                                       CutoffTimePointNumber = 4,
                                       save = T,
                                       outputPrefix = "BridgeR_5"){

  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.numeric(CutoffTimePointNumber))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.logical(save))
  stopifnot(is.character(outputPrefix))

  # file infor prep
  time_points <- length(hour)
  group_number <- length(group)
  input_matrix <- inputFile

  # halflife infor header prep
  halflife_infor_header <- c("Model", "R2", "half_life")

  # half calc function


  half_calc_3model <- function(time_exp_table){
    data_point <- length(time_exp_table$exp)
    if(!is.null(time_exp_table)){
      if(data_point >= CutoffTimePointNumber){
        #print(time_exp_table$exp[1])
        if(as.numeric(as.vector(as.matrix(time_exp_table$exp[1]))) > 0){
          # model1
          optim1 <- function(x){
            mRNA_exp <- exp(-x * time_exp_table$hour)
            sum((time_exp_table$exp - mRNA_exp)^2)
          }

          model1_pred <- function(x){
            mRNA_exp <- exp(-x * time_exp_table$hour)
            (cor(mRNA_exp, time_exp_table$exp, method="pearson"))^2
          }

          model1_half <- function(x){
            mRNA_half <- exp(-a_1 * x)
            (mRNA_half - 0.5)^2
          }

          # model2
          optim2 <- function(x){
            mRNA_exp <- (1.0 - x[2]) * exp(-x[1] * time_exp_table$hour) + x[2]
            sum((time_exp_table$exp - mRNA_exp)^2)
          }

          model2_pred <- function(a,b){
            mRNA_exp <- (1.0 - b) * exp(-a * time_exp_table$hour) + b
            (cor(mRNA_exp, time_exp_table$exp, method="pearson"))^2
          }

          model2_half <- function(x){
            mRNA_half <- (1.0 - b_2) * exp(-a_2 * x) + b_2
            (mRNA_half - 0.5)^2
          }

          # model3
          optim3 <- function(x){
            mRNA_exp <- x[3] * exp(-x[1] * time_exp_table$hour) + (1.0 - x[3]) * exp(-x[2] * time_exp_table$hour)
            sum((time_exp_table$exp - mRNA_exp)^2)
          }

          model3_pred <- function(a,b,c){
            mRNA_exp <- c * exp(-a * time_exp_table$hour) + (1.0 - c) * exp(- b * time_exp_table$hour)
            (cor(mRNA_exp, time_exp_table$exp, method="pearson"))^2
          }

          model3_half <- function(x){
            mRNA_half <- c_3 * exp(-a_3 * x) + (1.0 - c_3) * exp(-b_3 * x)
            (mRNA_half - 0.5)^2
          }

          # Model1: f(t) = exp(-a * t)
          out1 <- suppressWarnings(optim(1,optim1))
          min1 <- out1$value
          a_1 <- out1$par[1]
          half1_1 <- log(2) / a_1
          cor1 <- model1_pred(a_1)
          out1 <- suppressWarnings(optim(1,model1_half))
          half1_2 <- out1$par
          mRNA_pred <- exp(-a_1 * time_exp_table$hour)
          s2 <- sum((time_exp_table$exp - mRNA_pred)^2) / data_point
          AIC1 <- data_point * log(s2) + (2 * 0)

          # Model2: f(t) = (1 - b)exp(-a * t) + b
          out2 <- suppressWarnings(optim(c(1,0),optim2))
          min2 <- out2$value
          a_2 <- out2$par[1]
          b_2 <- out2$par[2]
          cor2 <- model2_pred(a_2, b_2)
          out2 <- suppressWarnings(optim(1,model2_half))
          half2 <- out2$par
          if(b_2 >= 0.5){
            half2 <- Inf
          }
          mRNA_pred <- (1.0 - b_2) * exp(-a_2 * time_exp_table$hour) + b_2
          s2 <- sum((time_exp_table$exp - mRNA_pred)^2) / data_point
          AIC2 <- data_point * log(s2) + (2 * 1)

          ###Model3: f(t) = c * exp(-a * t) + (1 - c) * exp(-b * t)###
          out3 <- suppressWarnings(optim(c(1,1,0.1),optim3))
          min3 <- out3$value
          a_3 <- out3$par[1]
          b_3 <- out3$par[2]
          c_3 <- out3$par[3]
          cor3 <- model3_pred(a_3,b_3,c_3)
          out3 <- suppressWarnings(optim(1,model3_half))
          half3 <- out3$par
          mRNA_pred <- c_3 * exp(-a_3 * time_exp_table$hour) + (1.0 - c_3) * exp(-b_3 * time_exp_table$hour)
          s2 <- sum((time_exp_table$exp - mRNA_pred)^2) / data_point
          AIC3 <- data_point * log(s2) + (2 * 2)

          # AIC selection
          half_table <- data.frame(half=c(as.numeric(half1_2),
                                          as.numeric(half2),
                                          as.numeric(half3)),
                                   AIC=c(as.numeric(AIC1),
                                         as.numeric(AIC2),
                                         as.numeric(AIC3)))
          AIC_list <- order(half_table$AIC)

          AIC_flg <- 0
          model <- NULL
          R2 <- NULL
          halflife <- NULL
          for(min_AIC_index in AIC_list){
            # Model
            if(min_AIC_index == 1){
              selected_half <- half1_2
              if(selected_half > 24){
                selected_half <- 24
              }
              if(a_1 > 0){
                model <- "model1"
                R2 <- cor1
                halflife <- selected_half
                AIC_flg <- 1
                break
              }
            }
            # Model2
            if(min_AIC_index == 2){
              selected_half <- half2
              if(selected_half == "Inf"){
                selected_half <- 24
              }else if(selected_half > 24){
                selected_half <- 24
              }
              if(a_2 > 0 && b_2 > 0 && b_2 < 1){
                model <- "model2"
                R2 <- cor2
                halflife <- selected_half
                AIC_flg <- 1
                break
              }
            }
            # Model3
            if(min_AIC_index == 3){
              selected_half <- half3
              if(selected_half > 24){
                selected_half <- 24
              }
              if(a_3 > 0 && b_3 > 0 && c_3 > 0 && c_3 < 1){
                model <- "model3"
                R2 <- cor3
                halflife <- selected_half
                AIC_flg <- 1
                break
              }
            }
          }

          if(AIC_flg == 0){
            if(a_1 < 0){
              model <- "model1"
              R2 <- cor1
              halflife <- 24
            }else{
              model <- "no_good_model"
              R2 <- "NA"
              halflife <- 24
            }
          }

          return(c(model,
                   R2,
                   halflife))
        }else{
          return(rep("NA", 3))
        }
      }else{
        return(rep("NA", 3))
      }
    }else{
      return(rep("NA", 3))
    }
  }

  # main function
  halflife_calc_3model <- function(data){
    data_vector <- NULL
    # infor data
    gene_infor <- data[infor_st:infor_ed]

    # exp data
    exp <- data[exp_st:exp_ed]
    #print(exp)
    data_vector <- c(gene_infor, exp)

    # outlier
    if (all(is.nan(exp)) || all(exp == "NaN")) {
      result_list <- rep("NA", 3)
      data_vector <- c(data_vector, result_list)
      #print("NG")
      return(data_vector)
    }

    # exp data as numeric
    exp <- as.numeric(exp)

    # raw data prep
    time_point_exp_raw <- data.frame(hour,exp)

    # Cutoff RPKM: 0
    time_point_exp_base <- time_point_exp_raw[time_point_exp_raw$exp > 0, ]

    # result
    model3_result <- half_calc_3model(time_point_exp_base)
    data_vector <- c(data_vector, model3_result)
    return(data_vector)
  }

  output_matrix <- NULL
  for (group_index in 1:group_number) {
    # header prep
    colname_st <- 1 + (group_index - 1) * (inforColumn + time_points)
    colname_ed <- group_index * (inforColumn + time_points)
    header_label <- colnames(input_matrix)[colname_st:colname_ed]

    # Infor data
    infor_st_ed <- generate_infor_st_ed(group_index,
                                        time_points,
                                        inforColumn)
    infor_st <- infor_st_ed[1]
    infor_ed <- infor_st_ed[2]

    # Exp data
    exp_st <- infor_ed + 1
    exp_ed <- infor_ed + time_points

    # output result
    result_matrix <- t(apply((input_matrix), 1, halflife_calc_3model))
    colnames(result_matrix) <- c(header_label, halflife_infor_header)
    output_matrix <- cbind(output_matrix, result_matrix)
  }
  output_matrix <- data.table(output_matrix)

  if (save == T) {
    write.table(output_matrix, quote = F, sep = "\t", row.names = F,
                file = paste(outputPrefix, "_halflife_calc_3models.txt", sep=""))
  }

  return(output_matrix)
}
