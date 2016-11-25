#' Calculate Fold-change of RNA half-life and p-value.
#'
#' \code{BridgeRPvalueEvaluation} returns the dataframe of
#' Fold-change of RNA half-life and p-value.
#'
#' @param inputFile The vector of tab-delimited matrix file.
#'
#' @param group The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param comparisonFile The vector of group names.
#'
#' @param inforColumn The number of information columns.
#'
#' @param CutoffTimePointNumber The number of minimum time points for calc.
#'
#' @param calibration Calibration of RNA half-life.
#'
#' @param save Whether to save the output matrix file.
#'
#' @param outputPrefix The prefix for the name of the output.
#'
#' @export
#'
#' @import data.table
#'
#' @importFrom BSDA tsum.test
#' @importFrom stats predict

BridgeRPvalueEvaluation <- function(inputFile,
                                    group = c("Control","Knockdown"),
                                    hour = c(0, 1, 2, 4, 8, 12),
                                    comparisonFile = c("Control","Knockdown"),
                                    inforColumn = 4,
                                    CutoffTimePointNumber = 4,
                                    calibration = FALSE,
                                    save = TRUE,
                                    outputPrefix = "BridgeR_6"){

  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.character(comparisonFile) && is.vector(comparisonFile))
  stopifnot(is.numeric(CutoffTimePointNumber))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.logical(calibration))
  stopifnot(is.logical(save))
  stopifnot(is.character(outputPrefix))

  # file infor prep
  time_points <- length(hour)
  group_number <- length(group)
  input_matrix <- inputFile
  comp_file_index <- as.vector(sapply(comparisonFile,
                                       function(test) which(group == test)))
  comp_label <- paste("log2(Relative half-life[",
                      comparisonFile[2], "/", comparisonFile[1],"])", sep="")

  # header/calibration rate prep
  halflife_table <- NULL
  correction_value <- NULL
  header_label <- NULL
  pvalue_label <- c(comp_label,
                    "p-value(Welch Modified Two-Sample t-test)")

  for (group_index in comp_file_index) {
    # header prep, +3: Model, R2, half_life
    colname_st <- 1 + (group_index - 1) * (inforColumn + time_points + 3)
    colname_ed <- group_index * (inforColumn + time_points + 3)
    header_label <- append(header_label,
                           colnames(input_matrix)[colname_st:colname_ed])

    # correction value prep
    if (calibration == TRUE) {
      # infor_st <- 1 + (group_index - 1) * (inforColumn + time_points + 3)
      infor_ed <- inforColumn + (group_index - 1) * (inforColumn + time_points + 3)
      halflife_index <- infor_ed + time_points + 3
      halflife_list <- as.numeric(input_matrix[[halflife_index]])
      halflife_table <- cbind(halflife_table, halflife_list)
    }
  }

  # calc correction value
  if (calibration == TRUE) {
    colnames(halflife_table) <- c("x", "y")
    halflife_table <- as.data.frame(halflife_table)
    halflife_table <- halflife_table[halflife_table$x < 24,]
    halflife_table <- halflife_table[halflife_table$y < 24,]
    test_lm <- lm(y ~ x + 0, data = halflife_table)
    correction_value <- round(as.numeric(test_lm$coefficients), digits = 3)
    # print(paste("Correction value: ", correction_value, sep = ""))
  }

  # calc p-value
  # function1
  halflife_upr_lwr_calc <- function(dataframe,
                                    predicted_exp_ul,
                                    halflife_value){
    predicted_exp_ul_table <-
      data.frame(hour=dataframe$hour, exp=predicted_exp_ul)
    predicted_exp_ul_model <-
      lm(predicted_exp_ul_table$exp ~ predicted_exp_ul_table$hour)
    predicted_exp_ul_model_summary <-
      summary(predicted_exp_ul_model)
    predicted_exp_ul_coef <-
      predicted_exp_ul_model_summary$coefficients[2]
    predicted_exp_ul_intercept <-
      predicted_exp_ul_model_summary$coefficients[1]
    predicted_exp_ul_R2 <-
      predicted_exp_ul_model_summary$r.squared #
    predicted_exp_ul_halflife <-
      (halflife_value - predicted_exp_ul_intercept) / predicted_exp_ul_coef
    return(predicted_exp_ul_halflife)
  }

  # function2
  pvalue_calc <- function(dataframe, flg){
    data_point <- length(dataframe$exp)
    # print(dataframe)
    if (!is.null(dataframe)) {
      if (data_point >= CutoffTimePointNumber) {
        if (as.numeric(as.vector(as.matrix(dataframe$exp[1]))) > 0) {
          # calc RNA half-life and R2
          fitted_model <- lm(log(dataframe$exp) ~ dataframe$hour - 1)
          model_summary <- summary(fitted_model)
          coef <- -(model_summary$coefficients[1])
          coef_error <- model_summary$coefficients[2]
          coef_p <- model_summary$coefficients[4]
          r_squared <- model_summary$r.squared
          adj_r_squared <- model_summary$adj.r.squared
          residual_standard_err <- model_summary$sigma

          half_life <- log(2)/coef
          half_life_w <- half_life
          halflife_value <- log(0.5)

          # half-life calibration
          if (calibration == TRUE) {
            if (flg == 1) {
              half_life_w <- half_life / correction_value
              half_life <- half_life / correction_value
              halflife_value <- -coef * half_life
            }
          }

          # >=24 hr: 24hr
          if (coef < 0) {
            half_life <- Inf
            half_life_w <- 24
          } else if (half_life > 24) {
            half_life_w <- 24
          }

          # predict fitting curve from time points
          predict_data <- data.frame(hour = dataframe$hour)
          pred_conf <- predict(fitted_model, predict_data, se.fit = TRUE)

          # fitting SD and value SD
          predicted_exp <- pred_conf$fit
          predicted_exp_SE <- pred_conf$se.fit    # Residual for model ??
          predicted_exp_df <- pred_conf$df    # Degree of freedom
          predicted_exp_residual <-
            pred_conf$residual.scale    # Residual for exp data ??

          # residual <- sqrt(predicted_exp_SE^2 + predicted_exp_residual^2) *
          #                    qt(0.975,df) #95% Prediction interval
          SE_residual <- sqrt(predicted_exp_SE^2 +
                                predicted_exp_residual^2)    # SE
          SD_residual <- SE_residual * sqrt(predicted_exp_df)
          predicted_exp_lwr <- predicted_exp - SD_residual
          predicted_exp_upr <- predicted_exp + SD_residual

          # predicted lwr/upr RNA half-life
          predicted_exp_lwr_halflife <- halflife_upr_lwr_calc(dataframe,
                                                              predicted_exp_lwr,
                                                              halflife_value)
          predicted_exp_upr_halflife <- halflife_upr_lwr_calc(dataframe,
                                                              predicted_exp_upr,
                                                              halflife_value)

          halflife_SD_minus <- half_life - predicted_exp_lwr_halflife
          halflife_SD_plus <- predicted_exp_upr_halflife - half_life

          #Half-life, The degree of freedom, SD_minus, SD_plus
          return(c(half_life_w, half_life, predicted_exp_df,
                   halflife_SD_minus, halflife_SD_plus))
        }
      }
    }
  }

  # function3
  pvalue_estimation <- function(data){
    data_vector <- NULL

    # data infor prep
    flg <- 0
    halflife_w <- NULL
    mean_for_p <- NULL
    N_for_p <- NULL
    SD_minus_for_p <- NULL
    SD_plus_for_p <- NULL
    flg_na <- 0

    # check each sample
    for (group_index in comp_file_index) {
      # infor data prep
      infor_st_ed <- generate_infor_st_ed(group_index,
                                          time_points,
                                          inforColumn)
      infor_st <- 1 + (group_index - 1) * (inforColumn + time_points + 3)
      infor_ed <- inforColumn + (group_index - 1) * (inforColumn + time_points + 3)
      gene_infor <- data[infor_st:infor_ed]

      # exp data prep
      exp_st <- infor_ed + 1
      exp_ed <- infor_ed + time_points
      exp <- as.numeric(data[exp_st:exp_ed])
      time_point_exp_raw <- data.frame(hour,exp)

      # output prep
      colname_st <- 1 + (group_index - 1) * (inforColumn + time_points + 3)
      colname_ed <- group_index * (inforColumn + time_points + 3)
      data_vector <- append(data_vector, data[colname_st:colname_ed])

      # model label prep
      model_index <- exp_ed + 1 #1: Model, #2:R2, #3:half_life
      model <- data[model_index]
      # print(model)

      # check exp data
      ttest_infor <- NULL
      if (model == "NA" || is.na(model)) {
        flg_na <- 1
        flg <- 1
        next
      } else if (model == "Raw"){
        time_point_exp_base <- time_point_exp_raw[time_point_exp_raw$exp > 0, ]
        ttest_infor <- pvalue_calc(time_point_exp_base, flg)
      } else {
        check <- parse_rm_hr_infor(model)    # util.R
        check_index <- sapply(check,
                              function(t) which(time_point_exp_raw$hour == t))
        time_point_exp_del <- time_point_exp_raw[-check_index,]
        time_point_exp_del <- time_point_exp_del[time_point_exp_del$exp > 0,]
        ttest_infor <- pvalue_calc(time_point_exp_del, flg)
      }

      # result
      halflife_w <- append(halflife_w, ttest_infor[1])
      mean_for_p <- append(mean_for_p, ttest_infor[2])
      N_for_p <- append(N_for_p, ttest_infor[3])
      SD_minus_for_p <- append(SD_minus_for_p, ttest_infor[4])
      SD_plus_for_p <- append(SD_plus_for_p, ttest_infor[5])
      flg <- 1
    }

    # pvalue estimation
    halflife_comp <- "NA"
    p_value <- "NA"

    if (calibration == TRUE){
      new_halflife_w_cond2 <- NULL
      if (flg_na == 1) {
        new_halflife_w_cond2 <- "NA"
      } else {
        new_halflife_w_cond2 <- halflife_w[2]
      }
      data_vector[length(data_vector)] <- new_halflife_w_cond2
    }

    if (flg_na == 1){
      data_vector <- append(data_vector, rep("NA", 2))
      return(data_vector)
    }

    if (mean_for_p[1] <= mean_for_p[2]) {    # Control_T1/2 <= KD_T1/2
      cond1_halflife_mean <- mean_for_p[1]
      cond2_halflife_mean <- mean_for_p[2]
      cond1_df <- N_for_p[1]
      cond2_df <- N_for_p[2]
      cond1_halflife_SD <- SD_plus_for_p[1]    # check!
      cond2_halflife_SD <- SD_minus_for_p[2]    # check!

      t_test <- tsum.test(mean.x = cond1_halflife_mean,
                          n.x = cond1_df,
                          s.x = cond1_halflife_SD,
                          mean.y = cond2_halflife_mean,
                          n.y = cond2_df,
                          s.y = cond2_halflife_SD)
      p_value <- t_test$p.value
    } else if (mean_for_p[1] > mean_for_p[2]){    # Control_T1/2 > KD_T1/2
      cond1_halflife_mean <- mean_for_p[1]
      cond2_halflife_mean <- mean_for_p[2]
      cond1_df <- N_for_p[1]
      cond2_df <- N_for_p[2]
      cond1_halflife_SD <- SD_minus_for_p[1]    # check!
      cond2_halflife_SD <- SD_plus_for_p[2]    # check!

      t_test <- tsum.test(mean.x=cond1_halflife_mean,
                          n.x=cond1_df,
                          s.x=cond1_halflife_SD,
                          mean.y=cond2_halflife_mean,
                          n.y=cond2_df,
                          s.y=cond2_halflife_SD)
      p_value <- t_test$p.value
    }

    # Fold-change (RNA half-lives: KD_T_half / Control_T_half )
    halflife_comp <- log2(halflife_w[2] / halflife_w[1])

    data_vector <- append(data_vector, c(halflife_comp, p_value))
    return(data_vector)
  }

  # result
  output_matrix <- t(apply((input_matrix), 1, pvalue_estimation))
  colnames(output_matrix) <- c(header_label, pvalue_label)
  output_matrix <- data.table(output_matrix)

  if (save == T) {
    write.table(output_matrix, quote = F, sep = "\t", row.names = F,
                file = paste(outputPrefix, "_halflife_pvalue_evaluation.txt", sep=""))
  }

  return(output_matrix)

}

# testing
# library(BSDA)
# pvalue_table <- BridgeRPvalueEvaluation(halflife_table, calibration = TRUE)
