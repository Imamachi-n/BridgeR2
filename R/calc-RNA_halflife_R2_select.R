#' Calculate RNA half-life for each gene.
#'
#' \code{BridgeRHalfLifeCalcR2Select} returns the dataframe of
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
#' @param R2_criteria The cutoff of R2 for R2 selection.
#'
#' @param TimePointRemoval1 The candicate_1 of time point removal.
#'
#' @param TimePointRemoval2 The candicate_2 of time point removal.
#'
#' @param ThresholdHalfLife1 The cutoff of TimePointRemoval1.
#'
#' @param ThresholdHalfLife2 The cutoff of TimePointRemoval2.
#'
#' @param save Whether to save the output matrix file.
#'
#' @param outputPrefix The prefix for the name of the output.

BridgeRHalfLifeCalcR2Select <- function(inputFile,
                                        group = c("Control","Knockdown"),
                                        hour = c(0, 1, 2, 4, 8, 12),
                                        inforColumn = 4,
                                        CutoffTimePointNumber = 4,
                                        R2_criteria = 0.90, # 0.90,
                                        TimePointRemoval1 = c(1,2),
                                        TimePointRemoval2 = c(8,12),
                                        ThresholdHalfLife1 = 3,
                                        ThresholdHalfLife2 = 12,
                                        save = T,
                                        outputPrefix = "BridgeR_5"){

  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.numeric(CutoffTimePointNumber))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.numeric(TimePointRemoval1) && is.vector(TimePointRemoval1))
  stopifnot(is.numeric(TimePointRemoval2) && is.vector(TimePointRemoval2))
  stopifnot(is.numeric(ThresholdHalfLife1))
  stopifnot(is.numeric(ThresholdHalfLife2))
  stopifnot(is.logical(save))
  stopifnot(is.character(outputPrefix))

  # file infor prep
  time_points <- length(hour)
  group_number <- length(group)
  input_matrix <- inputFile

  # timepoint removal prep
  TimePointRemoval1 <- sort(TimePointRemoval1, decreasing = T)
  TimePointRemoval2 <- sort(TimePointRemoval2, decreasing = F)
  TimePointRemoval1_length <- length(TimePointRemoval1)
  TimePointRemoval2_length <- length(TimePointRemoval2)
  # TimePointRemoval_combination <- choose(TimePointRemoval1_length, 2) +
  #                                 choose(TimePointRemoval2_length, 2)

  # halflife infor header prep
  halflife_infor_header <- c("Model", "R2", "half_life")

  # halflife calc
  ### Function1
  half_calc <- function(time_exp_table, label){
    data_point <- length(time_exp_table$exp)
    if(!is.null(time_exp_table)){
      if(data_point >= CutoffTimePointNumber){
        #print(time_exp_table$exp[1])
        if(as.numeric(as.vector(as.matrix(time_exp_table$exp[1]))) > 0){
          model <- lm(log(time_exp_table$exp) ~ time_exp_table$hour - 1)
          model_summary <- summary(model)
          coef <- -model_summary$coefficients[1]
          coef_error <- model_summary$coefficients[2]
          coef_p <- model_summary$coefficients[4]
          r_squared <- model_summary$r.squared
          adj_r_squared <- model_summary$adj.r.squared
          residual_standard_err <- model_summary$sigma
          half_life <- log(2) / coef
          if(coef < 0 || half_life >= 24){
            half_life <- 24
          }
          return(c(label,
                   r_squared,
                   half_life))
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

  # Function2
  test_R2 <- function(time_point_exp_raw, cutoff_data_point, halflife_Raw, R2_Raw){
    times_length <- length(cutoff_data_point) #c(12, 8) => 2, c(12) => 1
    times_index <- c(times_length)
    add_index <- times_length

    R2_list <- c(R2_Raw)
    half_list <- c(halflife_Raw)
    label_list <- c("Raw")
    for(counter in times_length:1){
      #excepted_time_points
      check_times <- cutoff_data_point[times_index] #c(24), c(24,12), c(24,12,8)
      time_point_exp_del <- NULL
      time_point_exp_del_label <- paste("Delete_",paste(check_times,collapse="hr_"),"hr",sep="")
      label_list <- append(label_list, time_point_exp_del_label) #
      time_point_exp_del <- time_point_exp_raw
      for(times_list in check_times){
        time_point_exp_del <- time_point_exp_del[time_point_exp_del$hour != as.numeric(times_list),]
      }
      time_point_exp_del <- time_point_exp_del[time_point_exp_del$exp > 0,]
      halflife_R2_result <- half_calc(time_point_exp_del, time_point_exp_del_label)

      R2_list <- append(R2_list, halflife_R2_result[2])    # r_squared
      half_list <- append(half_list, halflife_R2_result[3])    # half_life

      #Counter
      add_index <- add_index - 1
      times_index <- append(times_index, add_index)
    }

    R2_table <- data.frame(label=label_list, R2=R2_list, half=half_list)
    return(R2_table)
  }

  # Function3
  halflife_calc_R2_select <- function(data){
    data_vector <- NULL
    # infor data
    gene_infor <- data[infor_st:infor_ed]

    # exp data
    exp <- data[exp_st:exp_ed]
    #print(exp)
    data_vector <- c(gene_infor, exp)
    if (all(is.nan(exp)) || all(exp == "NaN")) {
      result_list <- rep("NA", 3)
      data_vector <- c(data_vector, result_list)
      #print("NG")
      return(data_vector)
    }
    exp <- as.numeric(exp)

    # raw data prep
    time_point_exp_raw <- data.frame(hour,exp)

    # Cutoff RPKM: 0
    time_point_exp_base <- time_point_exp_raw[time_point_exp_raw$exp > 0, ]

    # linear regression fitting
    LR_result <- half_calc(time_point_exp_base,"Exponential_Decay_Model")
    half_life_raw <- LR_result[3]
    R2_raw <- LR_result[2]

    # Re-calculation of RNA half-life
    R2_list <- NULL
    half_list <- NULL
    label_list <- NULL
    if (half_life_raw == "NA") {
      result_list <- rep("NA", 3)
      data_vector <- c(data_vector, result_list)
      return(data_vector)
    } else if (as.numeric(half_life_raw) < ThresholdHalfLife1) {
      R2_table <- test_R2(time_point_exp_raw,
                          TimePointRemoval2,
                          as.numeric(half_life_raw),
                          as.numeric(R2_raw))
      R2_table <- R2_table[R2_table$half != "NA",]
      sortlist <- order(as.numeric(as.vector(R2_table$R2)), decreasing = T)
      R2_table <- R2_table[sortlist,]
      # print(R2_table$half)
      # print(R2_table$half[1])
      # print(half_life_raw)
      # print((as.numeric(as.vector(R2_table$half[1])) > as.numeric(half_life_raw)))
      # print(R2_table)
      if (as.numeric(as.vector(R2_table$half[1])) > as.numeric(half_life_raw)) {
        data_vector <- c(data_vector, "Raw", R2_raw, half_life_raw)
        return(data_vector)
      }
    } else if (as.numeric(R2_raw) >= R2_criteria) {
      data_vector <- c(data_vector, "Raw", R2_raw, half_life_raw)
      return(data_vector)
    } else if (as.numeric(half_life_raw) >= ThresholdHalfLife2){
      R2_table <- test_R2(time_point_exp_raw,
                          TimePointRemoval1,
                          as.numeric(half_life_raw),
                          as.numeric(R2_raw))
    } else {
      data_vector <- c(data_vector, "Raw", R2_raw, half_life_raw)
      return(data_vector)
    }

    # R2 selection
    if (half_life_raw != "NA") {
      R2_table <- R2_table[R2_table$R2 != "NA",]
      sortlist <- order(as.numeric(as.vector(R2_table$R2)), decreasing = T)
      R2_table <- R2_table[sortlist,]
      result <- as.vector(as.matrix(R2_table[1,]))
      data_vector <- c(data_vector, result)
      return(data_vector)
    }
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
    result_matrix <- t(apply((input_matrix), 1, halflife_calc_R2_select))
    colnames(result_matrix) <- c(header_label, halflife_infor_header)
    output_matrix <- cbind(output_matrix, result_matrix)
  }
  output_matrix <- data.table(output_matrix)

  if (save == T) {
    write.table(output_matrix, quote = F, sep = "\t", row.names = F,
                file = paste(outputPrefix, "_halflife_calc_R2_selection.txt", sep=""))
  }

  return(output_matrix)
}

# Testing
inputFile <- normalized_table[121,]
group = c("Control","Knockdown")
hour = c(0, 1, 2, 4, 8, 12)
inforColumn = 4
CutoffTimePointNumber = 4
R2_criteria = 0.90
TimePointRemoval1 = c(1,2)
TimePointRemoval2 = c(8,12)
ThresholdHalfLife1 = 3
ThresholdHalfLife2 = 12
save = T
outputPrefix = "BridgeR_5"

halflife_table <- BridgeRHalfLifeCalcR2Select(normalized_table)

