generate_infor_st_ed <- function(group_index,
                                 time_points,
                                 inforColumn){
  infor_st <- 1 + (group_index - 1)*(time_points + inforColumn)
  infor_ed <- (inforColumn)*group_index + (group_index - 1)*time_points
  return(c(infor_st, infor_ed))
}


generate_hour_label <- function(group,
                                hour,
                                sample_index){
  hour_label <- NULL
  for (t in hour) {
    if (t == 0) {
      next
    }
    label <- t
    if (t < 10) {
      label <- paste("0", t, sep="")
    }
    hour_label <- append(hour_label,
                         paste(label, "hr_",group[sample_index], sep=""))
  }
  return(hour_label)
}


generate_fig_log10_matrix <- function(exp_data,
                                      label){
  exp_data <- t(exp_data)    # Inverse
  exp_data <- factor(exp_data)
  exp_data <- as.numeric(as.character(exp_data))    # exp_data for fig
  exp_data <- log10(exp_data)

  gene_number <- length(exp_data[[1]])
  label_data <- rep(label, gene_number)

  fig_data <- data.frame(exp=exp_data, label=factor(label_data))
  # fig_data <- fig_data[!is.infinite(fig_data[,1]),]
  return(fig_data)
}


parse_rm_hr_infor <- function(model){
  check <- gsub("Delete_","",model)
  check <- gsub("hr","",check)
  check <- as.numeric(strsplit(check,"_")[[1]])
  return(check)
}
