#' Shinyapp reporting to draw RNA decay curve.
#'
#' \code{BridgeReport} returns a shinyapp to draw RNA decay curve.
#'
#' @param inputFile The vector of tab-delimited matrix file.
#'
#' @param group The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param comparisonFile The vector of group names.
#'
#' @param searchRowName Row name for searching.
#'
#' @param inforColumn The number of information columns.
#'
#' @param color color of line graph for two decay curve.
#'
#' @param TimePointRemoval1 The candicate_1 of time point removal.
#'
#' @param TimePointRemoval2 The candicate_2 of time point removal.
#'
#' @export
#'
#' @import shiny shinydashboard data.table ggplot2
#' @importFrom plotly ggplotly renderPlotly plotlyOutput

BridgeReport <- function(inputFile,
                         group = c("Control","Knockdown"),
                         hour = c(0, 1, 2, 4, 8, 12),
                         comparisonFile = c("Control","Knockdown"),
                         searchRowName = "symbol",
                         inforColumn = 4,
                         color = c("black","red"),
                         TimePointRemoval1 = c(1,2),
                         TimePointRemoval2 = c(8,12)){

  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.character(comparisonFile) && is.vector(comparisonFile))
  stopifnot(is.character(searchRowName))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.character(color) && is.vector(color))
  stopifnot(is.numeric(TimePointRemoval1))
  stopifnot(is.numeric(TimePointRemoval2))

  # data infor
  time_points <- length(hour)
  group_number <- length(group)
  comp_file_index <- as.vector(sapply(comparisonFile,
                                      function(test) which(group == test)))

  # key setting
  setkeyv(inputFile, searchRowName)

  # removed time points label prep
  TimePointRemoval1 <- sort(TimePointRemoval1, decreasing = T)
  TimePointRemoval2 <- sort(TimePointRemoval2, decreasing = F)

  create_delete_time_label <- function(test_time){
    times_length <- length(test_time) #c(12, 8) => 2, c(12) => 1

    times_index <- c(times_length)
    add_index <- times_length
    label_list <- NULL
    for(counter in times_length:1){
      # excepted_time_points
      check_times <- test_time[times_index] #c(24), c(24,12), c(24,12,8)
      time_point_exp_del_label <- paste("Delete_",paste(check_times,collapse="hr_"),"hr",sep="")
      label_list <- append(label_list, time_point_exp_del_label) #

      # counter
      add_index <- add_index - 1
      times_index <- append(times_index, add_index)
    }
    return(label_list)
  }

  rm_time_points1 <- create_delete_time_label(TimePointRemoval1)
  rm_time_points2 <- create_delete_time_label(TimePointRemoval2)
  select_model_label <- c("Raw", rm_time_points1, rm_time_points2)

  # table header infor
  header_st <- inforColumn + 1    # 0hr_exp
  header_ed <- inforColumn + time_points + 3    # half-life
  table_header <- colnames(inputFile[header_st:header_ed])

  # index search prep
  exp_index_table <- NULL
  model_index_vec <- NULL

  for (group_index in comp_file_index) {
    # exp index prep
    infor_st <- 1 + (group_index - 1) * (time_points + inforColumn + 3)
    infor_ed <- inforColumn * group_index + (group_index - 1) * (time_points + 3)
    exp_st <- infor_ed + 1
    exp_ed <- infor_ed + time_points
    exp_st_ed <- data.frame(exp_st = c(exp_st),
                            exp_ed = c(exp_ed))
    exp_index_table <- rbind(exp_index_table, exp_st_ed)

    # model index prep
    model_index <- infor_ed + time_points + 1
    model_index_vec <- c(model_index_vec, model_index)
  }

  # exp index example
  #   st ed
  # 1  5 10
  # 2 18 23

  # model index example
  #   model
  # 1    11
  # 2    24

  # function_4
  decay_calc_for_shiny <- function(data,
                                   comparisonFile,
                                   comp_file_index,
                                   select1,
                                   select2){
    # index search
    fig_data <- NULL
    predicted_data <- NULL
    table_data <- NULL
    data_label <- comparisonFile
    model_list <- c(select1, select2)
    # model_list <- c("Raw", "Raw")

    for (table_index in 1:length(comp_file_index)) {
      # index infor
      exp_index <- as.numeric(exp_index_table[table_index, ])
      model_index <- as.numeric(model_index_vec[table_index])

      # extract exp/model data
      exp <- as.numeric(as.vector(data[exp_index[1]:exp_index[2]]))
      # model <- data[model_index]
      model <- model_list[table_index]

      # data table
      time_point_exp_raw <- data.frame(hour, exp)
      time_point_exp_del <- NULL
      if (model == "Raw") {
        time_point_exp_del <- time_point_exp_raw
      } else {
        check <- parse_rm_hr_infor(model)    # util.R
        check_index <- sapply(check,
                              function(t) which(time_point_exp_raw$hour == t))
        time_point_exp_del <- time_point_exp_raw[-check_index,]
      }
      time_point_exp_del <- time_point_exp_del[as.numeric(as.vector(time_point_exp_del$exp)) > 0,]
      label <- rep(data_label[table_index], nrow(time_point_exp_del))

      # model fitting
      fitting <- lm(log(as.numeric(as.vector(exp))) ~ hour - 1, time_point_exp_del)
      fitting_model <- summary(fitting)

      # halflife calc
      coef <- -fitting_model$coefficients[1]
      half_life <- log(2) / coef
      if(coef < 0 || half_life >= 24){
        half_life <- 24
      }

      # R2 calc
      R2 <- fitting_model$r.squared

      # table data prep
      table_data <- rbind(table_data, c(exp, R2, half_life))

      # predicted data prep
      predicted <- exp(predict(fitting, data.frame(hour=time_point_exp_del$hour)))

      # result
      predicted_data <- rbind(predicted_data, data.frame(hour=time_point_exp_del$hour,
                                                         exp=predicted,
                                                         Condition = label))
      fig_data <- rbind(fig_data,
                        cbind(time_point_exp_del, Condition = label))
    }

    # table data
    hour_label <- sapply(hour, function(t) paste(t, "hr", sep=""))
    colnames(table_data) <- c(hour_label, "R2", "Half-life")
    rownames(table_data) <- comparisonFile

    return(list(fig_data, predicted_data, table_data))
  }

  # ggplotly wrapper
  ggplotly_decay_curve <- function(data,
                                   predicted,
                                   gene_name,
                                   y_range){
    p <- ggplot(data,aes(x = as.numeric(hour), y = as.numeric(as.vector(exp)), colour = factor(Condition)))
    p <- p + geom_point(size=2.5, shape=19)
    p <- p + scale_color_manual(values = c("black", "orange"))
    p <- p + geom_line(data = predicted, size=0.9)
    p <- p + ggtitle(gene_name)
    p <- p + xlab("Time (hr)")
    p <- p + ylab("Relative RPKM (Time0 = 1)")
    p <- p + scale_x_continuous(breaks=seq(0, 12, by=2), limits = c(0, 12))
    ybreaks <- c(0.01, 0.1, 1, 10)
    p <- p + scale_y_log10(breaks=ybreaks,
                           labels=ybreaks,
                           limits=c(y_range[1], y_range[2]))
    p <- p + theme(title = element_text(size=15), axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
    p <- p + theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))
    # p <- p + theme(legend.position = "none")
    return(p)
  }

  # plotly wrapper
  plotly_decay_curve <- function(data,
                                 predicted,
                                 gene_name){
    p <- plot_ly(data, x = ~hour, y = ~exp(exp)) %>%
      layout(yaxis = list(type = "log"))
    return(p)
  }

  # UI - shiny dashboard
  ui_body <- dashboardBody(
    fluidRow(
      box(title = tagList(icon("check-square"), "Inputs"), status = "primary", solidHeader = TRUE, width = 4, # background = "navy",
          helpText("Select your favorite gene symbol to examine."),
          textInput("text",
                    label = "Input gene symbol",
                    placeholder = "Search..."),
          sliderInput("range_x",
                      label = "X-axis(Time course):",
                      min = 0, max = 12, value = c(0, 12)),
          sliderInput("range_y",
                      label = "Y-axis(Relative RNA remaining):",
                      min = 0.001, max = 10, value = c(0.001,1.5)),
          uiOutput("selectInput1"),
          uiOutput("selectInput2"),
          submitButton(tagList(icon("search"), "Update View"))
      ),

      box(title = tagList(icon("line-chart"), "RNA decay"), status = "primary", solidHeader = TRUE, width = 8,
          plotlyOutput("plot1",
                       width = 500, height = 400)
      ),

      infoBoxOutput("controlBox"),
      infoBoxOutput("treatedBox")
    ),
    fluidRow(
      box(title = tagList(icon("gear"), "Detail information"), status = "primary", solidHeader = TRUE, width = 12,
          tableOutput("mytable1")
      )
    )
  )

  ui <- dashboardPage(
    dashboardHeader(title = "BridgeReport"),
    dashboardSidebar(disable = TRUE),
    ui_body
  )

  # Server- define server logic required to draw RNA decay curve
  server <- shinyServer(function(input, output){
    # select time points
    # condition_1
    output$selectInput1 <- renderUI({
      return(selectInput("select1",
                         label = paste(comparisonFile[1],
                                       "(Select time points)", sep=""),
                         choices = select_model_label,
                         selected = 1))
    })

    # condition_2
    output$selectInput2 <- renderUI({
      return(selectInput("select2",
                         label = paste(comparisonFile[2],
                                       "(Select time points)", sep=""),
                         choices = select_model_label,
                         selected = 1))
    })

    # Draw RNA decay curve
    output$plot1 <- renderPlotly({
      # Start
      if (input$text == '') {
        p <- ggplot()
        return(p)
      }

      # exp data prep
      data <- as.vector(as.matrix(inputFile[input$text]))
      gene_name <- as.character(input$text)

      # data <- as.vector(as.matrix(inputFile["GAPDH"]))
      # gene_name <- as.character("GAPDH")

      data_list <- decay_calc_for_shiny(data,
                                        comparisonFile,
                                        comp_file_index,
                                        input$select1,
                                        input$select2)

      # data prep
      fig_data <- data_list[[1]]
      predicted_data <- data_list[[2]]
      table_data <- data_list[[3]]

      # table
      output$mytable1 = renderTable({
        table_data
      })

      # plotting
      return(ggplotly(ggplotly_decay_curve(fig_data,
                                           predicted_data,
                                           gene_name,
                                           input$range_y)))

      # return(plotly_decay_curve(fig_data,
      #                           predicted_data,
      #                           gene_name))
    })

    # Information box - condition_1
    output$controlBox <- renderInfoBox({
      data <- as.vector(as.matrix(inputFile[input$text]))
      infoBox(
        comparisonFile[1], data[model_index_vec[1]], icon = icon("line-chart"),
        color = "navy", fill = TRUE
      )
    })

    # Information box - condition_2
    output$treatedBox <- renderInfoBox({
      data <- as.vector(as.matrix(inputFile[input$text]))
      infoBox(
        comparisonFile[2], data[model_index_vec[2]], icon = icon("line-chart"),
        color = "yellow", fill = TRUE
      )
    })
  })

  # Run the application
  shinyApp(ui = ui, server = server)
}



# library(shiny)
# library(shinydashboard)
# library(plotly)
# library(ggplot2)
# library(data.table)
#
# pvalue_table <- fread("C:/Users/Naoto/OneDrive/Shiny_app/For_Git/BridgeR2/tmp/BridgeR_6_halflife_pvalue_evaluation.txt", header = T)
# BridgeReport(pvalue_table)
#
# inputFile <- pvalue_table
# group = c("Control","Knockdown")
# hour = c(0, 1, 2, 4, 8, 12)
# comparisonFile = c("Control","Knockdown")
# searchRowName = "symbol"
# inforColumn = 4
# color = c("black","red")
# TimePointRemoval1 = c(1,2)
# TimePointRemoval2 = c(8,12)

