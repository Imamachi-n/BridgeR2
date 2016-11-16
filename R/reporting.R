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
#' @param SearchRowName Row name for searching.
#'
#' @param inforColumn The number of information columns.
#'
#' @param color color of line graph for two decay curve.
#'
#' @param TimePointRemoval1 The candicate_1 of time point removal.
#'
#' @param TimePointRemoval2 The candicate_2 of time point removal.
#'

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
  stopifnot(is.character(SearchRowName))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.character(color) && is.vector(color))
  stopifnot(is.numeric(ThresholdHalfLife1))
  stopifnot(is.numeric(ThresholdHalfLife2))

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
    infor_st <- 1 + (group_index - 1) * (time_points + InforColumn + 3)
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

  # plot wrapper
  ggplotly_decay_curve <- function(data){
    ggplot(data,
           aes(x = hour, y = exp, colour = factor(Condition)))
    p <- p + geom_point() + scale_color_manual(values = c("black", "orange"))
    p <- p + stat_smooth(method = lm, se=FALSE, fullrange=TRUE)
    return(p)
  }

  # UI - shiny dashboard
  ui_body <- tabItem(tabName = "charts",
                     fluidRow(
                       tabBox(width = 12,
                              # first tab
                              tabPanel('Data Table',
                                       fluidRow(
                                         box(title = tagList(icon("table"), "BRIC-seq results"), status = "primary", solidHeader = TRUE, width = 12,
                                             DT::dataTableOutput('table_BRIC')
                                         )
                                       )
                              ),
                              # second tab
                              tabPanel('Plot',
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
                                             plotOutput("plot1",
                                                        width = 400, height = 400,
                                                        click = "plot1_click",
                                                        brush = brushOpts(id = "plot1_brush"))
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
                       )
                     )
  )

  ui <- dashboardBody(
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
      # exp data prep
      data <- as.vector(as.matrix(inputFile[input$text]))
      gene_name <- as.character(input$text)

      # index search
      fig_data <- NULL
      data_label <- comparisonFile
      for (table_index in 1:length(comp_file_index)) {
        # index infor
        exp_index <- as.numeric(exp_index_table[table_index, ])
        model_index <- as.numeric(model_index_vec[table_index])

        # extract exp/model data
        exp <- data[exp_index[1]:exp_index[2]]
        model <- data[model_index]

        # data table
        time_point_exp_raw <- data.frame(hour, exp)
        check <- parse_rm_hr_infor(model)    # util.R
        check_index <- sapply(check,
                              function(t) which(time_point_exp_raw$hour == t))
        time_point_exp_del <- time_point_exp_raw[-check_index,]
        time_point_exp_del <- time_point_exp_del[time_point_exp_del$exp > 0,]
        label <- rep(data_label[1], nrow(time_point_exp_del))

        # result
        fig_data <- rbind(fig_data,
                          cbind(time_point_exp_del, Condition = label))
      }

      # plotting
      ggplotly_decay_curve(fig_data)
    })

    # Information box - condition_1
    output$controlBox <- renderInfoBox({
      infoBox(
        comparisonFile[1], model_index_vec[1], icon = icon("line-chart"),
        color = "navy", fill = TRUE
      )
    })

    # Information box - condition_2
    output$treatedBox <- renderInfoBox({
      infoBox(
        comparisonFile[2], model_index_vec[2], icon = icon("line-chart"),
        color = "yellow", fill = TRUE
      )
    })

    # data table
    output$mytable1 <- renderTable({
      table_data <- NULL
      for (table_index in 1:length(comp_file_index)) {
        # index infor
        exp_index <- as.numeric(exp_index_table[table_index, ])
        R2_index <- as.numeric(model_index_vec[table_index]) + 1
        halflife_index <- R2_index + 1

        # extract exp/model data
        exp <- data[exp_index[1]:exp_index[2]]
        R2 <- data[R2_index]
        halflife <- data[halflife_index]
        table_data <- rbind(table_data, c(exp, R2, halflife))
      }
      hour_label <- sapply(hour, function(t) paste(t, "hr", sep=""))
      colnames(table_data) <- c(hour_label, "R2", "Half-life")
      rownames(table_data) <- comparisonFile

      return(table_data)
    })

  })
}



library(shiny)
library(shinydashboard)
library(plotly)
library(ggplot2)
library(data.table)
library(DT)

# ggplot2 testing
test_exp_table <- data.frame(hour=c(0,1,2,4,8,12, 0,1,2,4,8),
                             exp=c(1,0.7,0.5,0.3,0.2,0.1, 1,0.9,0.8,0.6,0.3),
                             label=c("black","black","black","black","black","black",
                                     "red","red","red","red","red","red"))

p <- ggplot(test_exp_table,
            aes(x = hour, y = exp, colour = factor(label)))
p <- p + geom_point() + scale_color_manual(values = c("black", "orange"))
p <- p + stat_smooth(method = lm, se=FALSE, fullrange=TRUE)
ggplotly(p)

# plotly testing
test <- data.frame(mtcars)

library(shiny)
library(plotly)

ui <- fluidPage(
  plotlyOutput("plot")
  # verbatimTextOutput("event")
)

server <- function(input, output) {

  # renderPlotly() also understands ggplot2 objects!
  output$plot <- renderPlotly({
    plot_ly(test, x = ~mpg, y = ~wt)
  })

  # output$event <- renderPrint({
  #   d <- event_data("plotly_hover")
  #   if (is.null(d)) "Hover on a point!" else d
  # })
}

shinyApp(ui, server)
