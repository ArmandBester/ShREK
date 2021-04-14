library(shiny)
library(tidyverse)
library(readxl)
library(reshape2)
library(minpack.lm)
library(broom)
library(rJava)
library(xlsxjars)
library(xlsx)
library(drc)
library(intubate)
library(magrittr)

# Define UI 
ui <- fluidPage(
  
  titlePanel("Enzyme kinetics"),
  
  sidebarLayout(
    
    sidebarPanel(
      
      #browse for a file
      fileInput(inputId = "file1", label = "Browse", 
                buttonLabel = "Browse", placeholder = "Choose a file"),
      
      # select whether your file has column names or not    
      checkboxInput("header", "Header", TRUE),
      
      #select a model
      selectInput("SelMod", "Select Model",
                  c("Michaelis–Menten" = "MM",
                    "MM Substrate inhibition" = "MMs",
                    "Sigmoidal kinetics" = "sig",
                    "EC50" = "ec50")),
      
      #Axis labels
      textInput(inputId = "xlab", label = "X-axis label", value = "[S]"),
      
      textInput(inputId = "ylab", label = "Y-axis label", value = "Vo"),
      
      numericInput(inputId = "fontSize", label = "Change label font sice",
                   value = 15, min = 1, max = 35, step = 1),
      
      #Download data
      h3(textOutput("DownloadData", container = span)),
      downloadButton('downloadData', 'Download Data'),
      
      
      #Download plots
      h3(textOutput("DownloadPlots", container = span)),
      
      selectInput("SelPlot", "Select Plot",
                  c("MM Plot" = "MMPlot", "Residuals" = "resid")),
      
      numericInput("width", "Width:", value = 12),
      numericInput("hight", "Hight:", value = 8),
      
      selectInput("res", "Resolution:", 
                  c("150" = 150,"300" = 300, "600" = 600), selected = "300"),
      
      downloadButton('downloadPlot', 'Download Plot')
    ),
    
    
    mainPanel(
      tabsetPanel(tabPanel("Data",
                           
        h3(textOutput("Data", container = span)),
        dataTableOutput("contents"),
        br(),
        
        h3(textOutput("ModelTitle", container = span)),
        dataTableOutput("models"),
        br()
      ),
      
      
      tabPanel("Graphics",
      h3(textOutput("MMplot", container = span)),
      plotOutput("scatter"),
      br(),
      
      h3(textOutput("Residuals", container = span)),
      plotOutput("residuals"),
      br(),
      
      h3(textOutput("LBplot", container = span)),
      plotOutput("LB"),
      br(),

      h3(textOutput("EHplot", container = span)),
      plotOutput("EH")
      ))
    )
  )
) 



##################################################################################################
# Define server logic 
server <- shinyServer(function(input, output, session) {
  
  data <- NULL
  predict_range <- NULL
  Model <- NULL
  modelList <- NULL
  modelParameters <- NULL
  all_coefs <- NULL
  sheets <- NULL
  wb <- NULL
  
  plot_color <- "Set1"

  
  in_data <- reactive({
    modelList <<- ""
    modelList <<- list()
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)     
    else filePath <- inFile$datapath
    
    # names(data)[1] <<- "Conc"
    sheets <<- excel_sheets(filePath)
    data <<- NULL
    predict_range <<- NULL

    for (sheet in sheets){
      new_data <- NULL
      new_data <- read_excel(filePath, sheet = sheet, col_names = input$header)
      names(new_data)[1] <- "Conc"
      
      #################################################################################
      longDataset <- melt(new_data, id.var = ("Conc"))[-2]  #this is for model building
      longDataset <- longDataset[complete.cases(longDataset), ] #remove NA caused by some 
      # observations having less replicates
      #################################################################################
      
      #add row avg
      Vo_Avg <- apply(new_data[,2:ncol(new_data)], 1, mean, na.rm=TRUE)
      SD <- apply(new_data[,2:ncol(new_data)], 1, sd, na.rm=TRUE)
      new_data$Vo_Avg <- Vo_Avg
      new_data$SD <- SD
      
      #### set model ####
      if (input$SelMod == "MM"){
        Model <<- nlsLM(value ~ (Vm * Conc) / (Km + Conc), longDataset, 
                        start = list(Vm = max(longDataset$value), 
                                     Km = mean(longDataset$Conc)))
      }
      if (input$SelMod == "MMs"){
        Model <<- nlsLM(value ~ (Vm * Conc) / (Km + Conc * (1 + Conc/Ki)), longDataset, 
                        start = list(Vm = max(longDataset$value), 
                                     Km = mean(longDataset$Conc), 
                                     Ki = max(longDataset$Conc)))
      }
      if (input$SelMod == "sig"){
        Model <<- nlsLM(value ~ (Vm * Conc^H) / (Kprime + Conc^H), longDataset, 
                        start = list(Vm = max(longDataset$value), 
                                     Kprime = mean(longDataset$Conc),
                                     H = 1))
      }
      if (input$SelMod == "ec50"){
       Model <- drm(value ~ Conc, data = longDataset, fct = LL.4(fixed=c(NA, NA, NA, NA)))
      }
      
      
      
      
      
      #### set model ####
      
      #### Add predicted values to original data ####
      new_data$predicted <- predict(Model, new_data[,1])
      new_data$set <- sheet
      data <<- rbind(data, new_data)
      #### END  ####
      
      #### Code for getting a range of predicted values for smoothing ####
      new_predict_range <- data.frame(Conc = seq(0, max(data$Conc), length.out = 100))
      new_predict_range$predicted <- predict(Model, newdata = new_predict_range)
      new_predict_range$set <- sheet
      predict_range <<- rbind(predict_range, new_predict_range)
      #### Code for getting a range of predicted values for smoothing ####
      
      #### Model list append
      if (input$SelMod == "ec50"){
      modelList[[sheet]] <<- Model$coefficients
      }
      if (input$SelMod != "ec50"){
        modelList[[sheet]] <<- Model
      }
    }
    
    data$resid <<- data$Vo_Avg - data$predicted
    
    
    return(NULL)
    
  }) 
    

    
  output$Data <- renderText({
    "Your data and some calculated values"
  })
  
  output$DownloadData <- renderText({
    "Download Data"
  })
  
  output$DownloadPlots <- renderText({
    "Download Plots"
  })
  output$contents <- renderDataTable({
    call.me=in_data()
    data
  }, options = list(lengthMenu = list(c(5, 15, -1), c('5', '15', 'All')),
                    pageLength = 5))
  
  output$ModelTitle <- renderText({
    "Model summary"
  })
  output$models <- renderDataTable({
    call.me = in_data()
    all_coefs <<- plyr::ldply(modelList, tidy, .id = "model")
    all_coefs
  })
  
  ####  MM Scatter plot ####
  output$scatter <- renderPlot({
    call.me = in_data()
    ggplot(data = data, aes(x = data$Conc, y = data$Vo_Avg, color = data$set))+
      scale_colour_brewer(palette = plot_color)+
      geom_point(alpha = 0.8) +
      geom_line(data = predict_range, aes(x = predict_range$Conc, 
                                          y = predict_range$predicted, 
                                          color = predict_range$set)) +
      geom_errorbar(aes(x = data$Conc, ymin = data$Vo_Avg - data$SD, ymax = data$Vo_Avg + data$SD))+
      theme_classic()+
      theme(legend.title=element_blank())+
      theme(text = element_text(size = as.integer(input$fontSize)))+
      xlab(input$xlab) +
      ylab(input$ylab)
  })
  
  
  ####  resid plot ####
  output$residuals <- renderPlot({
    call.me = in_data()
    ggplot(data = data, aes(x = data$Conc, y = data$resid, color = data$set))+
      scale_colour_brewer(palette = plot_color)+
      geom_point(alpha = 0.8) +
      geom_hline(yintercept = 0, linetype=3)+
      theme_classic()+
      theme(legend.title=element_blank())+
      theme(text = element_text(size = as.integer(input$fontSize)))+
      xlab(input$xlab) +
      ylab("Residuals")
  })
  
  
  ####  PREPARE DL PLOTS ####
  MMplotInput <- reactive({
    call.me = in_data()
    ggplot(data = data, aes(x = data$Conc, y = data$Vo_Avg, color = data$set))+
      scale_colour_brewer(palette = plot_color)+
      geom_point(alpha = 0.8) +
      geom_line(data = predict_range, aes(x = predict_range$Conc,
                                          y = predict_range$predicted,
                                          color = predict_range$set)) +
      geom_errorbar(aes(x = data$Conc, ymin = data$Vo_Avg - data$SD, ymax = data$Vo_Avg + data$SD))+
      theme_classic()+
      theme(legend.title=element_blank())+
      theme(text = element_text(size = as.integer(input$fontSize)))+
      xlab(input$xlab) +
      ylab(input$ylab)
})
  
  ResidplotInput <- reactive({
    call.me = in_data()
    ggplot(data = data, aes(x = data$Conc, y = data$resid, color = data$set))+
      scale_colour_brewer(palette = plot_color)+
      geom_point(alpha = 0.8) +
      geom_hline(yintercept = 0, linetype=3)+
      theme_classic()+
      theme(legend.title=element_blank())+
      theme(text = element_text(size = as.integer(input$fontSize)))+
      xlab(input$xlab) +
      ylab("Residuals")
  })
  
  

  ########### DOWNLOAD PLOT SECTION ######################

  output$downloadPlot <- downloadHandler(
    filename = function() { paste("plot", '.tiff', sep='') },
    content = function(file) {
      device <- function(..., width, height) grDevices::tiff(..., 
                                                             width = input$width, 
                                                             height = input$hight, 
                                                             res = as.numeric(input$res), 
                                                             units = "in")
      ggsave(file, plot = switch(input$SelPlot,
                                 MMPlot = MMplotInput(),
                                 resid = ResidplotInput()), 
                                 device = device)
    }
  )

  
  ########### DOWNLOAD DATA SECTION ######################
  output$downloadData <- downloadHandler(
    filename = function() {filename = "data.xlsx"},
    
    content = function(file){
      call.me = in_data()
      wb<-createWorkbook(type="xlsx")
      TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=16, 
                                         isBold=TRUE, underline=1)
      SUB_TITLE_STYLE <- CellStyle(wb) + 
        Font(wb,  heightInPoints=14,
             isItalic=TRUE, isBold=FALSE)
      TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)
      TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
        Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
        Border(color="black", position=c("TOP", "BOTTOM"), 
               pen=c("BORDER_THIN", "BORDER_THICK")) 
      
      for (sheet in sheets){
        sheetData <- data %>% 
          filter(data$set == sheet) 
        # %>% 
        #   select(-c(set))
        
        sheetPredict <- predict_range %>% 
          filter(set == sheet)
        # %>% 
        #   select(-c(set))
        
        modelParams <- modelList[[sheet]]
        modelParams <- tidy(modelParams)
        
        XLsheet <- createSheet(wb, sheetName = sheet)
        
        #helper function
        xlsx.addTitle<-function(XLsheet, rowIndex, title, titleStyle){
          rows <-createRow(XLsheet,rowIndex=rowIndex)
          sheetTitle <-createCell(rows, colIndex=1)
          setCellValue(sheetTitle[[1,1]], title)
          setCellStyle(sheetTitle[[1,1]], titleStyle)
        }
        
        xlsx.addTitle(XLsheet, rowIndex=1, title = sheet,
                      titleStyle = TITLE_STYLE)
        xlsx.addTitle(XLsheet, rowIndex=2, 
                      title="Your and some calculated values",
                      titleStyle = SUB_TITLE_STYLE)
        
        
        #add user and calculated data
        addDataFrame(sheetData, XLsheet, startRow=3, startColumn=1, 
                     colnamesStyle = TABLE_COLNAMES_STYLE,
                     rownamesStyle = TABLE_ROWNAMES_STYLE)
        #setColumnWidth(XLsheet, colIndex=c(1:ncol(sheetData)), colWidth=11)
        
        #add pridicted data for smoothing
        addDataFrame(sheetPredict, XLsheet, startRow=3, startColumn=ncol(sheetData) + 4, 
                     colnamesStyle = TABLE_COLNAMES_STYLE,
                     rownamesStyle = TABLE_ROWNAMES_STYLE)
        #setColumnWidth(XLsheet, colIndex=c(1:ncol(sheetData)), colWidth=11)
        
        #add parameters from modelList
        addDataFrame(modelParams, XLsheet, startRow=nrow(sheetData) + 6, startColumn=1, 
                     colnamesStyle = TABLE_COLNAMES_STYLE,
                     rownamesStyle = TABLE_ROWNAMES_STYLE)
      
      }
      saveWorkbook(wb, "temp.xlsx")
      file.rename('temp.xlsx', file) 
    } 
  )
  
  
})# en of server


# Run the application 
shinyApp(ui = ui, server = server)





