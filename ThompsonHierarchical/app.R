# to upload, run rsconnect::deployApp('ThompsonHierarchical')

library(shiny)
rm(list = ls())
source("ThompsonHierarchical.R")
source("ReadDataApp.R")


ui <- fluidPage(
#  titlePanel(h1("Optimal treatment assignment given covariates", align = "center")),
  sidebarLayout(
    sidebarPanel(
      width=6,
      fileInput("file1", "Choose CSV file of previous covariates, treatments, and outcomes",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      fileInput("file2", "Choose CSV file of covariates for new wave",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      hr(),
      actionButton(inputId = "calcbutton", label = "Calculate treatment assignment"),
      hr(),
      downloadButton("downloadData", "Download treatment assignment")
      
    ),
    
    mainPanel(
      width=6,
      includeMarkdown("instructions.md"),
      hr(),
      tableOutput("treatmentcounts"),
      hr(),
      tableOutput("designtable")
    )
  )  
)



server <- function(input, output, session) {
  
  v = reactiveValues()
  
  observeEvent(input$calcbutton,{
    req(input$file1)
    req(input$file2)
    
    #loading files
    #make sure that options are chosen right!
    priordata=ReadDataApp(input$file1$datapath)    
    newwave=ReadDataApp(input$file2$datapath,
                        priordata$key)    
    
    #calculating treatment assignment
    newwave$Dstar=as.integer(
                DtchoiceThompsonHierarchical(priordata$Y,priordata$D,priordata$X, #outcomes, treatments, and covariates thus far
                                              priordata$k,priordata$nx, #number of treatments and number of strata
                                              newwave$Xt))
    v$newwave =rename(newwave,
                       stratum =Xt,
                       treatment=Dstar)
    
    v$treatmentcounts=    v$newwave %>%
      group_by(stratum, treatment) %>%
      summarise(count=n())
  })
  
  output$designtable =  renderTable({
    v$newwave
   })
  
  output$treatmentcounts =  renderTable({
    v$treatmentcounts
  })
 
#download optimal design
  output$downloadData <- downloadHandler(
    filename = "treatmentassignment.csv",
    content = function(file) {
      write_csv(v$newwave, file)
    }
  )
  

 
}

# Run the app ----
shinyApp(ui = ui, server = server)