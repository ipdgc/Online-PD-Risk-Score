#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(ggplot2)
library(data.table)
library(optparse)
library(methods)
library(tools)
library(RColorBrewer)
library(shinydashboard)


#romeo <- readLines("romeo.txt")
#othello <- readLines("othello.txt")
#midsummers <- readLines("midsummers.txt")


###Dashboard Beginning for Text Analysis App ######
ui <- dashboardPage(skin = "blue",
                    dashboardHeader(title = "Polygenic Risk Score Calculator"),
                    ###Dashboard Sidebar Menu####
                    dashboardSidebar(
                      sidebarMenu(
                        ##Tab One
                        menuItem("File Upload",tabName = "file",icon = icon("file-text-o")),
                        ##Tab Two
                        menuItem("PRS Result",tabName = "text",icon = icon("bar-chart-o"))
                      )),
                    
                    ###Beginning of Dashboard Body <selection>####
                    dashboardBody(
                      tabItems(
                    ###File Upload Tab###
                        tabItem(tabName = "file",
                                fileInput("selection", "Upload Plink bfiles:",multiple = TRUE),
                                helpText(paste("Please upload Plink bfiles with the samples", 
                                               "you would like to analyze.")
                                         )
                                
                                ),
                        ###Text Output Tab####
                        tabItem(tabName = "text",
                                helpText(paste("This tab displays the PRS plot.")),
                                box(actionButton("prsTable","Create PRS Table"),
                                    br(),
                                    br(),
                                    downloadButton("downloadtwo", label="Download PRS Table"),
                                    downloadButton(outputId = "downloadsix",label = "Download PRS plot")
                                ),
                                br(),
                                br(),
                                plotOutput("plot2"))
                      ))
                    ###End of Dashboard Body####
)

###Dashboard End for Text Analysis App ######


# Define server logic required to run the Text Analysis App
server <- function(input, output, session) {
  
  ##Code for uploading Text File from User ##########
  
  ford <- reactive({ 
    
    req(input$selection) ## ?req #  require that the input is available
    
    inFile <- input$selection 
    
    df <- readLines(inFile$datapath)
    
    return(df)
    
  })
  
  
  
  ## Download code for wordcloud picture download ####
  
  output$download1 <- downloadHandler(
    filename = function() { paste("WordCloud",input$download3,sep = ".") },
    content = function(file) {
      if(input$download3=="png")
        png(file)
      else if (input$download3=="jpeg")
        jpeg(file)
      else if (input$download3=="bmp")
        bmp(file)
      else if (input$download3=="pdf")
        pdf(file)
      set.seed(1234)
      v <- terms()
      wordcloud(names(v),v, scale=c(6,0.5),
                min.freq = input$freq, max.words=input$max,
                rot.per=0.35,
                colors=brewer.pal(8, input$pal))
      dev.off()
    })
  
  
  
  ##Textbreakdown Download ###########
  
  output$downloadtwo <- downloadHandler(
    filename = function() { paste("TextBreakDown",input$name, sep='',".csv") },
    content = function(file) {
      write.csv(texterdf2(), file)
      
    })
  
  ##Emotion ggplot2 reactive download code for barplot###### 
  emotplot1 <- reactive({
    value<- ford()
    
    val_word <- get_tokens(value, pattern = "\\W")
    
    value <- get_nrc_sentiment(value)
    
    
    value <- as.data.frame(sort(colSums((prop.table(value[,1:8])))))
    
    colnames(value) <- "percentages"
    
    ggplot1<- ggplot(value, aes(x=sort(rownames(value),decreasing = FALSE), y=value$percentages)) +
      # plot the bars
      geom_bar(stat="identity", position="dodge",fill=input$colornow) +
      # create the label, "dodged" to fit the bars
      geom_text(aes(label=percent(value$percentages)), vjust=1, colour="white",
                position=position_dodge(.9), size=4)+labs(title="Emotional Sentiment",y = "Percentage",x="Emotion")+
      theme(panel.background = element_blank())
  })
  
  output$downloadseven <- downloadHandler(
    filename = function() { paste("Emotional Sentiment",'png',sep = ".") },
    content = function(file) {
      withProgress(message = 'Downloading BarPlot',
                   value = 0, {
                     for (i in 1:3) {
                       incProgress(1/3)
                       Sys.sleep(0.25)
                     }
                   },env = parent.frame(n=1))
      ggsave(file,emotplot1())})
  
  ##Positive vs Negative ggplot2 download code #######
  emotplot2 <- reactive({
    value<- ford()
    
    val_word <- get_tokens(value, pattern = "\\W")
    
    value <- get_nrc_sentiment(value)
    
    
    value <- as.data.frame(sort(colSums((prop.table(value[,9:10])))))
    
    colnames(value) <- "percentages"
    
    ggplot1<- ggplot(value, aes(x=sort(rownames(value),decreasing = FALSE), y=value$percentages))+
      # plot the bars
      geom_bar(stat="identity", position="dodge",fill=input$colornow2) +
      # create the label, "dodged" to fit the bars
      geom_text(aes(label=percent(value$percentages)), vjust=1, colour="white",
                position=position_dodge(.9), size=4)+labs(title="Positive vs. Negative Sentiment",y = "Percentage",x="Sentiment")+
      theme(panel.background = element_blank())
  })
  
  prs <- reactive({
    
    #for(i in 1:length(input$selection[,1])){
    #  lst[[i]] <- read.csv(input$selection[[i, 'datapath']])
    #}
    bfile <- input$selection[[, 'datapath']]
    plink_name <- input$selection[[, 'name']]
    
    file.rename(bfile[grep('bed',plink_name)],"test.bed")
    file.rename(bfile[grep('bim',plink_name)],"test.bim")
    file.rename(bfile[grep('fam',plink_name)],"test.fam")
    
    
    prs_per_sample <- system(paste0("Rscript PRSice.R --cov-file PCs.eigenvec_concate.txt --out \
$outfile -t test -b SPAIN3.meta --beta --snp MarkerName \
--A1 Allele1 --A2 Allele2 --stat Effect --se StdErr --pvalue P-value   \
--print-snp --score std --prsice . \
-n 24 --binary-target T --quantile 4 --prevalence 0.005 --bar-levels 5E-8,5E-7,5E-6,5E-5,5E-4,5E-3,5E-2,5E-1,1\
--fastscore
"))
    
    
    
  })
  
  observeEvent(input$wbdown,{output$wordbreakdown<-DT::renderDataTable({
    withProgress(message = 'Creating Word Breakdown',
                 value = 0, {
                   for (i in 1:3) {
                     incProgress(1/3)
                     Sys.sleep(0.25)
                   }
                 },env = parent.frame(n=1))
    
    #worddatabreakdown<- as.matrix.data.frame(texterdf3())  
    
    #wordatabreakdown <- worddatabreakdown[,1:2]
    #wordatabreakdown
    
    prs_per_sample
    
  })})
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
