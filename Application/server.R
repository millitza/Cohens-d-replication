library(shiny)
library(datasets)
library(xtable)
library(mvtnorm)
source("../function.R")

shinyServer(function(input, output) {
  
  # display original data set 
  origdata <- reactive({
    infile_orig <- input$dataset_orig        
    d_orig <- read.table(infile_orig$datapath, header = T)
    d_orig        
  })
  output$origdata <- DT::renderDataTable(DT::datatable(origdata(),options=list(paging=F,searching=F)))
  
  
  # display replication data set 
  repldata <- reactive({
    infile_repl <- input$dataset_repl       
    d_repl <- read.table(infile_repl$datapath, header = T)
    d_repl      
  })
  output$repldata <- DT::renderDataTable(DT::datatable(repldata(),options=list(paging=F,searching=F)))
  
  # set phi
  phi <- eventReactive(input$calc, {input$phi_user})
  # set iterations
  iter <- eventReactive(input$calc, {input$iter_user})
  
  # calculate BFs
  BF_output <- reactive({
    BFobject <- BF_repl_es(origdata(),repldata(),phi(),iter()+10000)
    BFobject
  })
  
  output$text1 <- renderUI({
  text <- paste("The presented Bayes factors were calculated using",BF_output()$repl_H1$iterations,"iterations including a burn-in period of ",BF_output()$repl_H1$burnin," iterations. For the second hypothesis, this number of iterations is used for each normal distribution in the prior on the effect size. The set value for phi is",BF_output()$repl_H3$phi,".")
    HTML(text)
  })
  
  output$BFtable <- renderUI({
    tab <- matrix(c(BF_output()$BF12,BF_output()$BF13),nrow=2,ncol=1)
    rownames(tab) <- c("BF12","BF13")
    colnames(tab) <- "Estimate"
    tab <- print(xtable(tab, align=rep("c", ncol(tab)+1)), 
                 floating=FALSE, tabular.environment="array", comment=FALSE, print.results=FALSE)
    html <- paste0( tab)
    list(
      withMathJax(HTML(html))
    )
  })
  
  output$plot <- renderPlot({
    phi <- BF_output()$repl_H3$phi
         xmin <- ((BF_output()$repl_H1$prior[2,1])-3*phi)
         xmax <- ((BF_output()$repl_H1$prior[2,1])+3*phi)
         x <- seq(xmin,xmax,.01)
         dens_H1 <- dnorm(x, mean=(BF_output()$repl_H1$prior[2,1]),sd=(BF_output()$repl_H1$prior[2,2]))
         plot(0:1,xlim=c(xmin,xmax),ylim=c(0,range(dens_H1)[2]),type="n",xlab="standardized mean difference in replication study",ylab="density",main="Prior densities on effect size under all three hypotheses")
         legend(x=xmax-phi,y=range(dens_H1)[2]-0.2,legend=c("H1","H2","H3"),lty=c(1,3,2),box.lty=0)
         lines(x,dens_H1,lty=1)
         
         dens_H2 <- 0.5*(dnorm(x, mean=(BF_output()$repl_H1$prior[2,1]-phi),sd=(BF_output()$repl_H1$prior[2,2]))+dnorm(x, mean=(BF_output()$repl_H1$prior[2,1]+phi),sd=(BF_output()$repl_H1$prior[2,2])))
         lines(x,dens_H2,lty=3)
         
         dens_H3 <- dnorm(x,mean=(BF_output()$repl_H1$prior[2,1]),sd=sqrt((BF_output()$repl_H1$prior[2,2])^2+phi^2))
         lines(x,dens_H3,lty=2)
  })
  

})