shinyUI(navbarPage(
  title = 'Replication of Cohen\'s d',
  tabPanel('Start',verticalLayout(
    fileInput("dataset_orig", "Upload your original data set"),
    fileInput("dataset_repl", "Upload your replication data set"),
    numericInput("phi_user", paste("Select a value of",expression(phi)),value=0.5,min=0,step=0.1),
    numericInput("iter_user",paste("Select the number of posterior draws / iterations after burn-in"),10000,min=0,step=1000),
    actionButton("calc", "Calculate evidence for effect size replication")
  )),
  tabPanel('Original data set',DT::dataTableOutput('origdata')),
  tabPanel('Replication data set',DT::dataTableOutput('repldata')),
   tabPanel('Results',
            verticalLayout(
              uiOutput("text1"),
              uiOutput("BFtable"),
              plotOutput("plot")
             ))

))
