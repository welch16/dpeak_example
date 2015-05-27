shinyUI(
  fluidPage(
    sidebarLayout(
      sidebarPanel(
        numericInput("nreads","Number of reads",500,min=50,max=5000),        
        numericInput("S","Genomic space (S)",1000),
        sliderInput("fl","Fragment length",min=25,max=300,value=150,step=1),
        sliderInput("p_D","Strand probability",value=.5,min=0,max=1,step=.01),
        textInput("text_w","BS weights (one for each weight except background, don't need to sum one)","10;50"),
        textInput("text_m","BS positions (need to be between 1 and S)", "500"),
        numericInput("delta","Delta",75,min=1,step=1),
        numericInput("sigma","Sigma",50,min=0),
        numericInput("beta","Beta",40,min=0)
                ),
      mainPanel(
        navbarPage(
          title = "Exploration",
          tabPanel("Peak",plotOutput("peak")),
          tabPanel("Scale-space",
                   actionButton("simulate","Update"),            
                   plotOutput("peakMat")
                   ),
          ## tabPanel("Strand imbalance",
          ##          actionButton("calc_fsr","Update"),
          ##          plotOutput("strandRatio")
          ##          ),
          tabPanel("Reads",verbatimTextOutput("printReads"))
                  )
                )
    )
  )
)
