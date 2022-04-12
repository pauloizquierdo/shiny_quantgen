ui <- pageWithSidebar( 
  
  headerPanel = headerPanel("Single Allele"),
  
  sidebarPanel(
    sliderInput("L",
                "Population size",
                value = 10000, min = 100, max = 10000),
    sliderInput("p",
                "allele frequency",
                value = 0.5, min = 0, max = 1),
    
    sliderInput(inputId = "H2",
                label = "heritability",
                value = 1, min = 0, max = 1),
    
    actionButton("goButton", "GO"),
  ), 
 
  mainPanel(type = "tabs",
            tabsetPanel(
  
              tabPanel("Only additive effects",
                       helpText(""),
                       
                       plotOutput(outputId = 'par_off')  
              ),
              tabPanel("single allele dominant",
                       helpText(""),
                       
                       plotOutput(outputId = 'out2')  
              ),
              tabPanel("single allele non-effect",
                       helpText(""),
                       
                       plotOutput(outputId = 'out3')  
              )
            )
  )
)


server <- function(input, output){
  
  rand <- eventReactive(input$goButton, {

    return(input)
  })
  
  output$par_off <- renderPlot({
    my.params<- rand()
    a1 <- 1000 # additivity
    d1 <- 0 #dominance
    X1 <- rbinom(n=my.params$L,size=2,p=my.params$p) # sampling genotypes under HWE
    # Genetic effect (0 = aa, 1 = Aa, 2 = AA)
    G <- ifelse(X1==2,a1,ifelse(X1==1,d1,-a1)) 
    vG <- var(G) # an estimate of the genetic variance
    vE <- vG*(1-my.params$H2)/my.params$H2 # with this error variance we will have the desired 
    E <- rnorm(n=n,sd=sqrt(vE)) # environmental effects for each sample
    Y <- G+E # Phenotype!!!

    hist(Y, breaks = 100, main = "Phenotype",
         col='darkolivegreen1')
    boxplot(Y~X1, col='darkolivegreen1',
            main = "Allele effect: Additive")
    layout(1:2)
    
    hist(Y, breaks = 50, main = "Phenotype",
         col='darkolivegreen1')
    boxplot(Y~X1, col='darkolivegreen1',
            main = "Allele effect: Additive")
    
    
  })    
  
  output$out2 <- renderPlot({
    my.params<- rand()
   
    #locus 1
    a2 <- 1000 # additivity
    d2 <- 800 #dominance
    X2 <- rbinom(n=my.params$L,size=2,p=my.params$p) # sampling genotypes under HWE
    
    # Genetic effect (0 = aa, 1 = Aa, 2 = AA)
    Gd <- ifelse(X2==2,a2,ifelse(X2==1,d2,-a2)) 
    
    vGd <- var(Gd) # an estimate of the genetic variance
    vEd <- vGd*(1-my.params$H2)/my.params$H2 # with this error variance we will have the desired 
    
    Ed <- rnorm(n=my.params$L,sd=sqrt(vEd)) # environmental effects for each sample
    
    Yd <- Gd+Ed # Phenotype!!!
    
    layout(1:2)
    hist(Yd, breaks = 100, main = "Population",
         col="salmon")
    boxplot(Yd~X2, main = "Allele effect: dominant", col="gold")
    
  })
  
  output$out3 <- renderPlot({
    my.params<- rand()
   
    #locus 1
    a1 <- 1000 # additivity
    d1 <- 0 #dominance
    
    #locus 2
    a2 <- 800
    d2 <- 0
    
    X1 <- rbinom(n=my.params$L,size=2,p=my.params$p) 
    X2 <- rbinom(n=my.params$L,size=2,p=my.params$p) 
    X3 <- rbinom(n=my.params$L,size=2,p=my.params$p) 
    
    # Genetic effect (0 = aa, 1 = Aa, 2 = AA)
    G <- ifelse(X1==2,a1,ifelse(X1==1,d1,-a1)) +
      ifelse(X2==2,a2,ifelse(X2==1,d2,-a2)) 
    
    vG <- var(G) # an estimate of the genetic variance
    vE <- vG*(1-my.params$H2)/my.params$H2 # with this error variance we will have the desired 
    
    En <- rnorm(n=my.params$L,sd=sqrt(vE)) # environmental effects for each sample
    
    Gn <- ifelse(X1==2,a1,ifelse(X1==1,0,-a1)) +
      ifelse(X2==2,a2,ifelse(X2==1,0,-a2))
    
    
    Yn <- Gn+En # Phenotype!!!
    
    layout(1:2)
    boxplot(Yn~X3, col= 'darkgoldenrod1',
            main = "Allele with non-effect")
    hist(Yn, breaks = 50, 
         col="coral", main = "Phenotype")
    })
  

  
}



# Run the application 
shinyApp(ui = ui, server = server)