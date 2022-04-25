library(shiny)
library(rrBLUP) # GS ans GWAS
library(dplyr) # organize data
library(qqman) # Manhattan plots
library(reshape)
library(grDevices)


ui <- fluidPage(
    h1("Quantitative genetics Shiny App"),
    p(style = "font-family:Impact",
    a("Workshop website", href="https://pauloizquierdo.github.io/Quantitative_Genetics/")
    ),
    
    titlePanel("1. Single allele"),
    fluidRow(
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
        
        actionButton("goButton1", "GO"),
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
                  )))
    ),
    
    titlePanel("2. Multiple allele"),
    fluidRow(
        sidebarPanel(
            sliderInput("L2",  "Population size",
                        value = 500, min = 50, max = 500),
            sliderInput("qtlnum2","number of qtl",
                        value = 20, min = 2, max = 50),
            sliderInput(inputId = "h2",
                        label = "heritability",
                        value = 0.5, min = 0, max = 1),
            actionButton("goButton2", "GO"),
        ), 
        
        mainPanel(type = "tabs",
                  tabsetPanel(
                      tabPanel("Allele effect",
                               helpText("Please allow 5 secs to load the figures!"),
                               plotOutput(outputId = 'par_off2')
                      ),
                      tabPanel("GWA",
                               helpText("Please allow 5 secs to load the figures!"),
                               plotOutput(outputId = 'par_off21')
                      )))
    ),
    
    titlePanel("3. Genomic prediction"),
    fluidRow(
        sidebarPanel(
            sliderInput("L3",
                        "Population size",
                        value = 1000, min = 50, max = 1000),
            sliderInput("qtlnum3",
                        "number of qtl",
                        value = 550, min = 2, max = 550),
            sliderInput(inputId = "h3",
                        label = "heritability",
                        value = 0.99, min = 0, max = 1),
            actionButton("goButton3", "GO"),
        ), 
        mainPanel(type = "tabs",
                  tabsetPanel(
                      tabPanel("Prediction accuracy",
                               helpText(""),
                               
                               plotOutput(outputId = 'par_off3')  
                      ),
                      tabPanel("Observed vs Predicted",
                               helpText(""),
                               
                               plotOutput(outputId = 'out4')  
                      )))
        )
    )

server <- function(input, output) {
    
rand <- eventReactive(input$goButton1, {
        
        return(input)
    })
    
    output$par_off <- renderPlot({
        my.params<- rand()
        a1 <- 1000 # additivity
        d1 <- 0 #dominance
        X1 <- rbinom(n=my.params$L,size=2,p=my.params$p) 
        G <- ifelse(X1==2,a1,ifelse(X1==1,d1,-a1)) 
        vG <- var(G) 
        vE <- vG*(1-my.params$H2)/my.params$H2 
        E <- rnorm(n=my.params$L,sd=sqrt(vE)) 
        Y <- G+E 
        
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
    
    
rand2 <- eventReactive(input$goButton2, {
        #parameters
        
        #package data for plotting
        return(input)
    })

    output$par_off2 <- renderPlot({
        my.params2<- rand2()
        popsize=my.params2$L2
        qtlnum=my.params2$qtlnum2
        h2=my.params2$h2
        M <- matrix(0,1100,popsize) # matrix of genotypes dimensions
        
        for (i in 1:popsize) { # fill the genotype matrix
            M[,i] <- ifelse(runif(1100)<0.5,-1,1)
        }
        
        colnames(M) <- 1:popsize
        
        # create the data frame with chromosomes
        geno <- data.frame(SNP=1:1100, 
                           Chromosome=rep(1:11,each=100),
                           Position=rep(1:100, times=11),
                           M,check.names=FALSE)
        
        QTL <- sample(1:1100, qtlnum)
        
        u <- rep(0,1100) # QTL positions
        u[QTL] <- rep(0.5,qtlnum) # QTL effect. This simulation is going to assume additive effects in all QTL
        
        g <- as.vector(crossprod(M,u)) #  cross product between Markers and marker effects
        
        
        y <- g + rnorm(popsize,mean=0,sd=sqrt((1-h2)/h2*var(g))) # phenotype with QTL effects and size pop = popsize
        
        pheno <- data.frame(Sample=1:popsize,Trait=y) # get the data frame input for the GWA 
        
        scores <- GWAS(pheno,geno,plot=F) # run the mixed model
        
        scores <- scores %>% mutate(pvalue= 10^-Trait) # transform the -log10 to pvalue
        
        Y_RA <- pheno[,2] # phenotype to estimate marker effect
        M_qtl1 <- M[QTL[1],] # Firts QTL
        
        fmR1=lm(Y_RA~M_qtl1) # fits a linear model on genotypes 
        
        R1A1A1=fmR1$coef[1] #(i.e., the intercept)
        R1A2A2=fmR1$coef[1]+fmR1$coef[2] # Homocygous A2A2, E[X=2]
        
        exp_values_R1 = matrix(c(R1A1A1, R1A2A2), 
                               ncol = 2, byrow = T)
        
        colnames(exp_values_R1) = c("aa", "AA")
        row.names(exp_values_R1) = c("QTL1")
        
        exp_values_R1 = as.table(exp_values_R1)
        
        
        boxplot(Y_RA~M_qtl1, xaxt = "n",
                main = "QTL", 
                xlab = "Genotype",
                ylab = "Trait",
                col=c("darkolivegreen1", "cornflowerblue")) 
        axis(1, at = c(1, 2),
             labels = c("aa", "AA")) 
        abline(h=c(R1A1A1,R1A2A2), col=c("red", "blue"), lty=c(1,2), lwd=c(2, 2))
        
        
    }) 
    ####
    output$par_off21 <- renderPlot({
        my.params21<- rand2()
        popsize=my.params21$L2
        qtlnum=my.params21$qtlnum2
        h2=my.params21$h2
        M <- matrix(0,1100,popsize) # matrix of genotypes dimensions
        
        for (i in 1:popsize) { # fill the genotype matrix
            M[,i] <- ifelse(runif(1100)<0.5,-1,1)
        }
        
        colnames(M) <- 1:popsize
        
        # create the data frame with chromosomes
        geno <- data.frame(SNP=1:1100, 
                           Chromosome=rep(1:11,each=100),
                           Position=rep(1:100, times=11),
                           M,check.names=FALSE)
        
        QTL <- sample(1:1100, qtlnum)
        
        u <- rep(0,1100) # QTL positions
        u[QTL] <- rep(0.5,qtlnum) # QTL effect. This simulation is going to assume additive effects in all QTL
        
        g <- as.vector(crossprod(M,u)) #  cross product between Markers and marker effects
        
        
        y <- g + rnorm(popsize,mean=0,sd=sqrt((1-h2)/h2*var(g))) # phenotype with QTL effects and size pop = popsize
        
        pheno <- data.frame(Sample=1:popsize,Trait=y) # get the data frame input for the GWA 
        
        scores <- GWAS(pheno,geno,plot=F) # run the mixed model
        
        scores <- scores %>% mutate(pvalue= 10^-Trait) # transform the -log10 to pvalue
        
        Y_RA <- pheno[,2] # phenotype to estimate marker effect
        M_qtl1 <- M[QTL[1],] # Firts QTL
        
        fmR1=lm(Y_RA~M_qtl1) # fits a linear model on genotypes 
        
        R1A1A1=fmR1$coef[1] #(i.e., the intercept)
        R1A2A2=fmR1$coef[1]+fmR1$coef[2] # Homocygous A2A2, E[X=2]
        
        exp_values_R1 = matrix(c(R1A1A1, R1A2A2), 
                               ncol = 2, byrow = T)
        
        colnames(exp_values_R1) = c("aa", "AA")
        row.names(exp_values_R1) = c("QTL1")
        
        exp_values_R1 = as.table(exp_values_R1)
        
        manhattan(scores, main = "Manhattan plot", # manhattan plot 
                  chr="Chromosome",
                  bp="Position", snp="SNP", p="pvalue", 
                  highlight = QTL, 
                  col = c("blue4", "orange3"),
                  suggestiveline = F, 
                  genomewideline = -log10(0.05/1100), # Bonferroni
                  cex = 0.6)
        
    }) 
    
    rand3 <- eventReactive(input$goButton3, {
        #parameters
        
        #package data for plotting
        return(input)
    })
    
    output$par_off3 <- renderPlot({
        my.params3<- rand3()
        popsize=my.params3$L3
        qtlnum=my.params3$qtlnum3
        h2=my.params3$h3
        M <- matrix(rep(0,popsize*550),550,popsize) # matrix of genotypes dimensions
        
        for (i in 1:popsize) { # fill the genotype matrix
            M[,i] <- ifelse(runif(550)<0.5,-1,1)
        }
        
        geno <- data.frame(SNP=1:550, # create the data frame with chromosomes
                           Chromosome=rep(1:11,each=50),
                           Position=rep(1:50, 
                                        times=11),
                           M,check.names=FALSE)
        
        QTL <- sample(1:550, qtlnum)
        
        u <- rep(0,550) 
        u[QTL] <- rep(0.5,qtlnum) # QTL effects
        g <- as.vector(crossprod(M,u)) # get the cross 
        
        y <- g + rnorm(popsize,sd=sqrt((1-h2)/h2*var(g))) # phenotype in a random distribution on 500 samples  with  the mean and sd desired
        
        pheno <- data.frame(line=1:popsize,y=y) # get the data frame input for the GWA 
        
        phenoGS <- as.matrix(pheno[,2]) # sacale data (subtract the mean and divide by sd). This could help to increase the prediction accuracy
        
        markers <- t(M) 
        cycles = 10 # number of cycles 
        accuracy = matrix(nrow = cycles, ncol=1) # empty matrix to store results
        
        for(r in 1:cycles) { # open loop to run GP model
            train = as.matrix(sample(1:popsize, 250)) # select 250 lines to train the model
            test = setdiff(1:popsize, train) #select remainder samples to validate the model
            
            yTRN  = phenoGS[train,] 
            xTRN = markers[train,]
            
            yTST  = phenoGS[test,]
            xTST = markers[test,]
            
            yNA = phenoGS # phenotype 
            yNA[test] <- NA 
            fmGS0 = mixed.solve(yTRN,  Z=xTRN)
            
            yHat0 = xTST %*% fmGS0$u 
            # store the correlation of each cycle
            accuracy[r,1] = cor(yHat0, yTST, use = "complete" ) 
            
        }
        predGS0 <- data.frame(yTST, yHat0) # create the data frame input for GP plot
        
        boxplot(accuracy, ylab= "Prediction accuracy",col="salmon")
        
    })    
    
    output$out4 <- renderPlot({
        my.params3<- rand3() 
        popsize=my.params3$L3
        qtlnum=my.params3$qtlnum3
        h2=my.params3$h3
        M <- matrix(rep(0,popsize*550),550,popsize) # matrix of genotypes dimensions
        
        for (i in 1:popsize) { # fill the genotype matrix
            M[,i] <- ifelse(runif(550)<0.5,-1,1)
        }
        
        geno <- data.frame(SNP=1:550, # create the data frame with chromosomes
                           Chromosome=rep(1:11,each=50),
                           Position=rep(1:50, 
                                        times=11),
                           M,check.names=FALSE)
        
        QTL <- sample(1:550, qtlnum)
        
        u <- rep(0,550) 
        u[QTL] <- rep(0.5,qtlnum) # QTL effects
        g <- as.vector(crossprod(M,u)) # get the cross 
        
        y <- g + rnorm(popsize,sd=sqrt((1-h2)/h2*var(g))) # phenotype in a random distribution on 500 samples  with  the mean and sd desired
        
        pheno <- data.frame(line=1:popsize,y=y) # get the data frame input for the GWA 
        
        phenoGS <- as.matrix(pheno[,2]) # sacale data (subtract the mean and divide by sd). This could help to increase the prediction accuracy
        
        markers <- t(M) # Transpose of M
        
        cycles = 10 # number of cycles 
        accuracy = matrix(nrow = cycles, ncol=1) # empty matrix to store results
        
        for(r in 1:cycles) { # open loop to run GP model
            train = as.matrix(sample(1:popsize, 250)) # select 250 lines to train the model
            test = setdiff(1:popsize, train) #select remainder samples to validate the model
            
            yTRN  = phenoGS[train,] #
            xTRN = markers[train,]
            
            yTST  = phenoGS[test,] 
            xTST = markers[test,]
            
            yNA = phenoGS # phenotype 
            yNA[test] <- NA 
            fmGS0 = mixed.solve(yTRN, Z=xTRN)
          
            yHat0 = xTST %*% fmGS0$u  
            accuracy[r,1] = cor(yHat0, yTST, use = "complete" ) 
        }
        predGS0 <- data.frame(yTST, yHat0) # create the data frame input for GP plot
        plot(predGS0$yHat0, predGS0$yTST, 
             xlab = "Observed", ylab = "Predicted")
        abline(lm(predGS0$yTST ~ predGS0$yHat0),   
               col="red")
    })
      
}


shinyApp(ui = ui, server = server)
