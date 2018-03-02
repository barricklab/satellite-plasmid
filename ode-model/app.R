#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(deSolve)
library(ggplot2)
library(shiny)
library(tidyr)
library(dplyr)

## Explanation of variables


## Different subpopulations
##
## X1 = wild-type plasmid
## X2 = satellite plasmid
## X3 = deletion plasmid
##
## X4 = wild-type plasmid, integration
## X5 = satellite plasmid, integration
## X6 = deletion plasmid, integration
##
## X7 = no plasmid, integration
## X8 = no plasmid (dead if antibiotic present)


## These are absolute fitnesses, so they must be >1
## for the OD to work properly

w1 = 0.45
w2 = 0.81
w3 = 0.90

w4 = 0.45
w5 = 0.81
w6 = 0.90

w7 = 0.90
w8 = 1.0

## u_del = rate of mutation to deletion plasmid
## u_sat = rate of generating satellite plasmid
## u_int = rate of integrating plasmid
##
## u_wt_loss = rate of losing plasmid from cell with wild-type plasmid
## u_sat_loss = rate of losing plasmid from cell with satellite plasmid

u_del = 1E-5
u_sat = 10E-5
u_int = 1E-7
u_wt_del_loss = 0.5**40
u_sat_loss = 0.5**12

#u_int = 0
#u_wt_del_loss = 0
#u_sat_loss = 0

# Note = fitnesses must be converted to natural log ln(w) in the continuous time model
# This should be done outside of (before) the ODE function

satellite_plasmid_ode <- function(t, y, p) {
  with(as.list(c(y, p)), {
    dX1 <- w1*X1 - u_del*w1*X1 - u_sat*w1*X1 - u_int*w1*X1
    dX2 <- w2*X2 + u_sat*w1*X1 - u_int*w4*X2
    dX3 <- w3*X3 + u_del*w1*X1 - u_int*w5*X3
    dX4 <- w4*X4 + u_int*w1*X1 - u_wt_del_loss*w4*X4
    dX5 <- w5*X5 + u_int*w4*X2 - u_sat_loss*w5*X5
    dX6 <- w6*X6 + u_int*w5*X3 - u_wt_del_loss*w6*X6
    dX7 <- w7*X7 + u_wt_del_loss*w4*X4 + u_wt_del_loss*w6*X6 + u_sat_loss*w5*X5
    dX8 <- w8*X8 + u_wt_del_loss*w1*X1 + u_wt_del_loss*w3*X3 + u_sat_loss*w2*X2
    return(list(c(dX1, dX2, dX3, dX4, dX5, dX6, dX7, dX8)))
  })
}



plot_ode <- function(p) {
  
  
  fixed_no_plasmid_fitness = 0
  if (p$w8 != 0) {
    fixed_no_plasmid_fitness = log(p$w8)
  }
  
  #Convert to log fitnesses for continuous time model
  parms <- c(
    w1=log(p$w1), 
    w2=log(p$w2), 
    w3=log(p$w3), 
    w4=log(p$w4), 
    w5=log(p$w5), 
    w6=log(p$w6),
    w7=log(p$w7),
    w8=fixed_no_plasmid_fitness,
    u_del=p$u_del,
    u_sat=p$u_sat,
    u_int=p$u_int,
    u_wt_del_loss=p$u_wt_del_loss,
    u_sat_loss=p$u_sat_loss,
    end_generation=p$end_generation
    )
  
  #Run model
  out <- ode(y = c(X1=1, X2=0, X3=0, X4=0, X5=0, X6=0, X7=0, X8=0), times=seq(0, p$end_generation, 0.1), satellite_plasmid_ode, parms)
  
  #Normalize the 
  out = data.frame(out)
  out$total = out$X1 + out$X2 + out$X3 + out$X4 + out$X5 + out$X6 + out$X7 + out$X8
  out$X1 = out$X1 / out$total
  out$X2 = out$X2 / out$total
  out$X3 = out$X3 / out$total
  out$X4 = out$X4 / out$total
  out$X5 = out$X5 / out$total
  out$X6 = out$X6 / out$total
  out$X7 = out$X7 / out$total
  out$X8 = out$X8 / out$total
  
  out = gather(out, key = "type", value = "n" ,X1 ,X2, X3, X4, X5, X6, X7, X8) %>% 
    select(-total) %>% 
    rename(generation=time)
  out$type = factor(out$type)
  levels(out$type) = c("wt-plasmid", 
                       "satellite-plasmid", 
                       "deletion-plasmid", 
                       "wt-plasmid-integration", 
                       "satellite-plasmid-integration", 
                       "deletion-plasmid-integration", 
                       "no-plasmid-integration",
                       "no-plasmid"
  )
  ggplot(out, aes(x=generation,y=n, color=type)) + geom_line() + scale_y_log10(limits=c(1E-10, 1))
}


# Define UI for application
ui <- fluidPage(
  numericInput("w1", label = "wt-plasmid-fitness", value = w1),
  numericInput("w2", label = "satellite-plasmid-fitness", value = w2),
  numericInput("w3", label = "deletion-plasmid-fitness", value = w3),
  numericInput("w4", label = "wt-plasmid-integration-fitness", value = w4),
  numericInput("w5", label = "satellite-plasmid-integration-fitness", value = w5),
  numericInput("w6", label = "deletion-plasmid-integration-fitness", value = w6),
  numericInput("w7", label = "no-plasmid-integration-fitness", value = w7),
  numericInput("w8", label = "no-plasmid-fitness", value = w8),
  numericInput("u_del", label = "deletion-mutation-rate", value = u_del),
  numericInput("u_sat", label = "satellite-mutation rate", value = u_sat),
  numericInput("u_int", label = "integration-mutation-rate", value = u_int),
  numericInput("u_wt_del_loss", label = "wt-or-deletion-plasmid-loss-rate", value = u_wt_del_loss),
  numericInput("u_sat_loss", label = "satellite-plasmid-loss-rate", value = u_sat_loss),
  
  numericInput("end_generation", label = "end generation", value = 100),
  plotOutput("satellite_plasmid")
)

# Define server logic 
server <- function(input, output) {
  
#  output$satellite_plasmid <- renderPlot({
#    parms <- c(w1=input$w1, w2=input$w2)
#    out <- ode(y = c(X1=1, X2=1), times=seq(0, 10, .1), satellite_plasmid, parms)
#    matplot.0D(out)
#  })
  
  
  output$satellite_plasmid <- renderPlot({
    
    
    if (input$w8 != 0) {
      min_w = 0.5*min(input$w1, input$w2, input$w3, input$w4, input$w5, input$w6, input$w7, input$w8)
      fixed_w8 = input$w8 / min_w
    } else {
      min_w = 0.5*min(input$w1, input$w2, input$w3, input$w4, input$w5, input$w6, input$w7)
      fixed_w8 = 0
    }
    parms <- data.frame( 
      w1=input$w1/min_w, 
      w2=input$w2/min_w, 
      w3=input$w3/min_w,
      w4=input$w4/min_w,
      w5=input$w5/min_w,
      w6=input$w6/min_w,
      w7=input$w7/min_w,
      w8=fixed_w8,
      u_del=input$u_del,
      u_sat=input$u_sat,
      u_int=input$u_int,
      u_wt_del_loss=input$u_wt_del_loss,
      u_sat_loss=input$u_sat_loss,
      end_generation=input$end_generation)
    plot_ode(parms) 
  })
}

# Code for testing outside of Shiny
parms <- data.frame( 
  w1=w1, 
  w2=w2, 
  w3=w3,
  w4=w4,
  w5=w5,
  w6=w6,
  w7=w7,
  w8=0,
  u_del=u_del,
  u_sat=u_sat,
  u_int=u_int,
  u_wt_del_loss=u_wt_del_loss,
  u_sat_loss=u_sat_loss,
  end_generation=100
  )
plot_ode(parms) 

# Run the application 
shinyApp(ui = ui, server = server)

