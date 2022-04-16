#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(stringr)
library(BSgenome)
library(BiocManager)
library(parallel)
library(shinyFiles)
library(Rsamtools)
library(dplyr)
library(tidyverse)

options(shiny.maxRequestSize = 10000 * 1024 ^ 2) # Up to 10GB file size

# Define UI for application that draws a histogram
shinyUI(navbarPage("TSS ExtractR", collapsible = T, theme = shinytheme("cerulean"),
                 tabPanel("Quality Control",
                          fluidPage(
                            #titlePanel("TSS ExtractR"),
                            sidebarLayout(
                              sidebarPanel(
                                h3(strong("Preprocess and Quality Control")),
                                numericInput(
                                  "threads",
                                  "Number of Threads",
                                  value = detectCores() / 2,
                                  min = 1,
                                  max = detectCores()
                                ),
                                hr(style = "border-color: lightgrey"),
                                h4(strong("Trim Reads"), align = "center"),
                                radioButtons(
                                  "adaptors", "Select Adaptor", 
                                  choices = c("From the list", "Custom Adaptor"),
                                  selected = "From the list"
                                ),
                                uiOutput("adaptorComp"),
                                numericInput(
                                  "minlength",
                                  "Minimum Read's Length",
                                  value = 25,
                                  min = 0
                                ),
                                hr(style = "border-color: lightgrey"),
                                shinyFilesButton("fastqfile", "Select FastQ file", title = "Select FastQ file", multiple = F),
                                checkboxInput("pairend", "Pair-end Sequencing", value = FALSE),
                                selectInput("aligner", "Select Aligner routine", c("Rbowtie", "Rhisat2")),
                                numericInput("complx", "Complexity (in bits)", value = 0.4, min = 0, max = 1),
                                actionButton(inputId = "submit", label = div("Submit", icon("arrow-right")))),
                              
                              # Show a plot of the generated distribution
                              mainPanel(textOutput("al_stats"))
                            )
                          )), 
                 tabPanel("Process BAM", fluidPage(
                   sidebarLayout(
                     sidebarPanel(
                       h3(strong("BAM File Analysis")),
                       checkboxInput("testsmp", "Test Sample"),
                       uiOutput("testfile"),
                       checkboxInput("controlsmp", "Control Sample"),
                       uiOutput("controlfile"),
                       div(style = "margin-top: 50px"),
                       actionButton(inputId = "submit1", label = div("Submit", icon("arrow-right")))),
                     
                     # Show a plot of the generated distribution
                     mainPanel(dataTableOutput("testout"), dataTableOutput("controlout"))
                   )
                 )),
                 tabPanel("Cluster TSSs", fluidPage(
                   # Sidebar with a slider input for number of bins
                   sidebarLayout(
                     sidebarPanel(
                       h3(strong("TSS clustering and Data Analysis")), 
                       textOutput("dfckecker"),
                       div(style = "margin-top: 30px"),
                       fileInput("tssfinalFile", "Or you can upload Existing and Unclustered CSV File", multiple = FALSE,
                                 accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                       hr(style = "border-color: lightgrey"),
                       checkboxInput("header", "Header", TRUE),
                       radioButtons("sep", "Separator", choices = c(Comma = ",",
                                                Semicolon = ";",
                                                Tab = "\t"), selected = ","),
                       radioButtons("quote", "Quote", choices = c(None = "",
                                                "Double Quote" = '"',
                                                "Single Quote" = "'"), selected = '"'),
                       hr(style = "border-color: lightgrey")),
                     
                     mainPanel(h1(textOutput("clTitle")), dataTableOutput("clustered"))
                   )
                 )),
                 tabPanel("Satistical Analysis", fluidPage(
                   sidebarLayout(
                     sidebarPanel(
                       h3(strong("Statistical Analysis using Beta-Binomial Distribution")),
                       textOutput("cldfckecker"),
                       div(style = "margin-top: 30px"),
                       fileInput("clusteredFile", "Or you can upload Existing and Clustered CSV File", multiple = FALSE,
                                 accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                       hr(style = "border-color: lightgrey"),
                       checkboxInput("header1", "Header", TRUE),
                       radioButtons("sep1", "Separator", choices = c(Comma = ",",
                                                                    Semicolon = ";",
                                                                    Tab = "\t"), selected = ","),
                       radioButtons("quote1", "Quote", choices = c(None = "",
                                                                  "Double Quote" = '"',
                                                                  "Single Quote" = "'"), selected = '"'),
                       hr(style = "border-color: lightgrey"),
                       numericInput("pvalue", "P-value cutoff", value = 0.001, min = 0, max = 0.1, step = 0.001)),
                     
                     mainPanel(h1(textOutput("clTitle1")), dataTableOutput("clustered1"))
                   )
                 )),
                 tabPanel("Extract Knowledge", fluidPage(
                   sidebarLayout(
                     sidebarPanel(
                       h3(strong("Data Plots Visualization"))),
                     
                     mainPanel()
                   )
                 )),
                 tabPanel("TSS Annotation", fluidPage(
                   sidebarLayout(
                     sidebarPanel(
                       h3(strong("TSS Classification"))),
                     
                     mainPanel()
                   )
                 )),
                 tabPanel("Info", fluidPage(
                   h3(strong("Info"))
                 ))
                ))
