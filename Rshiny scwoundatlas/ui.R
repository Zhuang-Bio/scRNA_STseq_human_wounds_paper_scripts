library(shiny) 
library(shinyhelper) 
library(data.table) 
library(Matrix) 
library(DT) 
library(magrittr) 
sc1conf = readRDS("sc1conf.rds")
sc1def  = readRDS("sc1def.rds")



### Start server code 
shinyUI(fluidPage( 
### HTML formatting of error messages 
 
tags$head(tags$style(HTML(".shiny-output-error-validation {color: red; font-weight: bold;}"))), 
list(tags$style(HTML(".navbar-default .navbar-nav { font-weight: bold; font-size: 16px; }"))), 
 
   
### Page title 
titlePanel("Single-cell Atlas of Human Skin Wound Healing"),  
navbarPage( 
  NULL,  
 ### Tab1.a1: cellInfo vs geneExpr on dimRed 
  tabPanel( 
    HTML("♦CellType vs GeneExpression"), 
    h4("Exploring the cell information vs gene expression on reduced dimensions"), 
    "In this tab, users can visualise both cell information and gene ",  
    "expression side-by-side on low-dimensional representions.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, h4("Dimension Reduction"), 
        fluidRow( 
          column( 
            12, selectInput("sc1a1drX", "X-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                           selected = sc1def$dimred[1]), 
            selectInput("sc1a1drY", "Y-axis:", choices = sc1conf[dimred == TRUE]$UI, 
                        selected = sc1def$dimred[2])) 
        ) 
      ), # End of column (6 space) 
      column( 
        3, h4("Subsetting cells"), #actionButton("sc1a1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = 1,#"input.sc1a1togL % 2 == 1", 
          selectInput("sc1a1sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp2), # set the default selection choice
          uiOutput("sc1a1sub1.ui"), 
          actionButton("sc1a1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1a1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ) 
      ), # End of column (6 space) 
      column( 
        6, h4("General graphics controls"), #actionButton("sc1a1tog0", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = 1, #"input.sc1a1tog0 % 2 == 1", 
          fluidRow( 
            column( 
              6, sliderInput("sc1a1siz", "Point size:", 
                             min = 0, max = 4, value = 1.25, step = 0.25), 
              radioButtons("sc1a1psz", "Plot size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE), 
              radioButtons("sc1a1fsz", "Font size:", 
                           choices = c("Small", "Medium", "Large"), 
                           selected = "Medium", inline = TRUE) 
            ), 
            column( 
              6, radioButtons("sc1a1asp", "Aspect ratio:", 
                              choices = c("Square", "Fixed", "Free"), 
                              selected = "Square", inline = TRUE), 
              checkboxInput("sc1a1txt", "Show axis text", value = TRUE) 
            ) 
          ) 
        ) 
      )  # End of column (6 space) 
    ),   # End of fluidRow (4 space) 
    fluidRow( 
      column( 
        6, style="border-right: 2px solid black", h4("Cell information"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a1inp1", "Group by:", 
                           choices = sc1conf$UI, 
                           selected = sc1def$meta1) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Cell information to colour cells by", 
                     content = c("Select cell information to colour cells", 
                                 "- Categorical covariates have a fixed colour palette", 
                                 paste0("- Continuous covariates are coloured in a ",  
                                        "Blue-Yellow-Red colour scheme, which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a1tog1", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a1tog1 % 2 == 1", 
              radioButtons("sc1a1col1", "Colour (Continuous data):", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "Blue-Yellow-Red"), 
              radioButtons("sc1a1ord1", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Original", inline = TRUE), 
              checkboxInput("sc1a1lab1", "Show cell info labels", value = TRUE) 
            ) 
          ) 
        ), 
        fluidRow(column(12, uiOutput("sc1a1oup1.ui"))), 
        downloadButton("sc1a1oup1.pdf", "Download PDF"), 
        downloadButton("sc1a1oup1.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a1oup1.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a1oup1.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5))
      ), # End of column (6 space) 
      column( 
        6, h4("Gene expression"), 
        fluidRow( 
          column( 
            6, selectInput("sc1a1inp2", "Gene name:", choices=NULL) %>%  
              helper(type = "inline", size = "m", fade = TRUE, 
                     title = "Gene expression to colour cells by", 
                     content = c("Select gene to colour cells by gene expression", 
                                 paste0("- Gene expression are coloured in a ", 
                                        "White-Red colour scheme which can be ", 
                                        "changed in the plot controls"))) 
          ), 
          column( 
            6, actionButton("sc1a1tog2", "Toggle plot controls"), 
            conditionalPanel( 
              condition = "input.sc1a1tog2 % 2 == 1", 
              radioButtons("sc1a1col2", "Colour:", 
                           choices = c("White-Red","Blue-Yellow-Red","Yellow-Green-Purple"), 
                           selected = "White-Red"), 
              radioButtons("sc1a1ord2", "Plot order:", 
                           choices = c("Max-1st", "Min-1st", "Original", "Random"), 
                           selected = "Max-1st", inline = TRUE) 
            ) 
          ) 
        ) , 
        fluidRow(column(12, uiOutput("sc1a1oup2.ui"))), 
        downloadButton("sc1a1oup2.pdf", "Download PDF"), 
        downloadButton("sc1a1oup2.png", "Download PNG"), br(), 
        div(style="display:inline-block", 
            numericInput("sc1a1oup2.h", "PDF / PNG height:", width = "138px", 
                         min = 4, max = 20, value = 6, step = 0.5)), 
        div(style="display:inline-block", 
            numericInput("sc1a1oup2.w", "PDF / PNG width:", width = "138px", 
                         min = 4, max = 20, value = 8, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  ),     # End of tab (2 space) 
 

 ### Tab1.c1: violinplot / boxplot 
 tabPanel( 
    HTML("♦Single GeneExp ViolinPlot"),  
   h4("Cell information / gene expression violin plot / box plot"), 
   "In this tab, users can visualise the gene expression or continuous cell information ",  
   "across groups of cells (e.g. condition / cell types).", 
   br(),br(), 
   fluidRow( 
     column( 
       3, style="border-right: 2px solid black", 
       selectInput("sc1c1inp1", "Cell information (X-axis):", 
                   choices = sc1conf[grp == TRUE]$UI, 
                   selected = sc1def$grp1) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Cell information to group cells by",  
                content = c("Select categorical cell information to group cells by",  
                            "- Single cells are grouped by this categorical covariate",  
                            "- Plotted as the X-axis of the violin plot / box plot")),  
       selectInput("sc1c1inp2", "Gene name (Y-axis):", choices=NULL) %>%  
         helper(type = "inline", size = "m", fade = TRUE, 
                title = "Gene to plot", 
                content = c("Select gene to plot on Y-axis", 
                            "- Input the gene symbol")), 
       h4("Subsetting cells"), #actionButton("sc1c1togL", "Toggle to subset cells"), 
       conditionalPanel( 
         condition = 1,#"input.sc1c1togL % 2 == 1", 
         selectInput("sc1c1sub1", "Cell information to subset:", 
                     choices = sc1conf[grp == TRUE]$UI, 
                     selected = sc1def$grp2), 
         uiOutput("sc1c1sub1.ui"), 
         actionButton("sc1c1sub1all", "Select all groups", class = "btn btn-primary"), 
         actionButton("sc1c1sub1non", "Deselect all groups", class = "btn btn-primary") 
       ), br(),
       radioButtons("sc1c1typ", "Plot type:", 
                    choices = c("violin", "boxplot"), 
                    selected = "violin", inline = TRUE), 
       checkboxInput("sc1c1pts", "Show data points", value = FALSE),  br(), 
       actionButton("sc1c1tog", "Toggle graphics controls"), 
       conditionalPanel( 
         condition = "input.sc1c1tog % 2 == 1", 
         sliderInput("sc1c1siz", "Data point size:",  
                     min = 0, max = 4, value = 1.25, step = 0.25),  
         radioButtons("sc1c1psz", "Plot size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE), 
         radioButtons("sc1c1fsz", "Font size:", 
                      choices = c("Small", "Medium", "Large"), 
                      selected = "Medium", inline = TRUE)) 
     ), # End of column (6 space) 
     column(9, uiOutput("sc1c1oup.ui"),  
            downloadButton("sc1c1oup.pdf", "Download PDF"),  
            downloadButton("sc1c1oup.png", "Download PNG"), br(), 
            div(style="display:inline-block", 
                numericInput("sc1c1oup.h", "PDF / PNG height:", width = "138px", 
                             min = 4, max = 20, value = 8, step = 0.5)), 
            div(style="display:inline-block", 
                numericInput("sc1c1oup.w", "PDF / PNG width:", width = "138px", 
                             min = 4, max = 20, value = 10, step = 0.5)) 
     )  # End of column (6 space) 
   )    # End of fluidRow (4 space) 
 ),     # End of tab (2 space) 
 

  ### Tab1.d1: Multiple gene expr 
  tabPanel( 
    HTML("♦Multiple GeneExp BubblePlot"), 
    h4("Multiple gene expressions bubbleplot / heatmap"), 
    "In this tab, users can visualise the gene expression patterns of ", 
    "multiple genes grouped by categorical cell information (e.g. condition / cell types).", br(), 
    "The normalised expression are averaged, log-transformed and then plotted.", 
    br(),br(), 
    fluidRow( 
      column( 
        3, style="border-right: 2px solid black", 
        textAreaInput("sc1d1inp", HTML("List of gene names <br /> 
                                          (Max 50 genes, separated <br /> 
                                           by , or ; or newline):"), 
                      height = "200px", 
                      value = paste0(sc1def$genes, collapse = ", ")) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "List of genes to plot on bubbleplot / heatmap", 
                 content = c("Input genes to plot", 
                             "- Maximum 50 genes (due to ploting space limitations)", 
                             "- Genes should be separated by comma, semicolon or newline")), 
        selectInput("sc1d1grp", "Group by:", 
                    choices = sc1conf[grp == TRUE]$UI, 
                    selected = sc1conf[grp == TRUE]$UI[3]) %>% 
          helper(type = "inline", size = "m", fade = TRUE, 
                 title = "Cell information to group cells by", 
                 content = c("Select categorical cell information to group cells by", 
                             "- Single cells are grouped by this categorical covariate", 
                             "- Plotted as the X-axis of the bubbleplot / heatmap")), 
        h4("Subsetting cells"),#actionButton("sc1d1togL", "Toggle to subset cells"), 
        conditionalPanel( 
          condition = 1,#"input.sc1d1togL % 2 == 1", 
          selectInput("sc1d1sub1", "Cell information to subset:", 
                      choices = sc1conf[grp == TRUE]$UI, 
                      selected = sc1def$grp1), 
          uiOutput("sc1d1sub1.ui"), 
          actionButton("sc1d1sub1all", "Select all groups", class = "btn btn-primary"), 
          actionButton("sc1d1sub1non", "Deselect all groups", class = "btn btn-primary") 
        ), br(), 
        radioButtons("sc1d1plt", "Plot type:", 
                     choices = c("Bubbleplot", "Heatmap"), 
                     selected = "Bubbleplot", inline = TRUE), 
        checkboxInput("sc1d1scl", "Scale gene expression", value = TRUE), 
        checkboxInput("sc1d1row", "Cluster rows (genes)", value = TRUE), 
        checkboxInput("sc1d1col", "Cluster columns (samples)", value = FALSE), 
        br(), 
        actionButton("sc1d1tog", "Toggle graphics controls"), 
        conditionalPanel( 
          condition = "input.sc1d1tog % 2 == 1", 
          radioButtons("sc1d1cols", "Colour scheme:", 
                       choices = c("White-Red", "Blue-Yellow-Red", 
                                   "Yellow-Green-Purple"), 
                       selected = "Blue-Yellow-Red"), 
          radioButtons("sc1d1psz", "Plot size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE), 
          radioButtons("sc1d1fsz", "Font size:", 
                       choices = c("Small", "Medium", "Large"), 
                       selected = "Medium", inline = TRUE)) 
      ), # End of column (6 space) 
      column(9, h4(htmlOutput("sc1d1oupTxt")), 
             uiOutput("sc1d1oup.ui"), 
             downloadButton("sc1d1oup.pdf", "Download PDF"), 
             downloadButton("sc1d1oup.png", "Download PNG"), br(), 
             div(style="display:inline-block", 
                 numericInput("sc1d1oup.h", "PDF / PNG height:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)), 
             div(style="display:inline-block", 
                 numericInput("sc1d1oup.w", "PDF / PNG width:", width = "138px", 
                              min = 4, max = 20, value = 10, step = 0.5)) 
      )  # End of column (6 space) 
    )    # End of fluidRow (4 space) 
  )      # End of tab (2 space) 
   ,    
br(), 
p("", style = "font-size: 125%;"), 
p(em("This tool was built on"), a("ShinyCell", 
  href = "https://github.com/SGDDNB/ShinyCell",target="_blank")), 
tags$footer(
        HTML(
          '
           <p align="center" style="margin:10px;font-size:15px">© 2023-2024<br>
                <a href="https://www.xulandenlab.com/tools" target="_blank">Xu Land&#233;n Laboratory | Karolinska Institutet</a> 
           </p>
          '
          )
    )
))) 
 
 
 
 