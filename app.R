# Load packages ----
set.seed(1)
library(shiny)
library(bslib)
library(tidyverse)

# Load data ----
gtex_hr <- read_tsv('data/GTEx_HR_combined.tsv') # GTEx HR expression data
df <- readxl::read_excel('data/hormone_receptor_list.xlsx') # Load hormone and HR info
df$Hormone_name1 <- recode(df$Hormone_name1, 'orphan' = 'Unknown')
df$Hormone_chemical_classes <- recode(df$Hormone_chemical_classes, 'NA' = 'Unknown')
df$Hormone1_tissue <- recode(df$Hormone1_tissue, 'NA' = 'Unknown')
df <- df %>% separate_rows(Hormone1_tissue, sep = ', ')
df <- df %>% filter(Hormone1_tissue != 'Unknown') # Remove orphan receptors
vec_receptor <- df$Gene_symbol
vec_hormone <- sort(df$Hormone_name1 %>% unique())
vec_tissue <- sort(gtex_hr$SMTS %>% unique())

# Source helper functions -----
source("sankey.R")

# Some settings
options(dplyr.summarise.inform = FALSE)

# User interface ----
ui <- fluidPage(
  theme = bs_theme(
    bootswatch = "flatly",
    base_font = font_google("Roboto"),
    primary = "#007BFF",
    secondary = "#6C757D"
  ),
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "style.css")
  ),
  
  fluidRow(
    column(
      width = 12,
      div(class = "header",
    img(src = "logo.png", alt = "Logo"),
    h1("HOHC: Human Organ Hormonal Communication")
    )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3, # 调整侧边栏宽度
      tabsetPanel(
        type = "tabs",
        tabPanel("Basic Filters",
                 radioButtons("sex", label = "Choose sex to display:",
                              choices = c("Male", "Female"), selected = "Male"),
                 selectInput("hormone_class", 
                             label = "Hormone class:",
                             choices = c("All", "Amino acid-derived", "Lipid-derived", "Peptide"),
                             selected = "All"),
                 selectInput('hormone_name', 
                             label = 'Hormone name:', 
                             choices = vec_hormone, 
                             selected = c("Androgen", "Glucagon-Like Peptide 1"), 
                             multiple = TRUE)
        ),
        tabPanel("Organ Filters",
                 selectInput('hormone_tissue', 
                             label = 'Outgoing organ (hormone organ):', 
                             choices = vec_tissue, 
                             selected = NULL, 
                             multiple = TRUE),
                 selectInput('receptor_tissue', 
                             label = 'Incoming organ (receptor organ):', 
                             choices = vec_tissue, 
                             selected = NULL, 
                             multiple = TRUE)
        ),
        tabPanel("Receptor Filters",
                 selectInput("receptor_class", 
                             label = "Receptor class:", 
                             choices = c("All", "Membrane HR", "Nuclear HR"), 
                             selected = "All"),
                 selectInput('receptor_name', 
                             label = 'Receptor name (HGNC gene symbol):', 
                             choices = vec_receptor, 
                             selected = NULL, 
                             multiple = TRUE)
        ),
        tabPanel("Advanced Options",
                 sliderInput("number_HRtissue", 
                             label = "Number of organs expressing a hormone receptor:", 
                             min = 1, max = 10, value = 3),
                 sliderInput("tpm_log1p", 
                             label = "Receptor expression range (log1p TPM):", 
                             min = 0, max = 6.35, value = c(0.1, 6.35))
        )
      ),
      hr(),
      div(
        style = "display: flex; justify-content: space-between; align-items: center; margin-top: 10px;",
        actionButton("reset_filters", "Reset Filters", class = "btn-secondary"),
        downloadButton("download_data", "Download Data", class = "btn-primary")
      )
    ),
    mainPanel(
      width = 9, # 增加主面板宽度
      div(
        class = "main-panel",
        textOutput("text"),
        plotOutput("plot", height = "500px")
      )
    )
  ),
  
  # Footer with Copyright
  fluidRow(
    column(
      width = 12,
      div(
        class = "footer",
        style = "background-color: #f8f9fa; border-top: 1px solid #cccccc; padding: 10px; text-align: center;",
        p("© 2024 HOHC. All rights reserved.", style = "color: #6c757d; font-size: 12px; margin: 0;")
      )
    )
  )
)

# Server logic ----
server <- function(input, output, session) {
  
  # Reset filters to default
  observeEvent(input$reset_filters, {
    updateRadioButtons(session, "sex", selected = "Male")
    updateSelectInput(session, "hormone_class", selected = "All")
    updateSelectInput(session, "hormone_name", selected = c("Androgen", "Glucagon-Like Peptide 1"))
    updateSelectInput(session, "hormone_tissue", choices = vec_tissue, selected = NULL)
    updateSelectInput(session, "receptor_tissue", choices = vec_tissue, selected = NULL) # 修复 receptor_tissue
    updateSelectInput(session, "receptor_class", selected = "All")
    updateSelectInput(session, "receptor_name", choices = vec_receptor, selected = NULL)
    updateSliderInput(session, "number_HRtissue", value = 3)
    updateSliderInput(session, "tpm_log1p", value = c(0.1, 6.35))
  })
  
  # Other existing server logic remains unchanged
  output$text <- renderText({
    res <- plot_sankey(input$sex, input$hormone_class, input$hormone_name, input$hormone_tissue,
                       input$receptor_class, input$receptor_name, input$receptor_tissue,
                       input$number_HRtissue, input$tpm_log1p[1], input$tpm_log1p[2],
                       gtex_hr, df, vec_receptor)
    
    if (is.character(res)) {
      return(res)
    } else {
      return(NULL)
    }
  })
  
  output$plot <- renderPlot({
    res <- plot_sankey(input$sex, input$hormone_class, input$hormone_name, input$hormone_tissue,
                       input$receptor_class, input$receptor_name, input$receptor_tissue,
                       input$number_HRtissue, input$tpm_log1p[1], input$tpm_log1p[2],
                       gtex_hr, df, vec_receptor)
    
    if (is.character(res)) {
      return(NULL)
    } else {
      return(res)
    }
  })
  
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("HOHC_data-", input$sex, ".tsv")
    },
    content = function(file) {
      write.table(filter_data(input$sex, input$hormone_class, input$hormone_name, input$hormone_tissue,
                              input$receptor_class, input$receptor_name, input$receptor_tissue,
                              input$number_HRtissue, input$tpm_log1p[1], input$tpm_log1p[2],
                              gtex_hr, df, vec_receptor), file, row.names = FALSE, sep = '\t', quote = FALSE)
    }
  )
}

# Run app ----
shinyApp(ui, server)