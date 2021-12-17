library(shiny)
library(myvariant)
library(BiocManager)
options(repos = BiocManager::repositories())
library(biomaRt)
library(knitr)
library(dplyr)

# Load the R files with the functions and the required files
source("load_files.R")
source("classiServer.R")

ui <- fluidPage(
    theme = bslib::bs_theme(bootswatch = "lumen"),
    titlePanel(title=div("Genetic Variant Classication",img(src = "logo-uocs.png"), align="center")),
    sidebarLayout(
        sidebarPanel(width=5,
                     tabsetPanel(
                        
                         # Tab for the variant ID and the manuals inputs
                          tabPanel("Variant classification",
                                  textInput("variant", strong("HGVS id or rsID"),"chr1:g.100380997del"),
                                  numericInput("mis.cut", "Missense  Z score cut-off (default=3.09)",3.09,-10,10),
                                  numericInput("af.cut", "Allele frequency expected for the disorder (default=0.01)",0.01,0,1),
                                  selectInput("PS2", div(HTML("If the variant is <em>de novo</em> (parental samples test negative):")),
                                              c("None selected","Paternity confirmed", "Paternity non confirmed"),selected="None selected"),
                                  selectInput("PS3", "There are functional studies about the variant:",
                                              c("None selected","Studies supportive of a damaging effect", "Studies shows no damaging effect"),selected="None selected"),
                                  selectInput("PM3", div(HTML("For recessive disorders, detected in <em>trans</em> with a pathogenic variant:")),
                                              c("Yes", "No"), selected = c("No")),
                                              selectInput("BP2", div(HTML("Observed in <em>trans</em> with a pathogenic variant for a fully penetrant dominant disorder; 
                                                  or observed in <em>cis</em> with a pathogenic variant in any inheritance pattern")),
                                                          c("Yes", "No"), selected = c("No")),
                                  selectInput("PP1", "Segregation analysis:",
                                              c("None selected","Co-segregation with disease in multiple affected family members", 
                                                "Lack of segregation in affected members of a family"),selected="None selected"),
                                  selectInput("PP4","Patient’s phenotype or family history is highly specific for a disease with a single genetic etiology:",
                                              c("Yes", "No"), selected = c("No")),
                                  selectInput("BP5","Variant found in a case with an alternate molecular basis for disease:",
                                              c("Yes", "No"), selected = c("No")),
                                  selectInput("BP7","A synonymous variant for which splicing prediction algorithms predict no impact 
                                              to the splice consensus sequence and the nucleotide is not highly conserved",
                                              c("Yes", "No"), selected = c("No")),
                                  actionButton("button1","Submit"),
                                  actionButton("reset","Reset")),
                        
                           # A tab for the AF calculator
                          tabPanel("AF calculator",
                                  h3("Maximum credible population allele frequency calculator",align="center"),
                                  numericInput("prevalence", "Prevalence = 1 in ... (people)",0, min=0),
                                  numericInput("heterogeneity", "Heterogeneity (0-1):",0, min=0, max=1, step=0.05),
                                  numericInput("penetrance", "Penetrance (0-1):",0, min=0, max=1,step=0.05),
                                  textOutput("calculator")))),
        
       # The main panel will show the  calculated criteria and the final classification
         mainPanel(
                  id = "mainpanel", width=7,
                  fluidRow(
                     uiOutput("manual")),
            fluidRow(dataTableOutput("data"))))
        )


server <- function(input, output, session) {
    
    # Apply the function to calculate the criteria for the variant 
    df<-eventReactive(input$button1,{
        validate(need(input$variant != "", "Please enter an ID"))
        first.table(input$variant,input$mis.cut,input$af.cut,input$PS2,input$PS3,
                             input$PM3,input$BP2,input$PP1,input$PP4,input$BP5,input$BP7)
        })  
    
    # Render a UI with the caculated criteria and gives the otion to change it
    observeEvent(input$button1,
                 output$manual<-renderUI({
                     fluidRow(column(1),column(3,
                                               checkboxGroupInput("patho", "Evidence of pathogenicity", group.pat, selected=check.pat(df())),
                                               actionButton("button2","Submit")),
                              column(3,checkboxGroupInput("benign", "Evidence of benign impact", group.ben, selected=check.ben(df()))))
               }))
    
    # Generate the final table with the classification and other data
    observeEvent(input$button2,{
        tab.final<-classification.func(df(),check.pat(df()),check.ben(df()),input$patho,input$benign)
        output$data<-renderDataTable({
            if( isTRUE(is.null(input$patho) && is.null(input$benign))){
              validate("Please select a criteria")
                }
            tab.final},escape=F)})
    # Apply the AF calculator formula and print the result
    output$calculator<-renderText({
        calc<-af.calc(input$prevalence, input$heterogeneity, input$penetrance)
        paste0("Maximum credible population allele frequency: ",calc)})
    
    # Reset all the values and data
    observeEvent(input$reset, {
        output$manual<-renderUI(NULL)
        output$data<-renderDataTable(NULL)
        updateTextInput(session, "variant", value = "")
        updateNumericInput(session,"mis.cut", value = 3.09)
        updateNumericInput(session,"af.cut", value = 0.01)
        updateSelectInput(session, "PS2", selected="None selected")
        updateSelectInput(session, "PS3", selected="None selected")
        updateSelectInput(session, "BP2", selected="No")
        updateSelectInput(session, "PP1", selected="None selected")
        updateSelectInput(session, "PP4", selected="No")
        updateSelectInput(session, "BP5", selected="No")
        updateSelectInput(session, "BP7", selected="No")
   })
}

#Run the application 
shinyApp(ui = ui, server = server)
