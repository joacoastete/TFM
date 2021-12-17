library(shiny)
library(myvariant)
library(biomaRt)
library(knitr)
library(purrr)
library(dplyr, warn.conflicts = FALSE)

source("/Users/Joaco/Desktop/Joaco/Medicina/UOC/TFM/TFM/load_files.R")
source("/Users/Joaco/Desktop/Joaco/Medicina/UOC/TFM/TFM/classiServer.R")

#rmsk

    ui <- fluidPage(
        theme = bslib::bs_theme(bootswatch = "lumen"),
        titlePanel("Genetic Variant Classication"),
        sidebarLayout(
            sidebarPanel(
            tabsetPanel(
            tabPanel("Import data",
                     textInput("variant", "HGVS","chr1:g.100380997del"),
                     numericInput("mis.cut", "Missense  Z score cut-off (default=3.09)",3.09,-10,10),
                     numericInput("af.cut", "Allele frequency expected for the disorder (default=0.01)",0.01,0,1),
                     selectInput("PS2", div(HTML("If the variant is <em>de novo</em> (parental samples test negative):")),
                                 c("None selected","Paternity confirmed", "Paternity non confirmed")),
                     selectInput("PS3", "There are functional studies about the variant:",
                                  c("None selected","Studies supportive of a damaging effect", "Studies shows no damaging effect")),
                     selectInput("PM3", div(HTML("For recessive disorders, detected in <em>trans</em> with a pathogenic variant:")),
                                  c("Yes", "No"), selected = c("No")),
                     selectInput("BP2", div(HTML("Observed in <em>trans</em> with a pathogenic variant for a fully penetrant dominant disorder; 
                                                  or observed in <em>cis</em>with a pathogenic variant in anyinheritance pattern")),
                                  c("Yes", "No"), selected = c("No")),
                     selectInput("PP1", "Segregation analysis:",
                                 c("None selected","Co-segregation with disease in multiple affected family members", "Lack of segregation in affected members of a family")),
                     selectInput("PP4","Patient’s phenotype or family history is highly specific for a disease with asingle genetic etiology:",
                                 c("Yes", "No"), selected = c("No")),
                     selectInput("BP5","Variant found in a case with an alternate molecular basis for disease:",
                                 c("Yes", "No"), selected = c("No")),
                     selectInput("BP7","A synonymous variant for which splicing prediction algorithms predict no impact to the splice consensus sequence and the nucleotide is not highly conserved",
                                 c("Yes", "No"), selected = c("No")),
                     actionButton("button1","Submit")),
            tabPanel("af calculator",checkboxGroupInput("prueba", "test", c(3,2,5)), textOutput("test")))),
        mainPanel(fluidRow(
                     uiOutput("manual")),
            fluidRow(tableOutput("data"))))
        )


server <- function(input, output) {
       df<-eventReactive(input$button1,{first.table(input$variant,input$mis.cut,input$af.cut,input$PS2,input$PS3,
                                                input$PM3,input$BP2,input$PP1,input$PP4,input$BP5,input$BP7)})
    #output$test<-renderText({paste0(str(input$patho))})
    output$manual<-renderUI({fluidRow(column(1),column(5,
                                             checkboxGroupInput("patho", "pathogenic", group.pat, selected=check.pat(df())),
                                             actionButton("button2","Submit")),
                                      column(5,checkboxGroupInput("benign", "benign", group.ben, selected=check.ben(df()))))})
   tab.final<-eventReactive(input$button2,{classification.func(df(),check.pat(df()),check.ben(df()),input$patho,input$benign)})
   output$data<-renderTable({tab.final()})
    #output$data <- renderTable({first.table(input$variant,input$mis.cut,input$af.cut,input$PS2,input$PS3,
                                           # input$PM3,input$BP2,input$PP1,input$PP4,input$BP5,input$BP7)})
}

# Run the application 
shinyApp(ui = ui, server = server)
