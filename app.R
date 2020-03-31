# This file is part of ConanVarvar, a tool for detection of copy number variants
#
# Copyright (C) 2020 Victor Chang Cardiac Research Institute
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# This program is provided without warranty of any kind.
#
# <http://www.gnu.org/licenses/>


base::suppressMessages(base::library(shiny))
base::suppressMessages(base::library(shinyjs))
base::suppressMessages(base::library(shinydashboard))
base::suppressMessages(base::library(markdown))

withConsoleRedirect <- function(containerId, expr) {

  txt <- utils::capture.output(results <- expr, type = "message")
  if (base::length(txt) > 0) {
    shiny::insertUI(base::paste0("#", containerId), where = "beforeEnd",
                    ui = base::paste0(txt, "\n", collapse = ""))
  }
  results

}

shiny::shinyApp(
  ui = shinydashboard::dashboardPage(
    shinydashboard::dashboardHeader(disable = F, title = "ConanVarvar", titleWidth = 290),
    shinydashboard::dashboardSidebar(
      width = 290,
      tags$head(
        tags$style(
          shiny::HTML("
            .main-header .logo {
                font-family: \"Goudy Bookletter 1911\", sans-serif;
                font-weight: bold;
                font-size: 24px;
            }

            .js-irs-0 .irs-bar {
                border-top-color: #f04b5d;
                border-bottom-color: #f04b5d;
            }
            .js-irs-0 .irs-bar-edge {
                border-color: #f04b5d;
            }
            .js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {
                background: #f04b5d;
            }
            .js-irs-1 .irs-bar {
                border-top-color: #f04b5d;
                border-bottom-color: #f04b5d;
            }
            .js-irs-1 .irs-bar-edge {
                border-color: #f04b5d;
            }
            .js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {
                background: #f04b5d;
            }
            .js-irs-2 .irs-bar {
                border-top-color: #f04b5d;
                border-bottom-color: #f04b5d;
            }
            .js-irs-2 .irs-bar-edge {
                border-color: #f04b5d;
            }
            .js-irs-2 .irs-single, .js-irs-2 .irs-bar-edge, .js-irs-2 .irs-bar {
                background: #f04b5d;
            }
            .js-irs-3 .irs-bar {
                border-top-color: #f04b5d;
                border-bottom-color: #f04b5d;
            }
            .js-irs-3 .irs-bar-edge {
                border-color: #f04b5d;
            }
            .js-irs-3 .irs-single, .js-irs-3 .irs-bar-edge, .js-irs-3 .irs-bar {
                background: #f04b5d;
            }
            .js-irs-4 .irs-bar {
                border-top-color: #f04b5d;
                border-bottom-color: #f04b5d;
            }
            .js-irs-4 .irs-bar-edge {
                border-color: #f04b5d;
            }
            .js-irs-4 .irs-single, .js-irs-4 .irs-bar-edge, .js-irs-4 .irs-bar {
                background: #f04b5d;
            }
            .js-irs-5 .irs-bar {
                border-top-color: #f04b5d;
                border-bottom-color: #f04b5d;
            }
            .js-irs-5 .irs-bar-edge {
                border-color: #f04b5d;
            }
            .js-irs-5 .irs-single, .js-irs-5 .irs-bar-edge, .js-irs-5 .irs-bar {
                background: #f04b5d;
            }
            .js-irs-6 .irs-bar {
                border-top-color: #f04b5d;
                border-bottom-color: #f04b5d;
            }
            .js-irs-6 .irs-bar-edge {
                border-color: #f04b5d;
            }
            .js-irs-6 .irs-single, .js-irs-6 .irs-bar-edge, .js-irs-6 .irs-bar {
                background: #f04b5d;
            }

            .skin-blue .main-header .logo {
                background-color: #b11116;
            }
            .skin-blue .main-header .logo:hover {
                background-color: #b11116;
            }
            .skin-blue .main-header .navbar {
                background-color: #999899;
            }
            .skin-blue .main-sidebar {
                background-color: #999899;
                color: #b11116;
            }
            .skin-blue .sidebar-menu > li.active > a,
            .skin-blue .sidebar-menu > li:hover > a {
                border-left-color: #b11116;
            }
            .skin-blue .main-sidebar .sidebar .sidebar-menu .active a {
                background-color: #f04b5d;
            }
            .skin-blue .main-sidebar .sidebar .sidebar-menu a {
                background-color: #f04b5d;
                color: #000000;
            }
            .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover {
                background-color: #b11116;
            }
            .skin-blue .main-header .navbar .sidebar-toggle:hover {
                background-color: #b11116;
            }
          ")),
        tags$script(
          shiny::HTML("
            $(document).ready(function() {
                var ids = ['results', 'settings', 'help'];
                for (i = 0; i < ids.length; i++) {
                    $('a[data-value=' + ids[i] + ']').addClass('my_subitem_class');
                }
                $('.my_subitem_class').on('click',function() {
                    $('.my_subitem_class').parent().removeClass('active');
                })
            })
          "))
      ),
      shiny::actionButton("Run", "Run", width = 260),
      shiny::conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                              tags$div("Please, wait...", id = "loadmessage")),
      shinydashboard::sidebarMenu(
        shinydashboard::menuItem("Results", icon = shiny::icon("table"), tabName = "results"),
        shinydashboard::menuItem("Settings", icon = shiny::icon("cog", lib = "glyphicon"), tabName = "settings",
                                 shiny::textInput("bam.files.dir",
                                                  "Path to BAM files",
                                                  value = "/path/to/data",
                                                  width = "100%"),
                                 shiny::textInput("counts.path",
                                                  "Path to read counts (optional)",
                                                  value = "/path/to/counts.rds",
                                                  width = "100%"),
                                 shiny::selectInput("reference",
                                                    "Reference genome",
                                                    choices = base::c("hg38", "hg19"),
                                                    width = "100%"),
                                 shiny::selectInput("format",
                                                    "Format of input sequences",
                                                    choices = base::c("UCSC", "NCBI"),
                                                    width = "100%"),
                                 shiny::selectInput("skip.sorting.BAM",
                                                    "Skip sorting?",
                                                    choices = base::c("Yes", "No"),
                                                    width = "100%"),
                                 shiny::selectInput("skip.indexing.BAM",
                                                    "Skip indexing?",
                                                    choices = base::c("Yes", "No"),
                                                    width = "100%"),
                                 shiny::textInput("output.dir",
                                                  "Output directory",
                                                  value = "/output/",
                                                  width = "100%"),
                                 shiny::selectInput("bin.size.text",
                                                    "Bin size",
                                                    choices = base::c("200kb", "100kb", "50kb", "10kb", "1kb"),
                                                    selected = "50kb",
                                                    width = "50%"),
                                 shiny::selectInput("centromeres.margin",
                                                    "Margin for centromeres",
                                                    choices = base::c("500kb", "100kb", "10kb", "0"),
                                                    width = "100%"),
                                 shiny::selectInput("seg.duplications.threshold",
                                                    "Lower threshold for segmental duplications",
                                                    choices = base::c("100kb", "50kb", "10kb", "0"),
                                                    width = "100%"),
                                 shiny::sliderInput("ncores", "Number of cores for parallelisation",
                                                    min = 4, max = 28,
                                                    value = 4, step = 4,
                                                    width = "100%"),
                                 shiny::selectInput("plot.results",
                                                    "Plot the results?",
                                                    choices = base::c("Yes", "No"),
                                                    width = "100%"),
                                 shinydashboard::menuItem("Finer Settings", tabName = "finer_settings",
                                                          shiny::sliderInput("min.mapq",
                                                                             "Minimum mapping quality of reads",
                                                                             min = 0.1, max = 1.0,
                                                                             value = 0.8, step = 0.1,
                                                                             width = "100%"),
                                                          shiny::sliderInput("rough.span",
                                                                             "First (rough) Loess span",
                                                                             min = 0.05, max = 0.50,
                                                                             value = 0.10, step = 0.05,
                                                                             width = "100%"),
                                                          shiny::sliderInput("final.span",
                                                                             "Second (final) Loess span",
                                                                             min = 0.05, max = 0.50,
                                                                             value = 0.30, step = 0.05,
                                                                             width = "100%"),
                                                          shiny::sliderInput("occurrence.margin",
                                                                             "Start/End margin",
                                                                             min = 0, max = 10,
                                                                             value = 6, step = 1,
                                                                             width = "100%"),
                                                          shiny::sliderInput("del.threshold",
                                                                             "Threshold for deletions",
                                                                             min = -2.50, max = -0.50,
                                                                             value = -0.70, step = 0.01,
                                                                             width = "100%"),
                                                          shiny::sliderInput("dup.threshold",
                                                                             "Threshold for duplications",
                                                                             min = 0.30, max = 1.50,
                                                                             value = 0.50, step = 0.01,
                                                                             width = "100%"),
                                                          shiny::selectInput("verbose",
                                                                             "Verbose",
                                                                             choices = base::c("Yes", "No"),
                                                                             width = "50%")),
                                 shiny::br()
                                 ),
        shinydashboard::menuItem("Help", icon = shiny::icon("question"), tabName = "help"))
    ),
    shinydashboard::dashboardBody(
      tags$head(
        tags$style(shiny::HTML("
                     .main-sidebar {
                         font-size: 20px;
                     }
                     .sidebar-menu .treeview-menu {
                         padding: 0 0 0 50px;
                     }
                     .loadmessage {
                         position: fixed;
                         top: 0px;
                         left: 0px;
                         width: 100%;
                         padding: 0px 0px 0px 0px;
                         text-align: center;
                         font-weight: bold;
                         font-size: 100%;
                         color: #000000; background-color: #ef4640;
                     }
                   "))
      ),
      shinydashboard::tabItems(
        shinydashboard::tabItem(tabName = "results",
                                shiny::fluidRow(
                                  shiny::h3("Console"), shiny::pre(id = "console"),
                                  shiny::downloadButton("download.table", "Variations"),
                                  shiny::downloadButton("download.plots", "Plots")
                                )),
        shinydashboard::tabItem(tabName = "settings"),
        shinydashboard::tabItem(tabName = "help",
                                shiny::fluidRow(
                                  shiny::includeMarkdown("info.md")
                                ))
      )
    ),
    shinyjs::useShinyjs(),
    ),
  server = function(input, output, session) {

    shiny::observeEvent(input$Run, {

      dict <- stats::setNames(base::c(500e3, 200e3, 100e3, 50e3, 10e3, 1e3, 0),
                              base::c('500kb', '200kb', '100kb', '50kb', '10kb', '1kb', '0'))

      base::suppressMessages(base::library(GenomeInfoDb))
      base::suppressMessages(base::library(rCGH))

      working.dir <<- base::getwd()
      output.dir <<- input$output.dir
      bam.files.dir <<- input$bam.files.dir

      verbose <<- base::ifelse(input$verbose == 'Yes', TRUE, FALSE)
      plot.results <<- base::ifelse(input$plot.results == 'Yes', TRUE, FALSE)

      ncores <<- base::as.numeric(input$ncores)

      path.to.counts <<- input$counts.path

      skip.sorting.BAM <<- base::ifelse(input$skip.sorting.BAM == 'Yes', TRUE, FALSE)
      skip.indexing.BAM <<- base::ifelse(input$skip.indexing.BAM == 'Yes', TRUE, FALSE)

      bin.size.text <<- input$bin.size.text
      bin.size <<- base::as.numeric(dict[bin.size.text])

      reference <<- input$reference
      seqnames <- base::readRDS(base::paste0(reference, '.seqnames.rds'))
      base::assign(
        x = "format",
        value = input$format,
        envir = .GlobalEnv
      )
      GenomeInfoDb::seqlevelsStyle(seqnames) <- input$format
      base::assign(
        x = "seqnames",
        value = seqnames,
        envir = .GlobalEnv
      )

      centromeres <<- base::readRDS(base::paste0(reference, '.centromeres.rds'))
      syndromes <<- base::readRDS(base::paste0(reference, '.syndromes.rds'))

      if (input$format == 'UCSC'){
        autosomes <<- base::c(base::paste0('chr', 1:22))
      } else {
        autosomes <<- base::c(base::paste0('', 1:22))
      }

      min.mapq <<- input$min.mapq

      rough.span <<- input$rough.span
      final.span <<- input$final.span

      del.threshold <<- input$del.threshold
      dup.threshold <<- input$dup.threshold

      seg.duplications <<- base::readRDS(base::paste0(reference, '.segmental.duplications.rds'))

      centromeres.margin <<- base::as.numeric(dict[input$centromeres.margin])
      seg.duplications.threshold <<- base::as.numeric(dict[input$seg.duplications.threshold])

      occurrence.margin <<- input$occurrence.margin

      base::withCallingHandlers(
        base::source("main.R"),
        message = function(m) {
          shinyjs::html("console", m$message, TRUE)
        }
      )
    })

    output$download.table <- shiny::downloadHandler(
      filename <- function() {
        base::paste0(output.dir, "/Variations.xlsx")
      },
      content <- function(file) {
        base::file.copy(base::paste0(output.dir, "/Variations.xlsx"), file)
      },
      contentType = "application/xlsx"
    )

    output$download.plots <- shiny::downloadHandler(
      filename <- function() {
        base::paste0(output.dir, "/Plots.zip")
      },
      content <- function(file) {
        utils::zip(zipfile = base::paste0(output.dir, '/Plots.zip'), files = base::dir(base::paste0(output.dir, "/Plots"), full.names = TRUE))
        base::file.copy(base::paste0(output.dir, "/Plots.zip"), file)
      },
      contentType = "application/zip"
    )

  }
)
