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


suppressMessages(library(shiny))
suppressMessages(library(shinythemes))
suppressMessages(library(markdown))

withConsoleRedirect <- function(containerId, expr) {

  txt <- capture.output(results <- expr, type = "message")
  if (length(txt) > 0) {
    insertUI(paste0("#", containerId), where = "beforeEnd",
             ui = paste0(txt, "\n", collapse = ""))
  }
  results

}

shinyApp(
    ui = fluidPage(
        shinyjs::useShinyjs(),
        theme = shinytheme("yeti"),
        sidebarLayout(
            position = "right",
            sidebarPanel(
              tags$head(tags$style(type = "text/css", "#loadmessage {
                position: fixed; top: 0px; left: 0px; width: 100%;
                padding: 5px 0px 5px 0px;
                text-align: center; font-weight: bold; font-size: 100%;
                color: #000000; background-color: #CCFF66;
                z-index: 105;
                }")),
              textInput("bam.files.dir",
                        "Full path to the directory with BAM files to be processed",
                        value = "/path/to/data",
                        width = "100%"),
              textInput("counts.path",
                        "If you have read counts that you want to re-use, specify full path to `counts.rds`",
                        value = "/path/to/counts.rds",
                        width = "100%"),
              selectInput("reference",
                          "Reference genome",
                          choices = c("hg38", "hg19"),
                          width = "50%"),
              selectInput("format",
                          "Format of input sequences",
                          choices = c("UCSC", "NCBI"),
                          width = "50%"),
              selectInput("skip.sorting.BAM",
                          "Skip sorting?",
                          choices = c("Yes", "No"),
                          width = "50%"),
              selectInput("skip.indexing.BAM",
                          "Skip indexing?",
                          choices = c("Yes", "No"),
                          width = "50%"),
              textInput("output.dir",
                        "Output directory",
                        value = getwd(),
                        width = "100%"),
              selectInput("bin.size.text",
                          "Bin size",
                          choices = c("200kb", "100kb", "50kb", "10kb", "1kb"),
                          width = "50%"),
              sliderInput("min.mapq",
                          "Minimum mapping quality of reads",
                          min = 0.1, max = 1.0,
                          value = 0.8, step = 0.1,
                          width = "100%"),
              sliderInput("rough.span",
                          "First (rough) Loess span",
                          min = 0.05, max = 0.50,
                          value = 0.10, step = 0.05,
                          width = "100%"),
              sliderInput("final.span",
                          "Second (final) Loess span",
                          min = 0.05, max = 0.50,
                          value = 0.30, step = 0.05,
                          width = "100%"),
              sliderInput("del.threshold",
                          "Threshold for deletions",
                          min = -2.50, max = -0.50,
                          value = -0.70, step = 0.01,
                          width = "100%"),
              sliderInput("dup.threshold",
                          "Threshold for duplications",
                          min = 0.30, max = 1.50,
                          value = 0.50, step = 0.01,
                          width = "100%"),
              selectInput("centromeres.margin",
                          "Centromeres margin",
                          choices = c("500kb", "100kb", "10kb", "0"),
                          width = "100%"),
              selectInput("seg.duplications.threshold",
                          "Lower threshold for segmental duplications",
                          choices = c("100kb", "50kb", "10kb", "0"),
                          width = "100%"),
              selectInput("verbose",
                          "Verbose",
                          choices = c("Yes", "No"),
                          width = "50%"),
              sliderInput("ncores", "Number of cores for parallelisation",
                          min = 4, max = 28,
                          value = 4, step = 4,
                          width = "100%"),
              selectInput("plot.results",
                          "Plot the results?",
                          choices = c("Yes", "No"),
                          width = "50%"),
              br(),
              actionButton("Run", "Run"),
              conditionalPanel(condition = "$('html').hasClass('shiny-busy')",
                               tags$div("Loading...", id = "loadmessage"))
            ),
            mainPanel(
                includeMarkdown("info.md"),
                h3("Console"), pre(id = "console"),
                downloadButton("download.table", "Variations"),
                downloadButton("download.plots", "Plots")
            )
        )
    ),
    server = function(input, output, session) {

      observeEvent(input$Run, {

        dict <- setNames(c(500e3, 200e3, 100e3, 50e3, 10e3, 1e3, 0),
                         c('500kb', '200kb', '100kb', '50kb', '10kb', '1kb', '0'))

        suppressMessages(library(GenomeInfoDb))
        suppressMessages(library(rCGH))

        working.dir <<- getwd()
        output.dir <<- input$output.dir
        bam.files.dir <<- input$bam.files.dir

        verbose <<- ifelse(input$verbose == 'Yes', TRUE, FALSE)
        plot.results <<- ifelse(input$plot.results == 'Yes', TRUE, FALSE)

        ncores <<- as.numeric(input$ncores)

        path.to.counts <<- input$counts.path

        skip.sorting.BAM <<- ifelse(input$skip.sorting.BAM == 'Yes', TRUE, FALSE)
        skip.indexing.BAM <<- ifelse(input$skip.indexing.BAM == 'Yes', TRUE, FALSE)

        bin.size.text <<- input$bin.size.text
        bin.size <<- as.numeric(dict[bin.size.text])

        reference <<- input$reference
        seqnames <- readRDS(paste0(reference, '.seqnames.rds'))
        assign(
          x = "format",
          value = input$format,
          envir = .GlobalEnv
        )
        seqlevelsStyle(seqnames) <- input$format
        assign(
          x = "seqnames",
          value = seqnames,
          envir = .GlobalEnv
        )

        centromeres <<- readRDS(paste0(reference, '.centromeres.rds'))
        syndromes <<- readRDS(paste0(reference, '.syndromes.rds'))

        if (input$format == 'UCSC'){
          autosomes <<- c(paste0('chr', 1:22))
        } else {
          autosomes <<- c(paste0('', 1:22))
        }

        min.mapq <<- input$min.mapq

        rough.span <<- input$rough.span
        final.span <<- input$final.span

        del.threshold <<- input$del.threshold
        dup.threshold <<- input$dup.threshold

        seg.duplications <<- readRDS(paste0(reference, '.segmental.duplications.rds'))

        centromeres.margin <<- as.numeric(dict[input$centromeres.margin])
        seg.duplications.threshold <<- as.numeric(dict[input$seg.duplications.threshold])

        withCallingHandlers(
          source("main.R"),
          message = function(m) {
            shinyjs::html("console", m$message, TRUE)
          }
        )
      })

      output$download.table <- downloadHandler(
        filename <- function() {
          paste0(output.dir, "/Variations.xlsx")
        },
        content <- function(file) {
          file.copy(paste0(output.dir, "/Variations.xlsx"), file)
        },
        contentType = "application/xlsx"
      )

      output$download.plots <- downloadHandler(
        filename <- function() {
          paste0(output.dir, "/Plots.zip")
        },
        content <- function(file) {
          zip(zipfile = paste0(output.dir, '/Plots.zip'), files = dir(paste0(output.dir, "/Plots"), full.names = TRUE))
          file.copy(paste0(output.dir, "/Plots.zip"), file)
        },
        contentType = "application/zip"
      )

    }
)
