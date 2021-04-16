# 19 Oct 2020.
# Devel version. Add miRNA panel plots and heatmaps.

# Shiny app to create timecourse panel plots.

library(magrittr)
library(data.table)
setDTthreads(1) # Do not allow data.table to run multiple threads.
library(ggplot2)
library(gridExtra)
library(DT)
library(viridis)

# library(reactlog) # For debugging. Disable for production.

# tell shiny to log all reactivity
# reactlog_enable()

# To do 27 Oct 2020.
# XXX 1. Fix miRNA heatmap plot display.

# XXX 3. Make mirbase id indexing case insensitive.
# 4. Add indexing for miRNA simplified names (e.g. find dre-miR-24b-5p given dre-miR-24).
# XXX 4a. ? Add message directing users to DataTable search box.
# 
# ! 5. Heatmap z-scores are incorrect. Need to compute z-scores using 
#    all 22 data points, not only the five timepoint-mean values.
# ! 6. Fix downloading buttons for data, png and pdf plot files.
# XXX 7. Check bug: 'Adding selected rows to gene set' doesn't work?
# 8. Add 'About This Dataset' text (Chelsea).
# 9. After fixes, push to github.
# 10. Update the live version on b2b.hci.utah.edu.

# To do list:
# XXX 1. Support entering gene names.
# XXX 2. Keep plot panels in the same order the gene names are entered in.
# XXX 3. Add title bar, info tab, instructions tab.

# Aux. to do list:
# XXX 4. Add tab for normalized count plots (alternative to rlog y-axis).
# XXX 5. Add DataTables-base gene searching (and maybe selecting).
# XXX 6. Make number of panel columns selectable.
# 7. Display warning if not all entered genes are found.
# XXX 8. Move ensembl gene id to 'subtitle' slot to make sure it fits.
# 9. Display gene ids of genes that don't match any id's in the data.
# 10. New tab with option for heatmap (blue-scale) like Expression atlas (https://www.ebi.ac.uk/gxa/home)

# XXX 11. Make DataTable selection add-able to input list.
# XXX 12. Fix DataTable decimal representation.
# XXX 13. Add "Download PDF" button. 
#         Download PNG" button.
# XXX 14. "Download Table" button to download diff expr table.
# 15. Feedback/Display which genes were selected. Show subset in DataTable??

# XXX 16. Put info/description in separate page or section.

# 17. To customize Diff Expr Table title, redefine h4 and h5 styles.

# 18. Fix warning "Warning: Error in pngfun: invalid quartz() device size"
#     This occurs when you plot ~50 genes.

# Items to add to "About This Dataset" tab.
# https://b2b.hci.utah.edu/gnomex/gnomexFlex.jsp?requestNumber=468R
# genome build: GRCz10 ensembl release 89.

# Load data.table of rlog values for all heart timecourse genes and samples.
tab = fread("ribozero_all_values_normcount_rlog_tpm_20190514.txt")
tab[, gene_label:=paste(external_gene_name, " - ", ensembl_gene_id, sep="")]

# Load mirna data.
mirna_tab = fread("mirna_datafiles/14806R_all_values_normcount_rlog_count_20200713.txt")
mirna_tab[, gene_label:=mirbase_id]

example_input_heart = paste(c("ENSDARG00000019096", # myl7
                              "ENSDARG00000098952", # gata4
                              "ENSDARG00000018004", # nkx2.5
                              "ENSDARG00000092060", # tbx5b
                              "vmhc","myh6",
                              "bmp10","tnnt2a",
                              "cmlc1", "tbx18",
                              "mb", "postna"), collapse="\n")

# To do: expand this list.
example_input_mirna = paste(c("dre-miR-499-5p",
                              "dre-miR-199-3p",
                              "dre-let-7a",
                              "dre-miR-99",
                              "dre-miR-430a-3p",
                              "dre-miR-125b-5p",
                              "dre-miR-21",
                              "dre-miR-196a-5p",
                              "dre-miR-143",
                              "dre-miR-1",
                              "dre-miR-7a",
                              "dre-miR-204-5p"), collapse="\n")
  
# Load differential expression results table.
de_res = fread("20180612_ribozero_DESeqres_rep-time.txt")

# mirna diff expr results.
de_res_mirna = fread("mirna_datafiles/14806R_deseq2_results_20200713.txt")


index_tab = fread("gene_string_index_table.txt")
# Similar table for mirna ids.
mirna_index_tab = fread("mirna_gene_string_table_temp.txt")

# To keep order of input gene ids, create a named vector (instead of data.table).
gene_string_index_vec = index_tab$ensembl_gene_id
names(gene_string_index_vec) = index_tab$gene_string

# mirna ids named vector.
mirna_string_index_vec = mirna_index_tab$mirbase_id
names(mirna_string_index_vec) = mirna_index_tab$gene_string

# Make a named vector (use as an ordered dictionary) to
# relate ensembl_gene_ids to gene_labels, and retaining order.
tmp_tab = tab[, list(ensembl_gene_id, gene_label)]
tmp_tab = tmp_tab[!duplicated(tmp_tab)]
gene_label_key = tmp_tab$gene_label
names(gene_label_key) = tmp_tab$ensembl_gene_id
rm(tmp_tab)

# mirna.
# Make a named vector (use as an ordered dictionary) to
# relate ensembl_gene_ids to gene_labels, and retaining order.
tmp_tab2 = mirna_tab[, list(mirbase_id, gene_label)]
tmp_tab2 = tmp_tab2[!duplicated(tmp_tab2)]
gene_label_key_mirna = tmp_tab2$gene_label
names(gene_label_key_mirna) = tmp_tab2$mirbase_id
rm(tmp_tab2)


ui <- navbarPage("Zebrafish Heart Gene Expression Explorer",
  tabPanel("Total RNA Visualization",
    fluidPage(
      tags$style(HTML("#diffexpr_table {border-top: 3px solid lightgrey;}")),
      tags$head(tags$style(HTML("#add_gene_ids {color: #fff;
                                                background-color: #337ab7;
                                                border-color: #2e6da4;}"))),
      tags$head(tags$style(HTML("#update_plots {color: #fff;
                                                background-color: #337ab7;
                                                border-color: #2e6da4;}"))),
      
      fluidRow(
        sidebarLayout(
          sidebarPanel(width=3,
                       textAreaInput(inputId="gene_id_input", 
                                     label="Input gene names or ensembl gene ids:",
                                     value=example_input_heart, rows=16),
                       br(),
                       actionButton(inputId="update_plots", "Update plots"),
                       br(),
                       br(),
                       selectInput(inputId="expr_units",
                                   label="Gene expression units:",
                                   selected="norm_count",
                                   choices=c("Rlog value"="rlog_value",
                                             "Normalized count"="norm_count", 
                                             "Transcripts per million"="tpm")),
                       
                       selectInput(inputId="plot_col_number",
                                   label="Number of plot columns",
                                   selected=4,
                                   choices=c("2"="2",
                                             "3"="3",
                                             "4"="4",
                                             "5"="5",
                                             "6"="6"))
            
                       # ,selectInput(inputId="heatmap_color_scale",
                       #             label="Select a heatmap color scale:",
                       #             choices=c("Blues"="blues_sequential", "Red-Blue"="red_blue_diverging"))
                       
          ),
          mainPanel(width=9,
            tabsetPanel(id="maintabset",
                        tabPanel(style = "overflow-y:scroll; min-height: 610px; position:relative;",
                                 title="Gene Panel Plots", 
                                 imageOutput(outputId="panel_plot")
                                 
                        ),
                        tabPanel(title="Heatmap Display",
                                 br(),
                                 imageOutput(outputId="heatmap_plot"),
                                 style = "overflow-y:scroll; min-height: 610px; position:relative;"
                        ),
                        tabPanel(title="Download Plots and Data",
                                 br(),
                                 downloadButton("download_pdf_plots", "Download panel plots (.pdf)"),
                                 br(),
                                 br(),
                                 downloadButton("download_png_plots", "Download panel plots (.png)"),
                                 br(),
                                 br(),
                                 downloadButton("download_current_point_data", "Download raw data for current gene list"),
                                 br(),
                                 br(),
                                 br(),
                                 br(),
                                 # downloadButton("download_full_raw_data", "Download full raw data table (.zip)"),
                                 # br(),
                                 # br(),
                                 downloadButton("download_full_deseq", "Download full DESeq2 results (.zip)")
                        )
            )
          )
        ),
        br()
      ),
      fluidRow(id="diffexpr_table", 
                
               column(width=4, 
                      tags$p(h4("Differential Expression Table")),
                      tags$p(h5(paste("DESeq2 results for all genes (n = ",
                                      length(unique(tab$ensembl_gene_id)),
                                      ")", sep="")))
               ),
               column(width=4, br(),
                      actionButton(inputId="add_gene_ids", "Add selected rows to gene list")
               )
      ),
      fluidRow(
               column(width=12, 
                      style = "min-height: 660px; position:relative;",
                      
                      DT::dataTableOutput(outputId = "geneIdTable")
               )
      )
      )
    ),
  tabPanel("miRNA Visualization",
    fluidPage(
      tags$style(HTML("#diffexpr_table {border-top: 3px solid lightgrey;}")),
      tags$head(tags$style(HTML("#add_gene_ids_mirna {color: #fff;
                                                background-color: #337ab7;
                                                border-color: #2e6da4;}"))),
      tags$head(tags$style(HTML("#update_plots_mirna {color: #fff;
                                                background-color: #337ab7;
                                                border-color: #2e6da4;}"))),
      
      fluidRow(
        sidebarLayout(
          sidebarPanel(width=3,
                       textAreaInput(inputId="gene_id_input_mirna", 
                                     label="Input gene names or ensembl gene ids:",
                                     value=example_input_mirna, rows=16),
                       br(),
                       actionButton(inputId="update_plots_mirna", "Update plots"),
                       br(),
                       br(),
                       selectInput(inputId="expr_units_mirna",
                                   label="Gene expression units:",
                                   selected="norm_count",
                                   choices=c("Rlog value"="rlog_value",
                                             "Normalized count"="norm_count")),
                       
                       selectInput(inputId="plot_col_number_mirna",
                                   label="Number of plot columns",
                                   selected=4,
                                   choices=c("2"="2",
                                             "3"="3",
                                             "4"="4",
                                             "5"="5",
                                             "6"="6")),
                       
                       tags$p(div(paste("Note: If your miR of interest does not appear,",
                                        "try searching in",
                                        "the Differential Expression Table below.",
                                        "Select the desired row(s),",
                                        "and then click 'Add selected rows to gene list'")))
            
                       # ,selectInput(inputId="heatmap_color_scale",
                       #             label="Select a heatmap color scale:",
                       #             choices=c("Blues"="blues_sequential", "Red-Blue"="red_blue_diverging"))
                       
          ),
          mainPanel(width=9,
            tabsetPanel(id="maintabset_mirna",
                        tabPanel(style = "overflow-y:scroll; min-height: 610px; position:relative;",
                                 title="Gene Panel Plots", 
                                 imageOutput(outputId="panel_plot_mirna")
                                 
                        ),
                        tabPanel(title="Heatmap Display",
                                 br(),
                                 imageOutput(outputId="heatmap_plot_mirna"),
                                 style = "overflow-y:scroll; min-height: 610px; position:relative;"
                        ),
                        tabPanel(title="Download Plots and Data",
                                 br(),
                                 downloadButton("download_pdf_plots_mirna", "Download miRNA panel plots (.pdf)"),
                                 br(),
                                 br(),
                                 downloadButton("download_png_plots_mirna", "Download miRNA panel plots (.png)"),
                                 br(),
                                 br(),
                                 downloadButton("download_current_point_data_mirna", "Download raw data for current miRNA list"),
                                 br(),
                                 br(),
                                 br(),
                                 br(),
                                 # downloadButton("download_full_raw_data", "Download full miRNA raw data table (.zip)"),
                                 # br(),
                                 # br(),
                                 downloadButton("download_full_deseq_mirna", "Download full miRNA DESeq2 results (.zip)")
                        )
            )
          )
        ),
        br()
      ),
      fluidRow(id="diffexpr_table_mirna", 
                
               column(width=4, 
                      tags$p(h4("Differential Expression Table")),
                      tags$p(h5(paste("DESeq2 results for all miRNAs (n = ", 
                                      length(unique(mirna_tab$mirbase_id)), 
                                      ")", sep="")))
               ),
               column(width=4, br(),
                      actionButton(inputId="add_gene_ids_mirna", "Add selected rows to gene list")
               )
      ),
      fluidRow(
               column(width=12, 
                      style = "min-height: 660px; position:relative;",
                      
                      DT::dataTableOutput(outputId = "geneIdTable_mirna")
               )
      )
      )
    ),
  tabPanel("About This Dataset", 
    fluidRow(
      column(width=3
      ),
      column(width=9,
             includeMarkdown("description.md")
      )
    )
  )
)

#browser()

server <- function(input, output, session) {
  
  new_gene_id_string = reactive({
    res = de_res[input$geneIdTable_rows_selected, external_gene_name]
    paste(res, collapse="\n")
  })
  
  # mirna.
  new_gene_id_string_mirna = reactive({
    res = de_res_mirna[input$geneIdTable_mirna_rows_selected, mirbase_id]
    paste(res, collapse="\n")
  })
      
  observeEvent(input$add_gene_ids, {
    updateTextAreaInput(session, "gene_id_input", 
                        value = paste(input$gene_id_input, "\n", new_gene_id_string(), sep="")
    )
    updateActionButton(session, inputId="update_plots")
  })
  
  # mirna.
  observeEvent(input$add_gene_ids_mirna, {
    updateTextAreaInput(session, "gene_id_input_mirna", 
                        value = paste(input$gene_id_input_mirna, "\n", new_gene_id_string_mirna(), sep="")
    )
    updateActionButton(session, inputId="update_plots_mirna")
  })
    
  # Reactive gene list. Split on any number of contiguous commas, spaces, tabs, newlines.
  gene_string_set <- eventReactive(input$update_plots, {
      strsplit(input$gene_id_input, split="[, \t\n]+")[[1]]
  }, ignoreNULL = FALSE)
  
  # mirna
  mirna_string_set <- eventReactive(input$update_plots_mirna, {
      strsplit(input$gene_id_input_mirna, split="[, \t\n]+")[[1]]
  }, ignoreNULL = FALSE)

  ens_gene_set <- reactive({
    res_list = list()
    gene_string_vec = gene_string_set()
    for (i in seq_along(gene_string_vec)) {
      tmp_gene_string = gene_string_vec[i]
      res_list[[i]] = index_tab[tolower(gene_string) %in% tolower(tmp_gene_string)]
    }
    res = rbindlist(res_list)
    res$ensembl_gene_id
  })
  
  # mirna
  mirna_set <- reactive({
    res_list = list()
    mirna_string_vec = mirna_string_set()
    for (i in seq_along(mirna_string_vec)) {
      tmp_gene_string = mirna_string_vec[i]
      res_list[[i]] = mirna_index_tab[tolower(gene_string) %in% tolower(tmp_gene_string)]
    }
    res = rbindlist(res_list)
    res$mirbase_id
  })
  
  
  # Reactive data subsetting.
  selected_tab <- reactive({
      tab[ensembl_gene_id %in% ens_gene_set()]
  })

  # mirna
  selected_tab_mirna <- reactive({
      mirna_tab[mirbase_id %in% mirna_set()]
  })
  
  panel_y_labels = c("rlog_value"="Rlog value",
                     "norm_count"="Normalized count", 
                     "tpm"="Transcripts per million")
  
  # mirna
  panel_y_labels_mirna = c("rlog_value"="Rlog value",
                           "norm_count"="Normalized count")
  
  plot_data <- reactive({
      plot_tab = selected_tab()

      means_tab = plot_tab[, list(mean_value=mean(get(input$expr_units)),
                                  ymax=max(get(input$expr_units)),
                                  ymin=min(get(input$expr_units))),
                       by=list(ensembl_gene_id, time_point, gene_label, external_gene_name)]

      panel_gene_ids = ens_gene_set()
      n_panels = length(panel_gene_ids)
      if (n_panels <= as.integer(input$plot_col_number)) {
        n_row = 1
        n_col = n_panels
      } else {
        n_row = ceiling(n_panels / as.integer(input$plot_col_number))
        n_col = as.integer(input$plot_col_number)
      }
      
      list(plot_tab=plot_tab,
           means_tab=means_tab,
           panel_gene_ids=panel_gene_ids,
           n_row=n_row,
           n_col=n_col)
  })
  
  # mirna
  plot_data_mirna <- reactive({
      plot_tab = selected_tab_mirna()

      means_tab = plot_tab[, list(mean_value=mean(get(input$expr_units_mirna)),
                                  ymax=max(get(input$expr_units_mirna)),
                                  ymin=min(get(input$expr_units_mirna))),
                       by=list(mirbase_id, time_point, gene_label)]

      panel_gene_ids = mirna_set()
      n_panels = length(panel_gene_ids)
      if (n_panels <= as.integer(input$plot_col_number_mirna)) {
        n_row = 1
        n_col = n_panels
      } else {
        n_row = ceiling(n_panels / as.integer(input$plot_col_number_mirna))
        n_col = as.integer(input$plot_col_number_mirna)
      }
      
      list(plot_tab=plot_tab,
           means_tab=means_tab,
           panel_gene_ids=panel_gene_ids,
           n_row=n_row,
           n_col=n_col)
  })

  
  # Reactive plot output.
  ggplot_output = reactive({
    
      plot_tab  = plot_data()$plot_tab
      means_tab = plot_data()$means_tab
      panel_gene_ids = plot_data()$panel_gene_ids
      
      panel_list = list()
      
      for (i in seq_along(panel_gene_ids)) {
        tmp_gene_id = panel_gene_ids[i]
        tmp_tab = plot_tab[ensembl_gene_id == tmp_gene_id]
        tmp_mean = means_tab[ensembl_gene_id == tmp_gene_id]
        panel_title = tmp_mean$external_gene_name[1]
        panel_subtitle = tmp_mean$ensembl_gene_id[1]
        
        tmp_plot = ggplot() +
                   geom_linerange(data=tmp_mean, 
                                  aes(x=time_point, ymin=ymin, ymax=ymax),
                                  color="grey50", size=0.6) +
                   geom_point(data=tmp_tab, 
                              aes(x=time_point, y=get(input$expr_units)),
                              size=1, colour="grey50") +
                   geom_line(data=tmp_mean, 
                             aes(x=time_point, y=mean_value, group=gene_label)) +
                   geom_point(data=tmp_mean,
                              aes(x=time_point, y=mean_value),
                              size=2, shape=21, fill="grey90") +
                   ylab(panel_y_labels[input$expr_units]) +
                   xlab("Timepoint") +
                   labs(title=panel_title, subtitle=panel_subtitle)
        
        panel_list[[i]] = tmp_plot
      }
      g = arrangeGrob(grobs=panel_list, 
                      ncol=plot_data()$n_col, 
                      nrow=plot_data()$n_row)
 
     
     
     pdf_file_name = paste("zebrafish_heart_gene_plots_", 
                            format(Sys.time(), "%Y%m%d_%H%M%S"), 
                            ".pdf", sep="")
     
     list(ggplot_obj=g,
          pdf_file_name=pdf_file_name)
  })

  # mirna
  ggplot_output_mirna = reactive({
    
      plot_tab  = plot_data_mirna()$plot_tab
      means_tab = plot_data_mirna()$means_tab
      panel_gene_ids = plot_data_mirna()$panel_gene_ids
      
      panel_list = list()
      
      for (i in seq_along(panel_gene_ids)) {
        tmp_gene_id = panel_gene_ids[i]
        tmp_tab = plot_tab[mirbase_id == tmp_gene_id]
        tmp_mean = means_tab[mirbase_id == tmp_gene_id]
        panel_title = tmp_mean$mirbase_id[1] # Replace with simplified miRNA names?
        panel_subtitle = "" # tmp_mean$ensembl_gene_id[1]
        
        tmp_plot = ggplot() +
                   geom_linerange(data=tmp_mean, 
                                  aes(x=time_point, ymin=ymin, ymax=ymax),
                                  color="grey50", size=0.6) +
                   geom_point(data=tmp_tab, 
                              aes(x=time_point, y=get(input$expr_units_mirna)),
                              size=1, colour="grey50") +
                   geom_line(data=tmp_mean, 
                             aes(x=time_point, y=mean_value, group=gene_label)) +
                   geom_point(data=tmp_mean,
                              aes(x=time_point, y=mean_value),
                              size=2, shape=21, fill="grey90") +
                   ylab(panel_y_labels_mirna[input$expr_units_mirna]) +
                   xlab("Timepoint") +
                   labs(title=panel_title, subtitle=panel_subtitle)
        
        panel_list[[i]] = tmp_plot
      }
      g = arrangeGrob(grobs=panel_list, 
                      ncol=plot_data_mirna()$n_col, 
                      nrow=plot_data_mirna()$n_row)
 
     
     
     pdf_file_name = paste("zebrafish_heart_miRNA_plots_", 
                            format(Sys.time(), "%Y%m%d_%H%M%S"), 
                            ".pdf", sep="")
     
     list(ggplot_obj=g,
          pdf_file_name=pdf_file_name)
  })

  
      
  observe({
    output$panel_plot <- renderPlot({
      grid::grid.draw(ggplot_output()$ggplot_obj)
    }, width=240 * plot_data()$n_col, height=200 * plot_data()$n_row, res=72)
  })
  
  # mirna
  observe({
    output$panel_plot_mirna <- renderPlot({
      grid::grid.draw(ggplot_output_mirna()$ggplot_obj)
    }, width=240 * plot_data_mirna()$n_col, height=200 * plot_data_mirna()$n_row, res=72)
  })

  
  # On b2b.hci.utah.edu, error presumably caused by ggsave lacking permissions.
  # "Warning: Error in <Anonymous>: cannot open file 'Rplots.pdf'"
  # Solution seems to be chown u0046369:shiny ./app_folder
  # followed by          chmod g+w ./app_folder
  output$download_pdf_plots = downloadHandler(
    filename = ggplot_output()$pdf_file_name,
    content = function(file) {
      ggsave(file, 
        plot=ggplot_output()$ggplot_obj, 
        width=(240 * plot_data()$n_col) / 72, 
        height=(200 * plot_data()$n_row) / 72)
    }
  )
  
  output$download_png_plots = downloadHandler(
    filename = gsub("pdf$", "png", ggplot_output()$pdf_file_name),
    content = function(file) {
      ggsave(file, 
             plot=ggplot_output()$ggplot_obj, 
             width=(240 * plot_data()$n_col) / 72, 
             height=(200 * plot_data()$n_row) / 72,
             units="in", dpi=)
    }
  )
  
  
  output$download_current_point_data = downloadHandler(
    filename = paste("zebrafish_heart_selected_genes_raw_data_", 
                     format(Sys.time(), "%Y%m%d_%H%M%S"), 
                     ".txt", sep=""),
    content = function(file) {
      fwrite(plot_data()$plot_tab, file=file, sep="\t")
    }
  )
  
  output$download_full_deseq = downloadHandler(
    filename = "zebrafish_heart_full_DESeq2_results.zip",
    content = function(file) {
      file.copy(from=file.path(".", "full_data_zips",
                               "zebrafish_heart_full_DESeq2_results.zip"),
                to=file)
    }
  )

  
  
  

  # mirna download handlers.
  output$download_pdf_plots_mirna = downloadHandler(
    filename = ggplot_output_mirna()$pdf_file_name,
    content = function(file) {
      ggsave(file, 
        plot=ggplot_output_mirna()$ggplot_obj, 
        width=(240 * plot_data()$n_col) / 72, 
        height=(200 * plot_data()$n_row) / 72)
    }
  )
  
  output$download_png_plots_mirna = downloadHandler(
    filename = gsub("pdf$", "png", ggplot_output_mirna()$pdf_file_name),
    content = function(file) {
      ggsave(file, 
             plot=ggplot_output_mirna()$ggplot_obj, 
             width=(240 * plot_data()$n_col) / 72, 
             height=(200 * plot_data()$n_row) / 72,
             units="in", dpi=)
    }
  )
  
  
  output$download_current_point_data_mirna = downloadHandler(
    filename = paste("zebrafish_heart_selected_mirnas_raw_data_", 
                     format(Sys.time(), "%Y%m%d_%H%M%S"), 
                     ".txt", sep=""),
    content = function(file) {
      fwrite(plot_data_mirna()$plot_tab, file=file, sep="\t")
    }
  )
  
  output$download_full_deseq_mirna = downloadHandler(
    filename = "zebrafish_heart_full_miRNA_DESeq2_results.zip",
    content = function(file) {
      file.copy(from=file.path(".", "full_data_zips",
                               "zebrafish_heart_full_miRNA_DESeq2_results.zip"),
                to=file)
    }
  )


  output$geneIdTable <- DT::renderDataTable(
    {
      DT::datatable(data = de_res[, list(external_gene_name,
                                         ensembl_gene_id,
                                         baseMean,
                                         log2FoldChange,
                                         pvalue,
                                         padj,
                                         chromosome_name,
                                         gene_biotype,
                                         `24hpf`,
                                         `36hpf`,
                                         `48hpf`,
                                         `60hpf`,
                                         `72hpf`)],
                    options = list(pageLength = 10),
                    rownames = FALSE,
                    filter="top",
                    colnames=c("\u00A0\u00A0Gene\u00A0name\u00A0\u00A0"="external_gene_name",
                               "Ensembl gene id"               ="ensembl_gene_id",
                               "Overall mean normalized counts"="baseMean",
                               "Log2 fold change"              ="log2FoldChange",
                               "P-value"                       ="pvalue",
                               "FDR adjusted p-value"          ="padj",
                               "Chromosome"                    ="chromosome_name",
                               "Gene biotype"                  ="gene_biotype",
                               "24hpf mean normalized counts"  ="24hpf",
                               "36hpf mean normalized counts"  ="36hpf",
                               "48hpf mean normalized counts"  ="48hpf",
                               "60hpf mean normalized counts"  ="60hpf",
                               "72hpf mean normalized counts"  ="72hpf")) %>%
        formatRound(columns=c("Overall mean normalized counts",
                              "24hpf mean normalized counts",
                              "36hpf mean normalized counts",
                              "48hpf mean normalized counts",
                              "60hpf mean normalized counts",
                              "72hpf mean normalized counts"), digits=2, mark="") %>%
        formatSignif(columns=c("Log2 fold change",
                               "P-value",
                               "FDR adjusted p-value"), digits=2)
    }, server=TRUE
  )

  # mirna
  output$geneIdTable_mirna <- DT::renderDataTable(
    {
      DT::datatable(data = de_res_mirna[, list(mirbase_id,
                                         # ensembl_gene_id,
                                         baseMean,
                                         log2FoldChange,
                                         pvalue,
                                         padj,
                                         # chromosome_name,
                                         # gene_biotype,
                                         `24hpf`,
                                         `36hpf`,
                                         `48hpf`,
                                         `60hpf`,
                                         `72hpf`)],
                    options = list(pageLength = 10),
                    rownames = FALSE,
                    filter="top",
                    colnames=c("\u00A0\u00A0miRNA\u00A0name\u00A0\u00A0"="mirbase_id",
                               # "Ensembl gene id"               ="ensembl_gene_id",
                               "Overall mean normalized counts"="baseMean",
                               "Log2 fold change"              ="log2FoldChange",
                               "P-value"                       ="pvalue",
                               "FDR adjusted p-value"          ="padj",
                               # "Chromosome"                    ="chromosome_name",
                               # "Gene biotype"                  ="gene_biotype",
                               "24hpf mean normalized counts"  ="24hpf",
                               "36hpf mean normalized counts"  ="36hpf",
                               "48hpf mean normalized counts"  ="48hpf",
                               "60hpf mean normalized counts"  ="60hpf",
                               "72hpf mean normalized counts"  ="72hpf")) %>%
        formatRound(columns=c("Overall mean normalized counts",
                              "24hpf mean normalized counts",
                              "36hpf mean normalized counts",
                              "48hpf mean normalized counts",
                              "60hpf mean normalized counts",
                              "72hpf mean normalized counts"), digits=2, mark="") %>%
        formatSignif(columns=c("Log2 fold change",
                               "P-value",
                               "FDR adjusted p-value"), digits=2)
    }, server=TRUE
  )

  #-------------------------------------------------
  # Heatmap plotting.
  
  observe({
    output$heatmap_plot <- renderPlot({
      
      plot_tab  = plot_data()$plot_tab
      
      heatmap_tab = plot_tab[, list(mean_norm_count=mean(norm_count)),
                             by=list(ensembl_gene_id, external_gene_name, gene_label, time_point)]
      heatmap_tab[, zscore:=(mean_norm_count - mean(mean_norm_count)) / sd(mean_norm_count),
                  by=list(ensembl_gene_id, external_gene_name, gene_label)]
      
      panel_gene_ids = plot_data()$panel_gene_ids
      
      # panel_gene_ids cannot contain duplicates for the heatmap plot.
      panel_gene_ids = panel_gene_ids[!duplicated(panel_gene_ids)]
      
      # Using the gene_label_key, create y-axis labels in the correct
      # order.
      heatmap_tab[, gene_label:=factor(as.character(gene_label), 
                                    levels=rev(gene_label_key[panel_gene_ids]))]

      heatmap_limits <- max(abs(heatmap_tab$zscore)) * c(-1, 1)
      
      # browser()
      
      heatmap_plot = ggplot(data=heatmap_tab, aes(x=time_point, y=gene_label, fill=zscore)) +
                     geom_tile(colour="white", size=1) +
                     #scale_fill_distiller(palette="Blues", direction=1) +
                     coord_cartesian(expand=FALSE) +
                     theme(axis.title.x=element_blank()) +
                     theme(axis.title.y=element_blank()) +
                     scale_fill_distiller(palette="RdBu", type="div", limit=heatmap_limits)
      
      grid::grid.draw(heatmap_plot)

      }, width=240 * 3, 
       height={pxh = 100 * sqrt(length(plot_data()$panel_gene_ids));
               if (pxh <= 560) {pxh} else {18 * length(plot_data()$panel_gene_ids)}
       },
       res=72)
  })

  
  # Chelsea heatmap code from manuscript Fig 3.
  # hm3 = ggplot(summary_signorm, aes(x = time_point, y = mirge_mir_id2, fill = zcount)) +
  # theme_bw() +
  # geom_tile() +
  # scale_fill_distiller(palette = "BrBG", name="z-score") +
  # #scale_x_discrete(labels = names(time_vec)) +
  # ylab("miRNA") +
  # xlab("Time point") +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10)) +
  # theme(axis.ticks = element_blank())
  
  # mirna
  observe({
    output$heatmap_plot_mirna <- renderPlot({
      
      plot_tab  = plot_data_mirna()$plot_tab
      
      heatmap_tab = plot_tab[, list(mean_norm_count=mean(norm_count)),
                             by=list(mirbase_id, gene_label, time_point)]
      heatmap_tab[, zscore:=(mean_norm_count - mean(mean_norm_count)) / sd(mean_norm_count),
                  by=list(mirbase_id, gene_label)]
      
      panel_gene_ids = plot_data_mirna()$panel_gene_ids
      
      # panel_gene_ids cannot contain duplicates for the heatmap plot.
      panel_gene_ids = panel_gene_ids[!duplicated(panel_gene_ids)]
      
      # Using the gene_label_key, create y-axis labels in the correct
      # order.
      heatmap_tab[, gene_label:=factor(as.character(gene_label), 
                                    levels=rev(gene_label_key_mirna[panel_gene_ids]))]

      heatmap_limits <- max(abs(heatmap_tab$zscore)) * c(-1, 1)
      
      #browser()
      
      heatmap_plot = ggplot(data=heatmap_tab, aes(x=time_point, y=gene_label, fill=zscore)) +
                     geom_tile(colour="white", size=1) +
                     #scale_fill_distiller(palette="Blues", direction=1) +
                     coord_cartesian(expand=FALSE) +
                     theme(axis.title.x=element_blank()) +
                     theme(axis.title.y=element_blank()) +
                     scale_fill_distiller(palette="RdBu", type="div", limit=heatmap_limits)
      
      grid::grid.draw(heatmap_plot)

      }, width=240 * 3, 
       height={pxh = 100 * sqrt(length(plot_data_mirna()$panel_gene_ids));
               if (pxh <= 560) {pxh} else {18 * length(plot_data_mirna()$panel_gene_ids)}
       },
       res=72)
  })

  
  

  #-------------------------------------------------

  
}

shinyApp(ui=ui, server=server)
















