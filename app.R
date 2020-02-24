# Shiny app to create timecourse panel plots.

library(magrittr)
library(data.table)
setDTthreads(1) # Do not allow data.table to run multiple threads.
library(ggplot2)
library(gridExtra)
library(DT)
library(viridis)

# To do list:
# XXX 1. Support entering gene names.
# XXX 2. Keep plot panels in the same order the gene names are entered in.
# XXX 3. Add title bar, info tab, instructions tab.


# Aux. to do list:
# XXX 4. Add tab for normalized count plots (alternative to rlog y-axis).
# XXX 5. Add DataTables-base gene searching (and maybe selecting).
# 6. Make number of panel columns selectable.
# 7. Display warning if not all entered genes are found.
# XXX 8. Move ensembl gene id to 'subtitle' slot to make sure it fits.
# 9. Display gene ids of genes that don't match any id's in the data.
# 10. New tab with option for heatmap (blue-scale) like Expression atlas (https://www.ebi.ac.uk/gxa/home)

# XXX 11. Make DataTable selection add-able to input list.
# XXX 12. Fix DataTable decimal representation.
# XXX 13. Add "Download PDF" button. 
#         Download PNG" button.
# 14. "Download Table" button to download diff expr table.
# 15. Feedback/Display which genes were selected. Show subset in DataTable??

# 16. Put info/description in separate page or section.

# 17. To customize Diff Expr Table title, redefine h4 and h5 styles.

# 18. Fix warning "Warning: Error in pngfun: invalid quartz() device size"
#     This occurs when you plot ~50 genes.

# Items to add to "About This Dataset" tab.
# https://b2b.hci.utah.edu/gnomex/gnomexFlex.jsp?requestNumber=468R
# genome build: GRCz10 ensembl release 89.

# Load data.table of rlog values for all heart timecourse genes and samples.
tab = fread("ribozero_all_values_normcount_rlog_tpm_20190514.txt")
tab[, gene_label:=paste(external_gene_name, " - ", ensembl_gene_id, sep="")]

example_input_heart = paste(c("ENSDARG00000019096", # myl7
                              "ENSDARG00000098952", # gata4
                              "ENSDARG00000018004", # nkx2.5
                              "ENSDARG00000092060", # tbx5b
                              "vmhc","myh6",
                              "bmp10","tnnt2a",
                              "cmlc1", "tbx18",
                              "mb", "postna"), collapse="\n")

# Load differential expression results table.
de_res = fread("20180612_ribozero_DESeqres_rep-time.txt")


index_tab = fread("gene_string_index_table.txt")

# To keep order of input gene ids, create a named vector (instead of data.table).
gene_string_index_vec = index_tab$ensembl_gene_id
names(gene_string_index_vec) = index_tab$gene_string

# Make a named vector (use as an ordered dictionary) to
# relate ensembl_gene_ids to gene_labels, and retaining order.
tmp_tab = tab[, list(ensembl_gene_id, gene_label)]
tmp_tab = tmp_tab[!duplicated(tmp_tab)]
gene_label_key = tmp_tab$gene_label
names(gene_label_key) = tmp_tab$ensembl_gene_id
rm(tmp_tab)

ui <- navbarPage("Zebrafish Heart Gene Expression Explorer",
  tabPanel("Visualization",
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
                      tags$p(h5("DESeq2 results for full dataset (n = 32250 genes)"))
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
  
      
  observeEvent(input$add_gene_ids, {
    updateTextAreaInput(session, "gene_id_input", 
                        value = paste(input$gene_id_input, "\n", new_gene_id_string(), sep="")
    )
    updateActionButton(session, inputId="update_plots")
  })
  
    
    # Reactive gene list. Split on any number of contiguous commas, spaces, tabs, newlines.
    gene_string_set <- eventReactive(input$update_plots, {
        strsplit(input$gene_id_input, split="[, \t\n]+")[[1]]
        
    }, ignoreNULL = FALSE)

    ens_gene_set <- reactive({
      res_list = list()
      gene_string_vec = gene_string_set()
      for (i in seq_along(gene_string_vec)) {
        tmp_gene_string = gene_string_vec[i]
        res_list[[i]] = index_tab[gene_string %in% tmp_gene_string]
      }
      res = rbindlist(res_list)
      res$ensembl_gene_id
    })

        # Reactive data subsetting.
    selected_tab <- reactive({
        tab[ensembl_gene_id %in% ens_gene_set()]
    })
 
    panel_y_labels = c("rlog_value"="Rlog value",
                       "norm_count"="Normalized count", 
                       "tpm"="Transcripts per million")
    
    plot_data <- reactive({
        plot_tab = selected_tab()
        

          means_tab = plot_tab[, list(mean_value=mean(get(input$expr_units)),
                                      ymax=max(get(input$expr_units)),
                                      ymin=min(get(input$expr_units))),
                           by=list(ensembl_gene_id, time_point, gene_label, external_gene_name)]
          
          
 
        # Make a table of centered z-scores for heatmap.
        
        #browser()
        
        #heatmap_means_tab = plot_tab[, list(mean_zscore=mean(zcount)), 
        #                      by=list(ensembl_gene_id, time_point, gene_label, external_gene_name)]
        #browser()
        
        panel_gene_ids = ens_gene_set()
        n_panels = length(panel_gene_ids)
        if (n_panels <= as.integer(input$plot_col_number)) {
          n_row = 1
          n_col = n_panels
        } else {
          n_row = ceiling(n_panels / as.integer(input$plot_col_number))
          n_col = as.integer(input$plot_col_number)
        }
        #browser()
        list(plot_tab=plot_tab,
             means_tab=means_tab,
             #heatmap_means_tab=heatmap_means_tab,
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

       
       
       pdf_file_name = paste("zebrafish_heart_gene_plots", 
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
    filename = paste("zebrafish_heart_selected_genes_raw_data", 
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

      heatmap_plot = ggplot(data=heatmap_tab, aes(x=time_point, y=gene_label, fill=zscore)) +
                     geom_tile(colour="white", size=1) +
                     #scale_fill_distiller(palette="Blues", direction=1) +
                     coord_cartesian(expand=FALSE) +
                     theme(axis.title.x=element_blank()) +
                     theme(axis.title.y=element_blank()) +
                     scale_fill_distiller(palette="RdBu", type="div", limit=heatmap_limits)
      
      # if(input$heatmap_color_scale == "blues_sequential") {
      #   heatmap_plot = heatmap_plot +
      #                  scale_fill_distiller(palette="Blues", direction=1)
      # } 
      
      
      
      
      # if(input$heatmap_color_scale == "red_blue_diverging") {
      #   heatmap_plot = heatmap_plot +
      #                  scale_fill_distiller(palette="RdBu", direction=1,
      #                                       values=seq(from=-2, to=2, length.out=7))
      #   
      # }
      
      # if(input$heatmap_color_scale == "red_blue_diverging") {
      #   heatmap_plot = heatmap_plot + 
      #                  scale_fill_gradientn(colours=c("#2166ac","#67a9cf","#d1e5f0","#f7f7f7",
      #                                                 "#fddbc7","#ef8a62","#b2182b"),
      #                                       values=seq(from=-2, to=2, length.out=7))
      # }
      grid::grid.draw(heatmap_plot)

      }, width=240 * 3, 
       height={pxh = 100 * sqrt(length(plot_data()$panel_gene_ids));
               if (pxh <= 560) {pxh} else {18 * length(plot_data()$panel_gene_ids)}
       },
       res=72)
  })
  

  #-------------------------------------------------

  
}

shinyApp(ui=ui, server=server)
















