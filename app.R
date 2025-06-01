# app.R

library(shiny)
library(shinythemes)
library(ggplot2)
library(plotly)
library(dplyr)
library(tibble)
library(readr)
library(tidyr)      # for pivot_longer
library(vegan)
library(pheatmap)
library(viridis)
library(ape)

# ---------------------------------------------------------------------------
# Helper function: clean sample names
# ---------------------------------------------------------------------------
fix_sample_names <- function(sample_names) {
  sample_names <- gsub("[-_/]", ".", sample_names)
  sample_names <- gsub("[\\[\\]\\^\\$]", "", sample_names)
  return(sample_names)
}

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  theme = shinytheme("flatly"),
  fluidRow(
    column(1, img(src = "hutch.png", height = "60px", width = "60px")),
    column(11, h2("Microbiome ITS & 16S analysis"))
  ),
  sidebarLayout(
    sidebarPanel(
      actionButton("show", "Upload files", width = "100%"),
      tags$br(), tags$br(),
      textOutput("selected_meta"),
      textOutput("selected_abund"),
      tags$hr(),
      selectInput(
        "analysis", "Select analysis:",
        choices = c(
          "Data",
          "Alpha Diversity",
          "Stacked Bar",
          "Reads per Sample",
          "Heatmap",
          "PCoA (Bray-Curtis)"
        )
      ),
      
      # Controls for Alpha Diversity
      conditionalPanel(
        "input.analysis=='Alpha Diversity'",
        selectInput("alpha_group", "Group by:", choices = NULL),
        selectInput(
          "alpha_index", "Diversity index:",
          choices = c(
            "Observed"  = "observed",
            "Shannon"   = "shannon",
            "Simpson"   = "simpson",
            "Evenness"  = "evenness"
          ),
          selected = "observed"
        )
      ),
      
      # Controls for Stacked Bar
      conditionalPanel(
        "input.analysis=='Stacked Bar'",
        sliderInput("num_stacks", "Number of taxa:", min = 2, max = 40, value = 10),
        selectInput("stack_group", "Divide plot by:", choices = NULL, selected = "None")
      ),
      
      # Controls for PCoA
      conditionalPanel(
        "input.analysis=='PCoA (Bray-Curtis)'",
        selectInput("pcoa_group", "Color/Group by:", choices = NULL, selected = NULL)
      ),
      
      width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Data",
          h4("Raw Taxon Table (head)"), tableOutput("head_raw_abund"),
          h4("Raw Metadata (head)"), tableOutput("head_raw_meta")
        ),
        tabPanel(
          "Alpha Diversity",
          plotlyOutput("alpha_boxplot"),
          downloadButton("downloadAlpha", "Download Alpha Table")
        ),
        tabPanel("Stacked Bar", plotlyOutput("stacked_bar")),
        tabPanel("Reads per Sample", plotlyOutput("reads_per_sample")),
        tabPanel("Heatmap", plotOutput("heatmap")),
        tabPanel("PCoA (Bray-Curtis)", plotOutput("pcoa"))
      ),
      width = 9
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {
  options(warn = -1)
  
  # Modal for file upload
  dataModal <- function(failed = FALSE) {
    modalDialog(
      title = "Upload metadata CSV and abundance CSV",
      fileInput("meta_file", "Choose metadata CSV", accept = ".csv"),
      fileInput("abund_file", "Choose abundance CSV", accept = ".csv"),
      if (failed) div(tags$b("Please upload both files", style = "color:red;")),
      footer = tagList(modalButton("Cancel"), actionButton("ok", "OK"))
    )
  }
  observeEvent(input$show, showModal(dataModal()))
  observeEvent(input$ok, {
    if (!is.null(input$meta_file) && !is.null(input$abund_file)) {
      removeModal()
    } else {
      showModal(dataModal(failed = TRUE))
    }
  })
  
  # 1) Read uploaded files
  raw_meta <- reactive({
    req(input$meta_file)
    df <- read_csv(input$meta_file$datapath, col_types = cols())
    if (ncol(df) < 2) {
      stop("Metadata must have at least two columns; first is sample ID.")
    }
    df
  })
  raw_abund <- reactive({
    req(input$abund_file)
    df <- read_csv(input$abund_file$datapath, col_types = cols())
    if (!"tax_name" %in% colnames(df)) {
      stop("Abundance file must have a 'tax_name' column first.")
    }
    df
  })
  
  # 2) Show heads and filenames
  output$head_raw_meta <- renderTable({ head(raw_meta(), 5) })
  output$head_raw_abund <- renderTable({ head(raw_abund(), 5) })
  output$selected_meta <- renderText({
    req(input$meta_file)
    paste("Metadata:", input$meta_file$name)
  })
  output$selected_abund <- renderText({
    req(input$abund_file)
    paste("Abundance table:", input$abund_file$name)
  })
  
  # 3) Prepare meta_df() and abund_df()
  meta_df <- reactive({
    df <- raw_meta()
    first_col <- colnames(df)[1]
    df <- df %>% rename(SampleID = !!first_col)
    df$SampleID <- fix_sample_names(df$SampleID)
    df
  })
  abund_df <- reactive({
    df <- raw_abund()
    keep_cols <- colnames(df)[!(colnames(df) %in% c("tax_id", "rank"))]
    df <- df[keep_cols]
    df2 <- df %>% column_to_rownames(var = "tax_name")
    df2[] <- lapply(df2, as.numeric)
    colnames(df2) <- fix_sample_names(colnames(df2))
    df2
  })
  
  # 4) Update “Group by” dropdowns
  observe({
    md <- meta_df()
    vars <- setdiff(colnames(md), "SampleID")
    updateSelectInput(session, "alpha_group", choices = vars, selected = vars[1])
    updateSelectInput(session, "stack_group", choices = c("None", vars), selected = "None")
    updateSelectInput(session, "pcoa_group", choices = vars, selected = vars[1])
  })
  
  # 5) Compute alpha diversity per sample
  alpha_df <- reactive({
    req(abund_df(), meta_df(), input$alpha_group, input$alpha_index)
    df_ab <- abund_df()
    df_meta <- meta_df()
    samples <- colnames(df_ab)
    if (!all(samples %in% df_meta$SampleID)) {
      stop("After cleaning, abundance samples do not match metadata SampleID.")
    }
    idx <- input$alpha_index
    values <- sapply(samples, function(s) {
      vec <- df_ab[[s]]
      if (idx == "observed") {
        sum(vec > 0, na.rm = TRUE)
      } else if (idx == "shannon") {
        vegan::diversity(vec, index = "shannon")
      } else if (idx == "simpson") {
        vegan::diversity(vec, index = "simpson")
      } else if (idx == "evenness") {
        H <- vegan::diversity(vec, index = "shannon")
        R <- sum(vec > 0, na.rm = TRUE)
        if (R > 0) H / log(R) else NA
      }
    })
    alpha_tab <- tibble(
      SampleID = samples,
      DiversityValue = as.numeric(values)
    )
    left_join(alpha_tab, df_meta, by = "SampleID")
  })
  
  # 6) Plot alpha diversity
  output$alpha_boxplot <- renderPlotly({
    req(input$analysis == "Alpha Diversity", alpha_df())
    df <- alpha_df()
    grp <- input$alpha_group
    df[[grp]] <- as.factor(df[[grp]])
    p <- ggplot(df, aes_string(x = grp, y = "DiversityValue", fill = grp, group = grp)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(aes_string(color = grp), width = 0.1, size = 2, alpha = 0.7) +
      labs(
        x = grp,
        y = paste("Value of", input$alpha_index),
        title = paste0("Alpha diversity (", input$alpha_index, ") by ", grp)
      ) +
      theme_bw() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
      )
    ggplotly(p, tooltip = c("y"))
  })
  output$downloadAlpha <- downloadHandler(
    filename = function() {
      paste0("alpha_diversity_", input$alpha_index, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write_csv(alpha_df(), file)
    }
  )
  
  # 7) Stacked Bar: one bar per sample, faceted if grouped
  output$stacked_bar <- renderPlotly({
    req(input$analysis == "Stacked Bar", abund_df(), meta_df())
    df_ab <- abund_df()
    df_md <- meta_df()
    grp <- input$stack_group
    n_tax <- input$num_stacks
    if (!all(colnames(df_ab) %in% df_md$SampleID)) {
      stop("After cleaning, abundance colnames do not match metadata SampleID.")
    }
    long0 <- df_ab %>%
      rownames_to_column("Taxon") %>%
      pivot_longer(
        cols = -Taxon,
        names_to = "SampleID",
        values_to = "Count"
      ) %>%
      left_join(df_md, by = "SampleID")
    total_per_sample <- long0 %>%
      group_by(SampleID) %>%
      summarise(TotalSample = sum(Count, na.rm = TRUE), .groups = "drop")
    long1 <- long0 %>%
      left_join(total_per_sample, by = "SampleID") %>%
      mutate(RelAbundance = ifelse(TotalSample > 0, Count / TotalSample, 0))
    total_ab_catalog <- long1 %>%
      group_by(Taxon) %>%
      summarise(TotalAcross = sum(RelAbundance, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(TotalAcross))
    top_taxa <- head(total_ab_catalog$Taxon, n_tax)
    if (grp == "None") {
      long2 <- long1 %>%
        mutate(Taxon2 = ifelse(Taxon %in% top_taxa, Taxon, "Other")) %>%
        group_by(SampleID, Taxon2) %>%
        summarise(RelAbund2 = sum(RelAbundance, na.rm = TRUE), .groups = "drop")
    } else {
      long2 <- long1 %>%
        mutate(Taxon2 = ifelse(Taxon %in% top_taxa, Taxon, "Other")) %>%
        group_by(SampleID, Taxon2, .data[[grp]]) %>%
        summarise(RelAbund2 = sum(RelAbundance, na.rm = TRUE), .groups = "drop")
      long2[[grp]] <- as.factor(long2[[grp]])
    }
    taxon_levels <- c(top_taxa, "Other")
    pal <- viridis::viridis(length(taxon_levels))
    colors_named <- setNames(pal, taxon_levels)
    if (grp == "None") {
      p <- ggplot(long2, aes(x = SampleID, y = RelAbund2, fill = Taxon2)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = colors_named) +
        labs(x = "Sample", y = "Relative Abundance", title = "Stacked bar (per sample)") +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.title = element_blank()
        )
    } else {
      p <- ggplot(long2, aes(x = SampleID, y = RelAbund2, fill = Taxon2)) +
        geom_bar(stat = "identity") +
        scale_fill_manual(values = colors_named) +
        facet_wrap(as.formula(paste("~", grp)), scales = "free_x", nrow = 1) +
        labs(x = grp, y = "Relative Abundance", title = paste0("Stacked bar (samples split by ", grp, ")")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.title = element_blank(),
          strip.text.x = element_text(size = 12)
        )
    }
    ggplotly(p)
  })
  
  # 8) Reads per Sample
  output$reads_per_sample <- renderPlotly({
    req(input$analysis == "Reads per Sample", abund_df())
    df <- abund_df()
    tot_reads <- tibble(
      Sample = colnames(df),
      Reads = colSums(df, na.rm = TRUE)
    )
    g2 <- ggplot(tot_reads, aes(x = Sample, y = Reads)) +
      geom_col(fill = "#2ca02c") +
      labs(x = "Sample", y = "Total Reads") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    ggplotly(g2)
  })
  
  # 9) Heatmap
  output$heatmap <- renderPlot({
    req(input$analysis == "Heatmap", abund_df())
    df <- abund_df()
    pheatmap::pheatmap(df, color = viridis::viridis(50), angle_col = 45)
  })
  
  # 10) PCoA (Bray-Curtis) – transpose abundance so samples are rows
  output$pcoa <- renderPlot({
    req(input$analysis == "PCoA (Bray-Curtis)", abund_df(), meta_df(), input$pcoa_group)
    df_ab <- abund_df()
    df_meta <- meta_df()
    grp <- input$pcoa_group
    df_t <- t(df_ab)  # now samples × taxa
    if (!all(rownames(df_t) %in% df_meta$SampleID)) {
      stop("After cleaning and transposing, rownames(df_t) do not match metadata SampleID.")
    }
    dist_mat <- vegan::vegdist(df_t, method = "bray")
    pc <- ape::pcoa(dist_mat)
    coords <- as.data.frame(pc$vectors[, 1:2]) %>%
      rownames_to_column("SampleID")
    coords <- left_join(coords, df_meta, by = "SampleID")
    
    # DEBUG: print sample IDs and grouping values
    message("=== In PCoA: coords$SampleID = ")
    message(paste(coords$SampleID, collapse = ", "))
    message("=== In PCoA: coords[[ grp ]] = ")
    message(paste(coords[[grp]], collapse = ", "))
    
    coords[[grp]] <- as.factor(coords[[grp]])
    cent <- coords %>%
      group_by_at(grp) %>%
      summarize(
        cx = mean(Axis.1, na.rm = TRUE),
        cy = mean(Axis.2, na.rm = TRUE),
        .groups = "drop"
      )
    segs <- left_join(coords, cent, by = grp)
    ggplot(coords, aes(x = Axis.1, y = Axis.2, color = .data[[grp]])) +
      geom_segment(
        data = segs,
        aes(x = Axis.1, y = Axis.2, xend = cx, yend = cy),
        alpha = 0.5
      ) +
      geom_point(size = 3) +
      stat_ellipse(aes(color = .data[[grp]]), type = "norm", level = 0.95) +
      labs(
        x = "PCoA 1",
        y = "PCoA 2",
        color = grp,
        title = "PCoA (Bray-Curtis)"
      ) +
      theme_minimal()
  })
}

# ---------------------------------------------------------------------------
# Launch the app
# ---------------------------------------------------------------------------
shinyApp(ui, server)
