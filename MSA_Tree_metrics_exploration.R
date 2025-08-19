# Load packages
library(readr)
library(ggplot2)
library(dplyr)
library(purrr)
library(shiny)
library(bslib)
library(corrplot)

# Define Tukey's Fences outlier detection function
detect_outliers_tukey <- function(x, k = 1.5) {
  # Check if input is numeric
  if (!is.numeric(x)) {
    return(rep(FALSE, length(x)))
  }
  
  # Calculate quartiles and IQR
  Q1 <- quantile(x, 0.25, na.rm = TRUE)
  Q3 <- quantile(x, 0.75, na.rm = TRUE)
  IQR <- Q3 - Q1
  
  # Calculate outlier bounds
  lower_bound <- Q1 - k * IQR
  upper_bound <- Q3 + k * IQR
  
  # Return logical vector indicating outliers
  is_outlier <- (x < lower_bound | x > upper_bound) & !is.na(x)
  
  return(is_outlier)
}

# Function to format numbers by removing trailing zeros
format_number_clean <- function(x, digits = 6) {
  formatted <- sprintf(paste0("%.", digits, "f"), x)
  # Remove trailing zeros and decimal point if not needed
  formatted <- sub("0+$", "", formatted)
  formatted <- sub("\\.$", "", formatted)
  return(formatted)
}

# Define UI ----
ui <- page_sidebar(
  title = "MSA and Tree metrics exploration and MSA simulation of empirical phylogenomic datasets",
  sidebar = sidebar(
    width = 300,
    helpText("Parameters/properties are divided into three groups: sequence-, substitution model-, gene tree-related ones"),
    
    # Choice between default files and upload
    checkboxGroupInput(
      "use_default_files",
      label = "Use default data files (from app directory):",
      choices = list(
        "DNA data (metrics.GTR.csv)" = "dna",
        "AA data (metrics.LG.csv)" = "aa"
      ),
      selected = c("dna", "aa")
    ),
    
    # File input for custom data (only one needed now)
    fileInput(
      "upload_file",
      label = "Or upload your own data (.csv):",
      accept = c(".csv"),
      multiple = FALSE, 
      buttonLabel = "Browse...",
      placeholder = "No file selected"
    ),
    
    selectInput(
      "var_type",
      label = "Sequence type of MSA alignments:",
      choices = list(
        "nucleotides (DNA)" = "dna",
        "amino acids (AA)" = "aa"
      ),
      selected = "dna"
    ),
    
    # Use a single selectInput for the actual variable to plot
    selectInput(
      "selected_variable",
      label = "Select Variable to Plot:",
      choices = NULL
    ),
    
    # New checkbox for correlation analysis
    checkboxInput(
      "enable_correlation",
      label = "Enable Correlation Analysis",
      value = FALSE
    ),
    
    # New checkbox for MSA simulation
    checkboxInput(
      "enable_msa_simulation",
      label = "Enable MSA Simulation",
      value = FALSE
    )
  ),
  
  # Main content area with navigation
  navset_card_tab(
    nav_panel(
      "Distribution Analysis",
      layout_sidebar(
        fillable = TRUE,
        sidebar = sidebar(
          width = 350,
          title = "Details & Settings",
          
          # Basic Statistics
          card(
            card_header("Basic Statistics"),
            verbatimTextOutput("basic_stats")
          ),
          
          # Histogram Settings
          card(
            card_header("Histogram Settings"),
            numericInput(
              "bins",
              label = "Number of bins:",
              value = 50,
              min = 1,
              step = 1
            ),
            numericInput(
              "xmin",
              label = "X-axis minimum:",
              value = NA
            ),
            numericInput(
              "xmax",
              label = "X-axis maximum:",
              value = NA
            )
          ),
          
          # Sample Exclusion Filter
          card(
            card_header("Sample Exclusion Filter"),
            
            # Tukey's Fences option
            checkboxInput(
              "use_tukey_filter",
              label = "Use Tukey's Fences to filter outliers",
              value = FALSE
            ),
            
            conditionalPanel(
              condition = "input.use_tukey_filter",
              numericInput(
                "tukey_k",
                label = "Tukey's k value:",
                value = 1.5,
                min = 0.1,
                max = 5,
                step = 0.1
              ),
              selectInput(
                "tukey_variable",
                label = "Apply Tukey's filter to which variable:",
                choices = NULL,
                selected = NULL
              ),
              radioButtons(
                "tukey_method",
                label = "Tukey's filtering method:",
                choices = list(
                  "Remove entire rows with outliers" = "remove_rows",
                  "Replace outliers with NA only" = "replace_na"
                ),
                selected = "remove_rows"
              )
            ),
            
            helpText("Exclude samples where the selected variable's value falls WITHIN the specified range."),
            selectInput(
              "filter_variable",
              label = "Select Variable to Exclude By:",
              choices = NULL,
              selected = NULL
            ),
            
            numericInput(
              "exclude_min", 
              label = "Exclude Value From (inclusive):",
              value = NA
            ),
            numericInput(
              "exclude_max", 
              label = "Exclude Value To (inclusive):",
              value = NA
            )
          ),
          
          # Plot Size Settings
          card(
            card_header("Plot Size Settings"),
            numericInput(
              "plot_width",
              label = "Plot width (inches):",
              value = 10,
              min = 1,
              max = 20,
              step = 0.5
            ),
            numericInput(
              "plot_height",
              label = "Plot height (inches):",
              value = 10,
              min = 1,
              max = 20,
              step = 0.5
            )
          ),
          
          # Download options
          card(
            card_header("Download Options"),
            downloadButton(
              "download_filtered_data",
              label = "Download Filtered Data (CSV)",
              class = "btn-primary"
            ),
            br(), br(),
            downloadButton(
              "download_plot_pdf",
              label = "Download Plot (PDF)",
              class = "btn-secondary"
            )
          )
        ),
        
        # Main content for the plot
        card(
          plotOutput("map_density")
        )
      )
    ),
    
    # Tab for correlation analysis
    nav_panel(
      "Correlation Analysis",
      conditionalPanel(
        condition = "input.enable_correlation",
        layout_sidebar(
          fillable = TRUE,
          sidebar = sidebar(
            width = 350,
            title = "Correlation Settings",
            
            # Variable Selection
            card(
              card_header("Variable Selection"),
              helpText("Select at least 2 variables for correlation analysis:"),
              checkboxGroupInput(
                "correlation_variables",
                label = "Select Variables:",
                choices = NULL
              )
            ),
            
            # Correlation Plot Settings
            card(
              card_header("Heatmap Settings"),
              radioButtons(
                "heatmap_type",
                label = "Display type:",
                choices = list(
                  "Full matrix" = "full",
                  "Lower triangle" = "lower",
                  "Upper triangle" = "upper"
                ),
                selected = "full"
              ),
              numericInput(
                "cluster_rectangles",
                label = "Number of cluster rectangles (optional):",
                value = NA,
                min = 1,
                step = 1
              ),
              helpText("Leave empty to show no cluster rectangles")
            ),
            
            # Heatmap Size Settings
            card(
              card_header("Heatmap Size Settings"),
              numericInput(
                "heatmap_width",
                label = "Heatmap width (inches):",
                value = 10,
                min = 1,
                max = 20,
                step = 0.5
              ),
              numericInput(
                "heatmap_height",
                label = "Heatmap height (inches):",
                value = 10,
                min = 1,
                max = 20,
                step = 0.5
              )
            ),
            
            # Download options for correlation
            card(
              card_header("Download Options"),
              downloadButton(
                "download_correlation_matrix",
                label = "Download Correlation Matrix (CSV)",
                class = "btn-primary"
              ),
              br(), br(),
              downloadButton(
                "download_correlation_plot",
                label = "Download Heatmap (PDF)",
                class = "btn-secondary"
              )
            )
          ),
          
          # Main content for correlation plot
          card(
            card_header("Correlation Heatmap"),
            plotOutput("correlation_plot", height = "600px")
          )
        )
      ),
      conditionalPanel(
        condition = "!input.enable_correlation",
        div(
          class = "text-center mt-5",
          h4("Correlation Analysis Disabled"),
          p("Please enable correlation analysis in the left sidebar to access this feature.")
        )
      )
    ),
    
    # Tab for MSA simulation
    nav_panel(
      "MSA Simulation",
      conditionalPanel(
        condition = "input.enable_msa_simulation",
        # add page description
        div(
          class = "alert alert-info mb-4",
          style = "margin: 5px;",
          p("NOTES: It is recommended to filter the data (e.g., using the csvtk tool) based on the exploratory analysis from Distribution Analysis and Correlation Analysis, re-upload the filtered csv file, and then proceed with the simulation. The variable names that are highly correlated with the simulation in the uploaded file should be consistent with the variable names demonstrated on the webpage, such as", 
            tags$code("GammaShape_alpha"), ", ", 
            tags$code("freqA"), ", ", 
            tags$code("rateAC"), ", ", 
            tags$code("num_sites"), "etc.", 
            style = "margin-bottom: 0;")
        ),
        layout_sidebar(
          fillable = TRUE,
          sidebar = sidebar(
            width = 350,
            title = "Simulation Settings",
            
            # Basic simulation settings
            card(
              card_header("Basic Settings"),
              textInput(
                "sim_prefix",
                label = "Sequence file prefix:",
                value = "sim",
                placeholder = "sim"
              ),
              numericInput(
                "num_alignments",
                label = "Number of MSA simulations:",
                value = 1000,
                min = 1,
                step = 1
              ),
              numericInput(
                "num_threads",
                label = "Number of threads (-nt):",
                value = 1,
                min = 1,
                step = 1
              )
            ),
            
            # Simulation strategy
            card(
              card_header("Simulation Strategy"),
              radioButtons(
                "sim_strategy",
                label = "Parameter sampling strategy:",
                choices = list(
                  "Complete empirical match (all parameters from same sample)" = "complete_match",
                  "Mixed empirical sampling (parameters from different samples)" = "mixed_empirical",
                  "Probability density function (PDF) estimation for some parameters" = "pdf_based"
                ),
                selected = "complete_match"
              )
            ),
            
            # Model components - modified for correct AA sequences
            card(
              card_header("Select model components:"),
              checkboxGroupInput(
                "model_components",
                label = NULL,
                choices = list(
                  "Base model (GTR for DNA/LG for AA)" = "base_model",
                  "Frequencies (+F)" = "frequencies",
                  "Gamma shape (+G4)" = "gamma",
                  "Invariant sites (+I)" = "invariant"
                ),
                selected = c("base_model", "frequencies", "gamma", "invariant")
              )
            ),
            
            # GTR parameter source for DNA (only show for DNA sequences)
            conditionalPanel(
              condition = "input.var_type == 'dna' && input.model_components.includes('base_model') && input.sim_strategy != 'complete_match'",
              card(
                card_header("GTR Rate Parameters Source"),
                radioButtons(
                  "gtr_source",
                  label = "GTR rate parameters source:",
                  choices = list(
                    "From same sample (maintain covariance)" = "same_sample",
                    "Mixed empirical sampling" = "mixed_empirical",
                    "Probability density function (PDF) estimation" = "pdf_based"
                  ),
                  selected = "same_sample"
                )
              )
            ),
            
            # Length settings - conditional display
            conditionalPanel(
              condition = "input.sim_strategy != 'complete_match'",
              card(
                card_header("Sequence Length Settings"),
                radioButtons(
                  "length_source",
                  label = "Sequence length source:",
                  choices = list(
                    "From empirical data (num_sites)" = "empirical",
                    "Fixed length" = "fixed"
                  ),
                  selected = "empirical"
                ),
                conditionalPanel(
                  condition = "input.length_source == 'fixed'",
                  numericInput(
                    "fixed_length",
                    label = "Fixed sequence length:",
                    value = 500,
                    min = 1,
                    step = 1
                  )
                )
              )
            ),
            
            # Tree settings (only for non-complete_match strategies)
            conditionalPanel(
              condition = "input.sim_strategy != 'complete_match'",
              card(
                card_header("Tree Settings"),
                radioButtons(
                  "tree_source",
                  label = "Tree source:",
                  choices = list(
                    "Random from empirical data (if available)" = "empirical",
                    "Use custom tree file" = "custom"
                  ),
                  selected = "empirical"
                ),
                conditionalPanel(
                  condition = "input.tree_source == 'custom'",
                  helpText("Tree parameter will be set to 'NEWICK_TREE' - replace with your tree file name")
                )
              )
            ),
            
            # PDF settings (only for pdf_based strategy)
            conditionalPanel(
              condition = "input.sim_strategy == 'pdf_based'",
              card(
                card_header("Probability density function (PDF) estimation Settings"),
                checkboxGroupInput(
                  "pdf_parameters",
                  label = "Use PDF estimation for these parameters:",
                  choices = list(
                    "Gamma shape (GammaShape_alpha)" = "gamma_pdf",
                    "Invariant sites (proportion_invariant)" = "invariant_pdf",
                    "Sequence length (num_sites)" = "length_pdf"
                  ),
                  selected = c()
                ),
                
                # GTR rate parameters for DNA when both PDF strategy and GTR PDF source are selected
                conditionalPanel(
                  condition = "input.var_type == 'dna' && input.model_components.includes('base_model') && input.gtr_source == 'pdf_based'",
                  hr(),
                  h6("GTR Rate Parameters PDF estimation:"),
                  checkboxGroupInput(
                    "gtr_pdf_parameters",
                    label = "Use PDF estimation for GTR rate parameters:",
                    choices = list(
                      "rateAC" = "rateAC_pdf",
                      "rateAG" = "rateAG_pdf", 
                      "rateAT" = "rateAT_pdf",
                      "rateCG" = "rateCG_pdf",
                      "rateCT" = "rateCT_pdf"
                    ),
                    selected = c()
                  )
                ),
                
                helpText("Other parameters will use empirical sampling. Note: If 'Sequence length' is selected here, it will override the 'Sequence Length Settings' above.")
              )
            ),
            
            # PDF estimation visualization options
            conditionalPanel(
              condition = "input.sim_strategy == 'pdf_based'",
              card(
                card_header("PDF Estimation Visualization"),
                helpText("Check variables to visualize PDF fit:"),
                uiOutput("pdf_visualization_controls")
              )
            ),
            
            # Download options
            card(
              card_header("Download Options"),
              downloadButton(
                "download_simulation_commands",
                label = "Download IQ-TREE Commands",
                class = "btn-primary"
              ),
              
              # PDF visualization download options
              conditionalPanel(
                condition = "input.sim_strategy == 'pdf_based' && output.has_pdf_vars",
                br(),
                h6("PDF Estimation Visualization Download:"),
                numericInput(
                  "pdf_plot_width",
                  label = "Plot width (inches):",
                  value = 12,
                  min = 1,
                  max = 20,
                  step = 0.5
                ),
                numericInput(
                  "pdf_plot_height",
                  label = "Plot height (inches):",
                  value = 8,
                  min = 1,
                  max = 20,
                  step = 0.5
                ),
                downloadButton(
                  "download_pdf_comparison_plot",
                  label = "Download PDF Comparison Plot (PDF)",
                  class = "btn-secondary"
                )
              )
            )
          ),
          
          # Main content for simulation preview
          div(
            card(
              card_header("Command Preview"),
              verbatimTextOutput("command_preview"),
              br(),
              h5("Statistics:"),
              verbatimTextOutput("simulation_stats")
            ),
            
            # PDF Comparison Plot (moved here from separate tab)
            conditionalPanel(
              condition = "input.sim_strategy == 'pdf_based'",
              card(
                card_header("PDF Estimation vs Empirical Data Comparison"),
                plotOutput("pdf_comparison_plot", height = "600px")
              )
            )
          )
        )
      ),
      conditionalPanel(
        condition = "!input.enable_msa_simulation",
        div(
          class = "text-center mt-5",
          h4("MSA Simulation Disabled"),
          p("Please enable MSA simulation in the left sidebar to access this feature.")
        )
      )
    )
  )
)

# Define server logic ----
server <- function(input, output, session) {
  
  # Set maximum request size
  options(shiny.maxRequestSize = 1024 * 1024 * 1024) # 1 GB
  
  # Reactive expression to read data
  current_data <- reactive({
    # Check if custom file is uploaded
    if (!is.null(input$upload_file)) {
      df <- tryCatch({
        read_csv(input$upload_file$datapath)
      }, error = function(e) {
        showNotification(paste("Error reading uploaded file:", e$message), type = "error")
        return(NULL)
      })
      return(df)
    }
    
    # Otherwise, use default files
    dfs <- list()
    
    if ("dna" %in% input$use_default_files) {
      if (file.exists("metrics.GTR.csv")) {
        dfs[["dna"]] <- tryCatch({
          read_csv("metrics.GTR.csv")
        }, error = function(e) {
          showNotification("Error reading metrics.GTR.csv", type = "warning")
          NULL
        })
      } else {
        showNotification("metrics.GTR.csv not found in app directory", type = "warning")
      }
    }
    
    if ("aa" %in% input$use_default_files) {
      if (file.exists("metrics.LG.csv")) {
        dfs[["aa"]] <- tryCatch({
          read_csv("metrics.LG.csv")
        }, error = function(e) {
          showNotification("Error reading metrics.LG.csv", type = "warning")
          NULL
        })
      } else {
        showNotification("metrics.LG.csv not found in app directory", type = "warning")
      }
    }
    
    # Filter and return based on selected type
    if (input$var_type %in% names(dfs) && !is.null(dfs[[input$var_type]])) {
      return(dfs[[input$var_type]])
    } else if (length(dfs) > 0) {
      # Return the first available dataset
      return(dfs[[1]])
    } else {
      return(NULL)
    }
  })
  
  # Update variable choices based on current data
  observe({
    current_df <- current_data()
    req(current_df)
    
    # Only consider numeric columns for analysis
    numeric_cols <- names(current_df)[sapply(current_df, is.numeric)]
    
    # Define predefined groups (only include numeric variables)
    msa_choices <- c(
      "num_taxa", "num_sites", "num_patterns", "num_parsimony_sites",
      "num_sites/num_taxa", "num_patterns/num_taxa", "num_parsimony_sites/num_taxa",
      "num_patterns/num_sites", "num_parsimony_sites/num_sites",
      "proportion_gaps", "proportion_invariant", "entropy", "bollback",
      "pattern_entropy", "rcfv", "average_pairwise_identity"
    )
    
    model_choices_dna <- c(
      "GammaShape_alpha",
      "freqA", "freqC", "freqG", "freqT",
      "rateAC", "rateAG", "rateAT", "rateCG", "rateCT", "rateGT"
    )
    
    model_choices_aa <- c(
      "GammaShape_alpha",
      "freqA", "freqR", "freqN", "freqD", "rateC",
      "rateQ", "rateE", "rateG", "rateH", "rateI",
      "freqL", "freqK", "freqM", "freqF", "freqP",
      "freqS", "freqT", "freqW", "freqY", "freqV"
    )
    
    gene_tree_choices <- c(
      "average_BS", "sd_BS", "total_tree_length", "average_internal_branch_length",
      "sd_internal_branch_length", "average_terminal_branch_length",
      "sd_terminal_branch_length", "tree_diameter", "average_patristic_distance",
      "sd_patristic_distance", "evo_rate", "treeness", "dvmc", "RF_distance"
    )
    
    # Filter choices based on available numeric columns
    filter_choices <- function(choices_list, available_cols) {
      choices_list[choices_list %in% available_cols]
    }
    
    filtered_msa_choices <- filter_choices(msa_choices, numeric_cols)
    filtered_gene_tree_choices <- filter_choices(gene_tree_choices, numeric_cols)
    
    model_choices_current <- if (input$var_type == "dna") {
      filter_choices(model_choices_dna, numeric_cols)
    } else {
      filter_choices(model_choices_aa, numeric_cols)
    }
    
    # Find other numeric variables (not in predefined groups)
    all_predefined <- c(msa_choices, model_choices_dna, model_choices_aa, gene_tree_choices)
    other_numeric_choices <- setdiff(numeric_cols, all_predefined)
    
    # Combine all choices
    all_choices <- list()
    if (length(filtered_msa_choices) > 0) {
      all_choices[["Sequence-related properties"]] <- as.list(filtered_msa_choices)
    }
    if (length(model_choices_current) > 0) {
      all_choices[["Substitution model-related parameters"]] <- as.list(model_choices_current)
    }
    if (length(filtered_gene_tree_choices) > 0) {
      all_choices[["Gene tree-related metrics"]] <- as.list(filtered_gene_tree_choices)
    }
    if (length(other_numeric_choices) > 0) {
      all_choices[["Other variables"]] <- as.list(other_numeric_choices)
    }
    
    # Update selected_variable
    current_selection <- input$selected_variable
    new_selected <- if (!is.null(current_selection) && current_selection %in% numeric_cols) {
      current_selection
    } else if (length(filtered_msa_choices) > 0) {
      filtered_msa_choices[[1]]
    } else if (length(numeric_cols) > 0) {
      numeric_cols[[1]]
    } else {
      NULL
    }
    
    updateSelectInput(session, "selected_variable",
                      choices = all_choices,
                      selected = new_selected)
    
    # Update correlation variables choices
    correlation_choices_list <- as.list(numeric_cols)
    names(correlation_choices_list) <- numeric_cols
    
    updateCheckboxGroupInput(session, "correlation_variables",
                             choices = correlation_choices_list,
                             selected = NULL)
    
    # Update filter variable choices (only numeric columns)
    filter_choices_list <- as.list(numeric_cols)
    names(filter_choices_list) <- numeric_cols
    
    current_filter_selection <- input$filter_variable
    new_filter_selected <- if (!is.null(current_filter_selection) && current_filter_selection %in% numeric_cols) {
      current_filter_selection
    } else if ("num_sites" %in% numeric_cols) {
      "num_sites"
    } else if (length(numeric_cols) > 0) {
      numeric_cols[[1]]
    } else {
      NULL
    }
    
    updateSelectInput(session, "filter_variable",
                      choices = filter_choices_list,
                      selected = new_filter_selected)
    
    # Update Tukey variable choices
    current_tukey_selection <- input$tukey_variable
    new_tukey_selected <- if (!is.null(current_tukey_selection) && current_tukey_selection %in% numeric_cols) {
      current_tukey_selection
    } else if ("num_sites" %in% numeric_cols) {
      "num_sites"
    } else if (length(numeric_cols) > 0) {
      numeric_cols[[1]]
    } else {
      NULL
    }
    
    updateSelectInput(session, "tukey_variable",
                      choices = filter_choices_list,
                      selected = new_tukey_selected)
  })
  
  # Generate PDF visualization controls with individual X-axis controls
  output$pdf_visualization_controls <- renderUI({
    df <- filtered_data()
    req(df)
    
    # Get all potential PDF variables
    all_pdf_vars <- c()
    
    # Regular PDF parameters
    if ("gamma_pdf" %in% input$pdf_parameters && "GammaShape_alpha" %in% names(df)) {
      all_pdf_vars <- c(all_pdf_vars, "GammaShape_alpha")
    }
    if ("invariant_pdf" %in% input$pdf_parameters && "proportion_invariant" %in% names(df)) {
      all_pdf_vars <- c(all_pdf_vars, "proportion_invariant")
    }
    if ("length_pdf" %in% input$pdf_parameters && "num_sites" %in% names(df)) {
      all_pdf_vars <- c(all_pdf_vars, "num_sites")
    }
    
    # GTR parameters for DNA
    if (input$var_type == "dna" && !is.null(input$gtr_pdf_parameters)) {
      gtr_vars <- c("rateAC", "rateAG", "rateAT", "rateCG", "rateCT")
      for (gtr_var in gtr_vars) {
        param_name <- paste0(gtr_var, "_pdf")
        if (param_name %in% input$gtr_pdf_parameters && gtr_var %in% names(df)) {
          all_pdf_vars <- c(all_pdf_vars, gtr_var)
        }
      }
    }
    
    if (length(all_pdf_vars) == 0) {
      return(p("No PDF parameters selected for visualization."))
    }
    
    # Create controls for each PDF variable
    controls <- lapply(all_pdf_vars, function(var) {
      div(
        style = "border: 1px solid #ddd; padding: 10px; margin-bottom: 10px; border-radius: 5px;",
        h6(paste("Variable:", var), style = "margin-bottom: 8px;"),
        checkboxInput(
          inputId = paste0("show_pdf_", var),
          label = paste("Show", var, "in plot"),
          value = FALSE
        ),
        conditionalPanel(
          condition = paste0("input.show_pdf_", var),
          div(
            style = "margin-left: 20px;",
            h6("X-axis Range for this variable:", style = "font-size: 12px; margin-bottom: 5px;"),
            fluidRow(
              column(6,
                     numericInput(
                       paste0("pdf_xmin_", var),
                       label = "Min:",
                       value = NA,
                       step = if(var == "num_sites") 50 else 0.1
                     )
              ),
              column(6,
                     numericInput(
                       paste0("pdf_xmax_", var),
                       label = "Max:",
                       value = NA,
                       step = if(var == "num_sites") 50 else 0.1
                     )
              )
            ),
            helpText("Leave empty for automatic range", style = "font-size: 11px; margin-top: 5px;")
          )
        )
      )
    })
    
    return(div(controls))
  })
  
  # Check if there are PDF variables to show (for conditional panels)
  output$has_pdf_vars <- reactive({
    df <- filtered_data()
    if (is.null(df)) return(FALSE)
    
    if (input$sim_strategy != "pdf_based") return(FALSE)
    
    # Check if any PDF parameters are selected
    has_any <- FALSE
    if ("gamma_pdf" %in% input$pdf_parameters && "GammaShape_alpha" %in% names(df)) has_any <- TRUE
    if ("invariant_pdf" %in% input$pdf_parameters && "proportion_invariant" %in% names(df)) has_any <- TRUE
    if ("length_pdf" %in% input$pdf_parameters && "num_sites" %in% names(df)) has_any <- TRUE
    
    # Check GTR parameters
    if (input$var_type == "dna" && !is.null(input$gtr_pdf_parameters)) {
      gtr_vars <- c("rateAC", "rateAG", "rateAT", "rateCG", "rateCT")
      for (gtr_var in gtr_vars) {
        param_name <- paste0(gtr_var, "_pdf")
        if (param_name %in% input$gtr_pdf_parameters && gtr_var %in% names(df)) {
          has_any <- TRUE
          break
        }
      }
    }
    
    return(has_any)
  })
  outputOptions(output, "has_pdf_vars", suspendWhenHidden = FALSE)
  
  # Apply Tukey's fences filtering first
  tukey_filtered_data <- reactive({
    df <- current_data()
    req(df)
    
    if (input$use_tukey_filter && !is.null(input$tukey_variable) && input$tukey_variable %in% names(df)) {
      
      if (input$tukey_method == "remove_rows") {
        # Remove entire rows where the selected variable has outliers
        outliers <- detect_outliers_tukey(df[[input$tukey_variable]], k = input$tukey_k)
        df_result <- df[!outliers, ]
        
        # Show notification about how many rows were removed
        n_removed <- sum(outliers, na.rm = TRUE)
        if (n_removed > 0) {
          showNotification(
            paste("Tukey's filter removed", n_removed, "rows with outliers in", input$tukey_variable),
            type = "message", duration = 3
          )
        }
        
      } else {
        # Replace outliers with NA only in the selected variable
        df_result <- df
        outliers <- detect_outliers_tukey(df[[input$tukey_variable]], k = input$tukey_k)
        df_result[[input$tukey_variable]][outliers] <- NA
        
        # Show notification about how many values were replaced
        n_replaced <- sum(outliers, na.rm = TRUE)
        if (n_replaced > 0) {
          showNotification(
            paste("Tukey's filter replaced", n_replaced, "outlier values with NA in", input$tukey_variable),
            type = "message", duration = 3
          )
        }
      }
      
      return(df_result)
    } else {
      return(df)
    }
  })
  
  # Apply additional range-based filtering
  filtered_data <- reactive({
    df <- tukey_filtered_data()
    req(df)
    
    filter_col_name <- input$filter_variable
    exclude_min_val <- input$exclude_min
    exclude_max_val <- input$exclude_max
    
    if (!is.null(filter_col_name) && filter_col_name %in% names(df) && is.numeric(df[[filter_col_name]])) {
      col_data <- df[[filter_col_name]]
      is_in_exclusion_range <- rep(FALSE, nrow(df))
      non_na_indices <- which(!is.na(col_data))
      
      if (length(non_na_indices) > 0) {
        temp_col_data <- col_data[non_na_indices]
        
        if (!is.na(exclude_min_val) && !is.na(exclude_max_val)) {
          is_in_exclusion_range[non_na_indices] <- (temp_col_data >= exclude_min_val & temp_col_data <= exclude_max_val)
        } else if (!is.na(exclude_min_val) && is.na(exclude_max_val)) {
          is_in_exclusion_range[non_na_indices] <- (temp_col_data >= exclude_min_val)
        } else if (is.na(exclude_min_val) && !is.na(exclude_max_val)) {
          is_in_exclusion_range[non_na_indices] <- (temp_col_data <= exclude_max_val)
        }
      }
      
      df <- df[!is_in_exclusion_range, ]
    }
    
    return(df)
  })
  
  # Get selected column data
  selected_column_data <- reactive({
    req(input$selected_variable)
    df_filtered <- filtered_data()
    req(df_filtered)
    
    data_col <- df_filtered[[input$selected_variable]]
    if (is.null(data_col)) {
      showNotification(
        paste("Selected variable '", input$selected_variable, "' not found in the filtered data."),
        type = "warning", duration = 5
      )
      return(numeric(0))
    }
    data_col
  })
  
  # Calculate basic statistics
  basic_statistics <- reactive({
    plot_data <- selected_column_data()
    
    if (length(plot_data) == 0 || all(is.na(plot_data))) {
      return("No data available for statistics after filtering.")
    }
    
    plot_data_clean <- na.omit(plot_data)
    
    if (length(plot_data_clean) == 0) {
      return("No valid data points for statistics after removing NAs and filtering.")
    }
    
    mean_val <- mean(plot_data_clean)
    median_val <- median(plot_data_clean)
    min_val <- min(plot_data_clean)
    max_val <- max(plot_data_clean)
    q25_val <- quantile(plot_data_clean, 0.25, na.rm = TRUE)
    q75_val <- quantile(plot_data_clean, 0.75, na.rm = TRUE)
    sd_val <- sd(plot_data_clean, na.rm = TRUE)
    
    paste0(
      "Mean:       ", format(mean_val, digits = 4), "\n",
      "Median:     ", format(median_val, digits = 4), "\n",
      "Min:        ", format(min_val, digits = 4), "\n",
      "Max:        ", format(max_val, digits = 4), "\n",
      "25th Pctl:  ", format(q25_val, digits = 4), "\n",
      "75th Pctl:  ", format(q75_val, digits = 4), "\n",
      "Std Dev:    ", format(sd_val, digits = 4), "\n",
      "N (without NA):  ", length(plot_data_clean), "\n",
      "N (total):  ", length(plot_data)
    )
  })
  
  output$basic_stats <- renderText({
    basic_statistics()
  })
  
  # Create the plot
  create_plot <- function() {
    df_to_plot <- filtered_data()
    plot_data_column <- selected_column_data()
    
    if (length(plot_data_column) == 0 || all(is.na(plot_data_column))) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "No data to display after filtering or selected variable is empty/all NA.") + 
               theme_void())
    }
    
    x_min <- input$xmin
    if (is.na(x_min)) {
      x_min <- min(plot_data_column, na.rm = TRUE)
    }
    
    x_max <- input$xmax
    if (is.na(x_max)) {
      x_max <- max(plot_data_column, na.rm = TRUE)
    }
    
    if (is.null(df_to_plot) || !(input$selected_variable %in% names(df_to_plot))) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "Data not loaded or selected variable is missing after filtering.") + 
               theme_void())
    }
    
    # Create exclusion text for plot title
    exclusion_text <- ""
    if (input$use_tukey_filter && !is.null(input$tukey_variable)) {
      if (input$tukey_method == "remove_rows") {
        exclusion_text <- paste0(exclusion_text, "\n(Tukey's fences applied to ", input$tukey_variable, " with k=", input$tukey_k, ", rows removed)")
      } else {
        exclusion_text <- paste0(exclusion_text, "\n(Tukey's fences applied to ", input$tukey_variable, " with k=", input$tukey_k, ", outliers → NA)")
      }
    }
    if (!is.na(input$exclude_min) && !is.na(input$exclude_max)) {
      exclusion_text <- paste0(exclusion_text, "\n(Excluding samples where ", input$filter_variable, " is between ", input$exclude_min, " and ", input$exclude_max, ")")
    } else if (!is.na(input$exclude_min)) {
      exclusion_text <- paste0(exclusion_text, "\n(Excluding samples where ", input$filter_variable, " is >= ", input$exclude_min, ")")
    } else if (!is.na(input$exclude_max)) {
      exclusion_text <- paste0(exclusion_text, "\n(Excluding samples where ", input$filter_variable, " is <= ", input$exclude_max, ")")
    }
    
    p <- ggplot(df_to_plot, aes(x = .data[[input$selected_variable]])) +
      geom_histogram(aes(y = after_stat(density)),
                     fill = "#D4EDDA",
                     color = "#B0B0B0",
                     alpha = 0.7,
                     bins = input$bins,
                     na.rm = TRUE) +
      geom_density(color = "#FF8C00", na.rm = TRUE) +
      theme_bw() +
      labs(x = input$selected_variable, 
           title = paste("Distribution of", input$selected_variable, exclusion_text))
    
    if (is.numeric(x_min) && is.numeric(x_max) && x_min < x_max) {
      p <- p + xlim(x_min, x_max)
    }
    
    return(p)
  }
  
  output$map_density <- renderPlot({
    create_plot()
  })
  
  # Create PDF comparison plot with individual X-axis controls
  create_pdf_comparison_plot <- function() {
    df <- filtered_data()
    req(df)
    
    # Get variables to plot
    vars_to_plot <- character(0)
    
    # Check regular PDF parameters
    if ("gamma_pdf" %in% input$pdf_parameters && "GammaShape_alpha" %in% names(df) && 
        !is.null(input$show_pdf_GammaShape_alpha) && input$show_pdf_GammaShape_alpha) {
      vars_to_plot <- c(vars_to_plot, "GammaShape_alpha")
    }
    if ("invariant_pdf" %in% input$pdf_parameters && "proportion_invariant" %in% names(df) && 
        !is.null(input$show_pdf_proportion_invariant) && input$show_pdf_proportion_invariant) {
      vars_to_plot <- c(vars_to_plot, "proportion_invariant")
    }
    if ("length_pdf" %in% input$pdf_parameters && "num_sites" %in% names(df) && 
        !is.null(input$show_pdf_num_sites) && input$show_pdf_num_sites) {
      vars_to_plot <- c(vars_to_plot, "num_sites")
    }
    
    # Check GTR parameters
    if (input$var_type == "dna" && !is.null(input$gtr_pdf_parameters)) {
      gtr_vars <- c("rateAC", "rateAG", "rateAT", "rateCG", "rateCT")
      for (gtr_var in gtr_vars) {
        param_name <- paste0(gtr_var, "_pdf")
        show_var_name <- paste0("show_pdf_", gtr_var)
        if (param_name %in% input$gtr_pdf_parameters && gtr_var %in% names(df) && 
            !is.null(input[[show_var_name]]) && input[[show_var_name]]) {
          vars_to_plot <- c(vars_to_plot, gtr_var)
        }
      }
    }
    
    if (length(vars_to_plot) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "No variables selected for PDF visualization") + 
               theme_void())
    }
    
    # Create comparison plots
    plot_list <- list()
    
    for (var in vars_to_plot) {
      data_col <- na.omit(df[[var]])
      if (length(data_col) < 2) next
      
      # Generate PDF samples
      samples <- generate_pdf_samples(var, 1000, df)
      if (is.null(samples)) next
      
      # Create data frame for plotting
      df1 <- data.frame(value = data_col, group = "empirical")
      df2 <- data.frame(value = samples, group = "simulated")
      combined_df <- rbind(df1, df2)
      
      # Create plot with individual x-axis range controls
      p <- ggplot(combined_df, aes(x = value, color = group)) +
        geom_density(alpha = 0.5, position = "identity") +
        scale_color_manual(values = c("empirical" = "blue", "simulated" = "red")) +
        labs(x = var, y = "Density", title = paste("PDF vs Empirical:", var)) +
        theme_minimal()
      
      # Apply individual x-axis range for this variable
      var_xmin_input <- paste0("pdf_xmin_", var)
      var_xmax_input <- paste0("pdf_xmax_", var)
      var_xmin <- input[[var_xmin_input]]
      var_xmax <- input[[var_xmax_input]]
      
      if (!is.null(var_xmin) && !is.null(var_xmax) && 
          !is.na(var_xmin) && !is.na(var_xmax) && 
          var_xmin < var_xmax) {
        p <- p + xlim(var_xmin, var_xmax)
      } else {
        # Adjust xlim based on variable type (automatic range)
        if (var == "num_sites") {
          # For sequence length, use appropriate range
          x_range <- range(combined_df$value, na.rm = TRUE)
          p <- p + xlim(max(0, x_range[1] - 100), x_range[2] + 100)
        } else {
          # For other parameters, use smaller range
          p <- p + xlim(0, 3)
        }
      }
      
      plot_list[[var]] <- p
    }
    
    if (length(plot_list) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "No valid data for PDF visualization") + 
               theme_void())
    }
    
    # Arrange plots in a grid
    if (length(plot_list) == 1) {
      return(plot_list[[1]])
    } else {
      # Use gridExtra for better plot arrangement
      if (requireNamespace("gridExtra", quietly = TRUE)) {
        return(gridExtra::arrangeGrob(grobs = plot_list, ncol = 2))
      } else {
        # Fallback to first plot if gridExtra is not available
        return(plot_list[[1]])
      }
    }
  }
  
  # PDF comparison plot output
  output$pdf_comparison_plot <- renderPlot({
    plot_obj <- create_pdf_comparison_plot()
    
    if (inherits(plot_obj, "grob")) {
      # If it's a grob (from gridExtra), use grid.draw
      grid::grid.draw(plot_obj)
    } else {
      # If it's a ggplot, print it normally
      print(plot_obj)
    }
  })
  
  # Correlation analysis reactive expressions
  correlation_data <- reactive({
    req(input$enable_correlation, input$correlation_variables)
    
    if (length(input$correlation_variables) < 2) {
      return(NULL)
    }
    
    df <- filtered_data()
    req(df)
    
    # Select only the chosen variables
    correlation_vars <- input$correlation_variables
    available_vars <- correlation_vars[correlation_vars %in% names(df)]
    
    if (length(available_vars) < 2) {
      return(NULL)
    }
    
    data_for_correlation <- df[, available_vars, drop = FALSE]
    
    # Remove rows where all selected variables are NA
    complete_rows <- complete.cases(data_for_correlation)
    data_for_correlation <- data_for_correlation[complete_rows, , drop = FALSE]
    
    if (nrow(data_for_correlation) == 0) {
      return(NULL)
    }
    
    return(data_for_correlation)
  })
  
  correlation_matrix <- reactive({
    data <- correlation_data()
    req(data)
    
    # Calculate correlation matrix
    cor_matrix <- cor(data, method = "spearman", use = "pairwise.complete.obs")
    return(cor_matrix)
  })
  
  # Create correlation plot
  create_correlation_plot <- function() {
    cor_matrix <- correlation_matrix()
    req(cor_matrix)
    
    # Prepare corrplot parameters
    addrect_value <- if (is.na(input$cluster_rectangles)) NULL else input$cluster_rectangles
    
    corrplot(cor_matrix,
             method = "circle",
             type = input$heatmap_type,
             order = "hclust",
             hclust.method = "ward.D2",
             addrect = addrect_value,
             tl.srt = 45)
  }
  
  output$correlation_plot <- renderPlot({
    if (!input$enable_correlation) {
      return(NULL)
    }
    
    if (is.null(input$correlation_variables) || length(input$correlation_variables) < 2) {
      plot.new()
      text(0.5, 0.5, "Please select at least 2 variables for correlation analysis", cex = 1.2)
      return()
    }
    
    cor_data <- correlation_data()
    if (is.null(cor_data)) {
      plot.new()
      text(0.5, 0.5, "No valid data available for correlation analysis", cex = 1.2)
      return()
    }
    
    create_correlation_plot()
  })
  
  # MSA Simulation Functions
  
  # Generate PDF samples with noise
  generate_pdf_samples <- function(variable_name, sample_size, filtered_df) {
    if (!(variable_name %in% names(filtered_df))) {
      return(NULL)
    }
    
    data_col <- filtered_df[[variable_name]]
    data_clean <- na.omit(data_col)
    
    if (length(data_clean) < 2) {
      return(NULL)
    }
    
    # Create histogram for PDF estimation
    hist_result <- hist(data_clean, breaks = "FD", plot = FALSE)
    
    # Create PDF estimates using linear interpolation
    pdf_estimates <- approxfun(hist_result$mids, hist_result$density, method = "linear", rule = 2)
    
    # Sample from the estimated PDFs
    samples <- sample(hist_result$mids, size = sample_size, prob = hist_result$density, replace = TRUE)
    
    # Add noise based on variable type
    if (variable_name == "num_sites") {
      # For length parameters, add noise in range (-5, 5)
      noise <- runif(sample_size, -3, 3)
      samples <- samples + noise
      samples <- round(samples, 0)  # Round to integers
      samples[samples < 1] <- 1     # Ensure positive integers
    } else {
      # For other parameters, add small noise in range (-0.0001, 0.0001)
      noise <- runif(sample_size, -0.0001, 0.0001)
      samples <- samples + noise
    }
    
    return(samples)
  }
  
  # Generate simulation commands
  generate_simulation_commands <- reactive({
    req(input$enable_msa_simulation)
    
    df <- filtered_data()
    req(df)
    
    if (nrow(df) == 0) {
      return(list(commands = character(0), stats = "No data available after filtering."))
    }
    
    num_sims <- input$num_alignments
    seq_type <- toupper(input$var_type)  # DNA or AA
    
    # Required variables for DNA
    dna_freq_vars <- c("freqA", "freqC", "freqG", "freqT")
    dna_rate_vars <- c("rateAC", "rateAG", "rateAT", "rateCG", "rateCT")
    
    # Required variables for AA (corrected variable names to match actual column names)
    aa_freq_vars <- c("freqA", "freqR", "freqN", "freqD", "freqC",
                      "freqQ", "freqE", "freqG", "freqH", "freqI",
                      "freqL", "freqK", "freqM", "freqF", "freqP",
                      "freqS", "freqT", "freqW", "freqY", "freqV")
    
    commands <- character(num_sims)
    
    # Check which variables are available
    has_gamma <- "GammaShape_alpha" %in% names(df)
    has_invariant <- "proportion_invariant" %in% names(df)
    has_num_sites <- "num_sites" %in% names(df)
    has_loci <- "loci" %in% names(df)
    
    if (seq_type == "DNA") {
      has_frequencies <- all(dna_freq_vars %in% names(df))
      has_rates <- all(dna_rate_vars %in% names(df))
      freq_vars <- dna_freq_vars
      rate_vars <- dna_rate_vars
    } else {
      has_frequencies <- all(aa_freq_vars %in% names(df))
      has_rates <- FALSE  # LG model doesn't need rate parameters
      freq_vars <- aa_freq_vars
      rate_vars <- character(0)
    }
    
    # Debug information for AA sequences
    if (seq_type == "AA") {
      missing_freq_vars <- aa_freq_vars[!aa_freq_vars %in% names(df)]
      if (length(missing_freq_vars) > 0) {
        showNotification(
          paste("Missing AA frequency variables:", paste(missing_freq_vars, collapse = ", ")),
          type = "warning", duration = 5
        )
      }
    }
    
    # Generate PDF samples if needed
    pdf_samples <- list()
    if (input$sim_strategy == "pdf_based") {
      if ("gamma_pdf" %in% input$pdf_parameters && has_gamma) {
        pdf_samples$gamma <- generate_pdf_samples("GammaShape_alpha", num_sims, df)
      }
      if ("invariant_pdf" %in% input$pdf_parameters && has_invariant) {
        pdf_samples$invariant <- generate_pdf_samples("proportion_invariant", num_sims, df)
      }
      if ("length_pdf" %in% input$pdf_parameters && has_num_sites) {
        pdf_samples$length <- generate_pdf_samples("num_sites", num_sims, df)
      }
    }
    
    # Generate PDF samples for GTR rates if needed
    gtr_pdf_samples <- NULL
    if (seq_type == "DNA" && "base_model" %in% input$model_components && 
        input$sim_strategy == "pdf_based" && 
        !is.null(input$gtr_source) && input$gtr_source == "pdf_based" && has_rates) {
      gtr_pdf_samples <- list()
      for (var in rate_vars) {
        param_name <- paste0(var, "_pdf")
        if (!is.null(input$gtr_pdf_parameters) && param_name %in% input$gtr_pdf_parameters) {
          gtr_pdf_samples[[var]] <- generate_pdf_samples(var, num_sims, df)
        }
      }
    }
    
    # Pre-select rows for complete match strategy
    selected_rows <- NULL
    if (input$sim_strategy == "complete_match") {
      # Filter rows that have all required data
      required_vars <- character(0)
      if ("base_model" %in% input$model_components && seq_type == "DNA" && has_rates) {
        required_vars <- c(required_vars, rate_vars)
      }
      if ("frequencies" %in% input$model_components && has_frequencies) {
        required_vars <- c(required_vars, freq_vars)
      }
      if ("gamma" %in% input$model_components && has_gamma) {
        required_vars <- c(required_vars, "GammaShape_alpha")
      }
      if ("invariant" %in% input$model_components && has_invariant) {
        required_vars <- c(required_vars, "proportion_invariant")
      }
      if (has_num_sites) {
        required_vars <- c(required_vars, "num_sites")
      }
      if (has_loci) {
        required_vars <- c(required_vars, "loci")
      }
      
      # Find rows with complete data
      if (length(required_vars) > 0) {
        complete_rows <- complete.cases(df[, required_vars, drop = FALSE])
        valid_row_indices <- which(complete_rows)
        
        if (length(valid_row_indices) == 0) {
          return(list(commands = character(0), stats = "No complete rows found for complete empirical match strategy."))
        }
        
        # Sample rows with replacement
        selected_rows <- sample(valid_row_indices, num_sims, replace = TRUE)
      }
    }
    
    for (i in 1:num_sims) {
      # Generate sequence file name
      seq_name <- paste0(input$sim_prefix, sprintf("%04d", i))
      
      # Initialize command
      cmd <- paste0("iqtree3 --alisim ", seq_name, " --seqtype ", seq_type)
      
      # Build model string
      model_parts <- character(0)
      
      # For complete match, use the pre-selected row
      if (input$sim_strategy == "complete_match" && !is.null(selected_rows)) {
        row_idx <- selected_rows[i]
      }
      
      # Base model
      if ("base_model" %in% input$model_components) {
        if (seq_type == "DNA") {
          if (has_rates) {
            if (input$sim_strategy == "complete_match") {
              rates <- as.numeric(df[row_idx, rate_vars])
            } else if (input$sim_strategy == "pdf_based" && !is.null(input$gtr_source) && 
                       input$gtr_source == "pdf_based" && !is.null(gtr_pdf_samples)) {
              # GTR rates from PDF sampling
              rates <- numeric(length(rate_vars))
              for (j in seq_along(rate_vars)) {
                param_name <- paste0(rate_vars[j], "_pdf")
                if (!is.null(input$gtr_pdf_parameters) && param_name %in% input$gtr_pdf_parameters && 
                    !is.null(gtr_pdf_samples[[rate_vars[j]]])) {
                  rates[j] <- gtr_pdf_samples[[rate_vars[j]]][i]
                } else {
                  # Fall back to empirical sampling
                  available_values <- na.omit(df[[rate_vars[j]]])
                  if (length(available_values) > 0) {
                    rates[j] <- sample(available_values, 1)
                  } else {
                    rates[j] <- 1.0
                  }
                }
              }
            } else if (!is.null(input$gtr_source) && input$gtr_source == "same_sample") {
              # GTR rates from same sample
              available_rows <- which(complete.cases(df[, rate_vars, drop = FALSE]))
              if (length(available_rows) > 0) {
                rate_row <- sample(available_rows, 1)
                rates <- as.numeric(df[rate_row, rate_vars])
              } else {
                rates <- rep(1.0, length(rate_vars))  # Default values
              }
            } else {
              # Mixed empirical sampling for GTR rates
              rates <- numeric(length(rate_vars))
              for (j in seq_along(rate_vars)) {
                available_values <- na.omit(df[[rate_vars[j]]])
                if (length(available_values) > 0) {
                  rates[j] <- sample(available_values, 1)
                } else {
                  rates[j] <- 1.0  # Default value
                }
              }
            }
            rate_string <- paste(sapply(rates, format_number_clean), collapse = "/")
            model_parts <- c(model_parts, paste0("GTR{", rate_string, "}"))
          } else {
            model_parts <- c(model_parts, "GTR")
          }
        } else {
          # For AA sequences, LG model without brackets
          model_parts <- c(model_parts, "LG")
        }
      }
      
      # Frequencies
      if ("frequencies" %in% input$model_components && has_frequencies) {
        if (input$sim_strategy == "complete_match") {
          freqs <- as.numeric(df[row_idx, freq_vars])
        } else {
          # For mixed strategies, frequencies must come from same row to ensure they sum to 1
          available_rows <- which(complete.cases(df[, freq_vars, drop = FALSE]))
          if (length(available_rows) > 0) {
            freq_row <- sample(available_rows, 1)
            freqs <- as.numeric(df[freq_row, freq_vars])
          } else {
            freqs <- rep(1.0 / length(freq_vars), length(freq_vars))  # Equal frequencies as default
          }
        }
        freq_string <- paste(sapply(freqs, format_number_clean), collapse = "/")
        model_parts <- c(model_parts, paste0("F{", freq_string, "}"))
      }
      
      # Gamma shape
      if ("gamma" %in% input$model_components && has_gamma) {
        if (input$sim_strategy == "pdf_based" && "gamma_pdf" %in% input$pdf_parameters && !is.null(pdf_samples$gamma)) {
          gamma_val <- pdf_samples$gamma[i]
        } else if (input$sim_strategy == "complete_match") {
          gamma_val <- df[row_idx, "GammaShape_alpha"]
        } else {
          available_values <- na.omit(df[["GammaShape_alpha"]])
          if (length(available_values) > 0) {
            gamma_val <- sample(available_values, 1)
          } else {
            gamma_val <- 1.0  # Default value
          }
        }
        model_parts <- c(model_parts, paste0("G4{", format_number_clean(gamma_val), "}"))
      }
      
      # Invariant sites
      if ("invariant" %in% input$model_components && has_invariant) {
        if (input$sim_strategy == "pdf_based" && "invariant_pdf" %in% input$pdf_parameters && !is.null(pdf_samples$invariant)) {
          inv_val <- pdf_samples$invariant[i]
        } else if (input$sim_strategy == "complete_match") {
          inv_val <- df[row_idx, "proportion_invariant"]
        } else {
          available_values <- na.omit(df[["proportion_invariant"]])
          if (length(available_values) > 0) {
            inv_val <- sample(available_values, 1)
          } else {
            inv_val <- 0.0  # Default value
          }
        }
        model_parts <- c(model_parts, paste0("I{", format_number_clean(inv_val), "}"))
      }
      
      # Add model to command
      if (length(model_parts) > 0) {
        model_string <- paste(model_parts, collapse = "+")
        cmd <- paste0(cmd, " -m \"", model_string, "\"")
      }
      
      # Tree file
      if (input$sim_strategy == "complete_match") {
        if (has_loci) {
          tree_name <- paste0(df[row_idx, "loci"], ".tre")
        } else {
          tree_name <- "NEWICK_TREE"
        }
      } else {
        # For mixed strategies
        if (input$tree_source == "custom") {
          tree_name <- "NEWICK_TREE"
        } else if (has_loci) {
          available_loci <- na.omit(df[["loci"]])
          if (length(available_loci) > 0) {
            tree_name <- paste0(sample(available_loci, 1), ".tre")
          } else {
            tree_name <- "NEWICK_TREE"
          }
        } else {
          tree_name <- "NEWICK_TREE"
        }
      }
      cmd <- paste0(cmd, " -t ", tree_name)
      
      # Sequence length - check PDF first to avoid conflicts
      if (input$sim_strategy == "pdf_based" && "length_pdf" %in% input$pdf_parameters && !is.null(pdf_samples$length)) {
        length_val <- pdf_samples$length[i]
      } else if (input$sim_strategy == "complete_match") {
        if (has_num_sites) {
          length_val <- df[row_idx, "num_sites"]
          if (is.na(length_val)) {
            length_val <- 500  # Default
          }
        } else {
          length_val <- 500
        }
      } else if (input$length_source == "fixed") {
        length_val <- input$fixed_length
      } else if (has_num_sites) {
        available_lengths <- na.omit(df[["num_sites"]])
        if (length(available_lengths) > 0) {
          length_val <- sample(available_lengths, 1)
        } else {
          length_val <- 500
        }
      } else {
        length_val <- 500
      }
      cmd <- paste0(cmd, " --length ", as.integer(length_val))
      
      # Fixed parameters
      cmd <- paste0(cmd, " --num-alignments 1 -af fasta -nt ", input$num_threads)
      
      commands[i] <- cmd
    }
    
    # Generate statistics
    stats_text <- paste0(
      "Total simulations: ", num_sims, "\n",
      "Sequence type: ", seq_type, "\n",
      "Strategy: ", switch(input$sim_strategy,
                           "complete_match" = "Complete empirical match",
                           "mixed_empirical" = "Mixed empirical sampling", 
                           "pdf_based" = "Probability density function (PDF) estimation"), "\n",
      "Model components: ", paste(input$model_components, collapse = ", "), "\n",
      if (seq_type == "DNA" && "base_model" %in% input$model_components && input$sim_strategy != "complete_match") {
        paste("GTR parameter source: ", switch(input$gtr_source,
                                               "same_sample" = "From same sample",
                                               "mixed_empirical" = "Mixed empirical sampling",
                                               "pdf_based" = "PDF estimation"), "\n")
      } else "",
      "Available variables in data: ", 
      paste(c(
        if(has_gamma) "GammaShape_alpha",
        if(has_invariant) "proportion_invariant", 
        if(has_num_sites) "num_sites",
        if(has_loci) "loci",
        if(has_frequencies) "frequencies",
        if(has_rates) "rates"
      ), collapse = ", "), "\n",
      "Data rows available: ", nrow(df), "\n",
      if (input$sim_strategy == "complete_match" && !is.null(selected_rows)) {
        paste("Complete rows used: ", length(unique(selected_rows)), " unique rows")
      } else ""
    )
    
    return(list(commands = commands, stats = stats_text))
  })
  
  # Command preview output
  output$command_preview <- renderText({
    if (!input$enable_msa_simulation) {
      return("MSA simulation disabled")
    }
    
    sim_result <- generate_simulation_commands()
    
    if (length(sim_result$commands) == 0) {
      return("No commands generated. Please check your data and settings.")
    }
    
    # Show first few commands as preview
    preview_count <- min(5, length(sim_result$commands))
    preview_commands <- sim_result$commands[1:preview_count]
    
    preview_text <- paste0(
      "Preview of first ", preview_count, " commands:\n\n",
      paste(paste0(1:preview_count, ". ", preview_commands), collapse = "\n"),
      if (length(sim_result$commands) > preview_count) {
        paste0("\n\n... and ", length(sim_result$commands) - preview_count, " more commands")
      } else ""
    )
    
    return(preview_text)
  })
  
  # Simulation statistics output
  output$simulation_stats <- renderText({
    if (!input$enable_msa_simulation) {
      return("")
    }
    
    sim_result <- generate_simulation_commands()
    return(sim_result$stats)
  })
  
  # Download handlers for original functionality
  output$download_filtered_data <- downloadHandler(
    filename = function() {
      paste0("filtered_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      # Make sure to use the fully filtered data (including Tukey's filtering)
      filtered_df <- filtered_data()
      write_csv(filtered_df, file)
    }
  )
  
  output$download_plot_pdf <- downloadHandler(
    filename = function() {
      paste0("plot_", input$selected_variable, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      ggsave(file, plot = create_plot(), device = "pdf", 
             width = input$plot_width, height = input$plot_height, units = "in")
    }
  )
  
  # Download handlers for correlation functionality
  output$download_correlation_matrix <- downloadHandler(
    filename = function() {
      paste0("correlation_matrix_", Sys.Date(), ".csv")
    },
    content = function(file) {
      cor_matrix <- correlation_matrix()
      req(cor_matrix)
      
      # Convert matrix to data frame with row names as first column
      cor_df <- as.data.frame(cor_matrix)
      cor_df <- cbind(Variable = rownames(cor_matrix), cor_df)
      
      write_csv(cor_df, file)
    }
  )
  
  output$download_correlation_plot <- downloadHandler(
    filename = function() {
      paste0("correlation_heatmap_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, width = input$heatmap_width, height = input$heatmap_height)
      create_correlation_plot()
      dev.off()
    }
  )
  
  # Download handler for PDF comparison plot
  output$download_pdf_comparison_plot <- downloadHandler(
    filename = function() {
      paste0("pdf_comparison_plot_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      plot_obj <- create_pdf_comparison_plot()
      
      if (inherits(plot_obj, "grob")) {
        # If it's a grob (from gridExtra), save using ggsave
        ggsave(file, plot = plot_obj, device = "pdf", 
               width = input$pdf_plot_width, height = input$pdf_plot_height, units = "in")
      } else {
        # If it's a ggplot, save using ggsave
        ggsave(file, plot = plot_obj, device = "pdf", 
               width = input$pdf_plot_width, height = input$pdf_plot_height, units = "in")
      }
    }
  )
  
  # Download handler for simulation commands
  output$download_simulation_commands <- downloadHandler(
    filename = function() {
      paste0("iqtree_simulation_commands_", Sys.Date(), ".txt")
    },
    content = function(file) {
      sim_result <- generate_simulation_commands()
      
      if (length(sim_result$commands) == 0) {
        writeLines("No commands generated. Please check your data and settings.", file)
        return()
      }
      
      # Write header with settings information
      header <- c(
        "# IQ-TREE MSA Simulation Commands",
        paste("# Generated on:", Sys.time()),
        paste("# Total simulations:", input$num_alignments),
        paste("# Sequence type:", toupper(input$var_type)),
        paste("# Strategy:", input$sim_strategy),
        paste("# Model components:", paste(input$model_components, collapse = ", ")),
        "#",
        "# NOTES:",
        "# 1. Indel Settings: Generated MSAs will be gapless by default.",
        if (input$sim_strategy == "complete_match") {
          c("# 2. For Complete empirical match mode: You can use scripts to transfer", 
            "#    gap patterns from original FASTA sequences to simulated MSAs.")
        } else {
          c("# 2. For mixed/PDF strategies: Add --indel and --indel-size parameters", 
            "#    manually for gapped MSA simulation.")
        },
        "# 3. Command Execution: These commands can be modified as needed and executed",
        "#    through loops or parallel processing for efficient batch processing.",
        "#"
      )
      
      header <- c(header, "# Commands:")
      
      writeLines(c(header, "", sim_result$commands), file)
    }
  )
}

# Run the app ----
shinyApp(ui = ui, server = server)
