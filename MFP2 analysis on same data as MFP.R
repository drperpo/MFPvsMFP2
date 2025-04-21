# Script: mfp2_analysis_load_subsets.R
# Description: Runs mfp2 models on pre-existing data subsets loaded from an .rds file.
#              This ensures the analysis uses the exact same data as a previous run
#              (e.g., one using the 'mfp' package).
# Version: 2 (Corrected formula handling for mfp2 call)

# --- 1. Load Libraries ---
# Ensure these packages are installed: install.packages(c("mfp2", "survival", "dplyr", "reshape2", "ggplot2"))
library(mfp2)
library(survival)
library(dplyr)    # Used in results extraction potentially
library(reshape2) # Used in results extraction
library(ggplot2)  # Potentially used if plotting is added later

message("Libraries loaded.")

# Set seed (optional, but good practice if any random operations were added later)
set.seed(123)

# --- 2. Load Pre-existing Data Subsets ---
subset_file <- "analysisMFP2/data/dat_subsets_list.rds" # The file saved by the first script ('New analysis mfp.R')

message(paste("Attempting to load data subsets from:", subset_file))
if (file.exists(subset_file)) {
  tryCatch({
    dat_subsets <- readRDS(subset_file)
    message("Successfully loaded data subsets list.")
    # Optional: Check the names and dimensions of loaded subsets
    print("Loaded subset names:")
    print(names(dat_subsets))
    # Example check on one subset:
    # if("dat_10K_1" %in% names(dat_subsets)) {
    #   print(paste("Dimensions of dat_10K_1:", nrow(dat_subsets$dat_10K_1), "rows,", ncol(dat_subsets$dat_10K_1), "columns"))
    # }
  }, error = function(e) {
    stop("Error loading ", subset_file, ". Please ensure the file exists and is a valid RDS file. Original error: ", e$message)
  })
} else {
  stop(paste("Data subset file not found:", subset_file, ". Please run the script that generates this file first ('New analysis mfp.R')."))
}

# --- 3. Define Model Structure (No separate formula variable needed for mfp2 call) ---
# The formula will be specified directly in the mfp2() function call below.
# Ensure the loaded data has the correct column types (factors, numeric).
# Preprocessing (grade dummies, factor levels) should be done *before* saving the .rds file.
message("Model structure defined (formula specified directly in mfp2 calls).")

# --- 4. Run mfp2 Models (AIC-like and BIC) ---
message("Fitting mfp2 models using loaded subsets...")

# Define model names corresponding to the expected subset names in the loaded list
# This should match the names used when the list was created.
model_names_map <- list(
  "M1" = "dat_Full", # Assumes the full dataset was saved as 'dat_Full' in the .rds file
  "M200K1" = "dat_200K_1", "M200K2" = "dat_200K_2",
  "M100K1" = "dat_100K_1", "M100K2" = "dat_100K_2",
  "M50K1" = "dat_50K_1", "M50K2" = "dat_50K_2",
  "M25K1" = "dat_25K_1", "M25K2" = "dat_25K_2", "M25K3" = "dat_25K_3",
  "M10K1" = "dat_10K_1", "M10K2" = "dat_10K_2", "M10K3" = "dat_10K_3",
  "M5K1" = "dat_5K_1", "M5K2" = "dat_5K_2", "M5K3" = "dat_5K_3",
  "M2K5_1" = "dat_2K5_1", "M2K5_2" = "dat_2K5_2", "M2K5_3" = "dat_2K5_3", # Adjusted name if needed
  "M1K1" = "dat_1K_1", "M1K2" = "dat_1K_2", "M1K3" = "dat_1K_3"
)

# Store models in lists
aic_models_mfp2 <- list()
bic_models_mfp2 <- list()

# --- 4a. AIC-like Models (select = 0.157) ---
message("\nFitting AIC-like models (select=0.157)...")
aic_select_level <- 0.157

for (model_name in names(model_names_map)) {
  data_name <- model_names_map[[model_name]]
  if (!data_name %in% names(dat_subsets)) {
    warning(paste("Data subset", data_name, "not found in the loaded list for model", model_name, ". Skipping."))
    next
  }
  current_data <- dat_subsets[[data_name]]
  
  # Define the full formula temporarily for checking required variables
  temp_formula_check <- Surv(futime, Death) ~ fp(AgeDgc, df = 4) + fp(tSizeMm, df = 4) + fp(EodNPos, df = 4) + YearDgc + grade2_3 + grade3 + factor(rtsurg) + factor(racegrp) + factor(Lateral) + factor(ER_imp) + factor(PR_imp) + factor(married) + factor(innerkwd)
  required_vars <- all.vars(temp_formula_check)
  
  # Basic check for required columns in the subset
  if (!all(required_vars %in% names(current_data))) {
    missing_cols <- setdiff(required_vars, names(current_data))
    warning(paste("Subset", data_name, "is missing required columns:", paste(missing_cols, collapse=", "), ". Skipping model", model_name))
    next
  }
  
  message(paste("Fitting", model_name, "(AIC) using loaded data", data_name, "(", nrow(current_data), "rows )"))
  aic_models_mfp2[[model_name]] <- tryCatch({
    # Specify the formula directly in the mfp2 call
    mfp2(Surv(futime, Death) ~
           fp(AgeDgc, df = 4) +
           fp(tSizeMm, df = 4) +
           fp(EodNPos, df = 4) +
           YearDgc +          # Initially linear, mfp2 can select FP
           grade2_3 + grade3 + # Assumed numeric 0/1 from preprocessing
           factor(rtsurg) +   # Explicitly factor
           factor(racegrp) +  # Explicitly factor (ensure ref level is correct in loaded data)
           factor(Lateral) +  # Explicitly factor
           factor(ER_imp) +   # Explicitly factor
           factor(PR_imp) +   # Explicitly factor
           factor(married) +  # Explicitly factor
           factor(innerkwd),  # Explicitly factor
         data = current_data,
         family = "cox",
         select = aic_select_level,
         verbose = FALSE) # Set verbose=TRUE for more details if needed
  }, error = function(e) {
    warning(paste("Error fitting", model_name, "(AIC):", e$message)); return(NULL)
  })
  if (!is.null(aic_models_mfp2[[model_name]])) { message(paste("Finished", model_name, "(AIC)")) }
}


# --- 4b. BIC Models (select based on events) ---
message("\nFitting BIC models (select = pchisq(log(events), 1, lower.tail = FALSE))...")

for (model_name in names(model_names_map)) {
  model_name_bic <- paste0(model_name, "b") # Add suffix for BIC model name
  data_name <- model_names_map[[model_name]]
  
  if (!data_name %in% names(dat_subsets)) {
    # Warning given in AIC loop
    next
  }
  current_data <- dat_subsets[[data_name]]
  
  # Define the full formula temporarily for checking required variables
  temp_formula_check <- Surv(futime, Death) ~ fp(AgeDgc, df = 4) + fp(tSizeMm, df = 4) + fp(EodNPos, df = 4) + YearDgc + grade2_3 + grade3 + factor(rtsurg) + factor(racegrp) + factor(Lateral) + factor(ER_imp) + factor(PR_imp) + factor(married) + factor(innerkwd)
  required_vars <- all.vars(temp_formula_check)
  
  # Check required vars again (could be skipped if AIC loop passed)
  if (!all(required_vars %in% names(current_data))) {
    # Warning given in AIC loop
    next
  }
  
  # Calculate BIC-based selection level (slc) based on events in THIS subset
  if (!"Death" %in% names(current_data) || !is.numeric(current_data$Death)) {
    warning(paste("Outcome variable 'Death' not found or not numeric in subset", data_name, "for model", model_name_bic, ". Skipping."))
    next
  }
  events <- sum(current_data$Death, na.rm = TRUE)
  if (events <= 0) {
    warning(paste("No events (Death > 0) found in", data_name, "for model", model_name_bic, ". Skipping."))
    next
  }
  slc <- round(pchisq(log(events), 1, lower.tail = FALSE), 4)
  
  message(paste("Fitting", model_name_bic, "using loaded data", data_name, "(", nrow(current_data), "rows,", events, "events) with BIC criterion (select=", slc, ")"))
  bic_models_mfp2[[model_name_bic]] <- tryCatch({
    # Specify the formula directly in the mfp2 call
    mfp2(Surv(futime, Death) ~
           fp(AgeDgc, df = 4) +
           fp(tSizeMm, df = 4) +
           fp(EodNPos, df = 4) +
           YearDgc +
           grade2_3 + grade3 +
           factor(rtsurg) +
           factor(racegrp) +
           factor(Lateral) +
           factor(ER_imp) +
           factor(PR_imp) +
           factor(married) +
           factor(innerkwd),
         data = current_data,
         family = "cox",
         select = slc, # Use calculated BIC-based selection level
         verbose = FALSE)
  }, error = function(e) {
    warning(paste("Error fitting", model_name_bic, ":", e$message)); return(NULL)
  })
  if (!is.null(bic_models_mfp2[[model_name_bic]])) { message(paste("Finished", model_name_bic)) }
}
message("Finished all mfp2 model fitting.")

# --- 5. Extract and Save mfp2 Results Tables ---
message("\nGenerating results tables for mfp2 models...")

# Use the same extraction function as in the original 'new analysis mfp2.R' script.
# Ensure this function definition is available here.
# (Copying the function from the provided 'new analysis mfp2.R' content)
# --- Function to Extract Results from mfp2 Models (REVISED V4 - Fix Init) ---
extract_mfp2_results <- function(model_list) {
  # Identify potential continuous variables based on typical usage in formula
  # fp() terms + linear terms that *could* become fp()
  fp_input_vars <- c("AgeDgc", "tSizeMm", "EodNPos")
  linear_vars <- c("YearDgc")
  all_cont_vars <- c(fp_input_vars, linear_vars)
  
  fp_summary_list <- list()    # Store results for FP/Linear terms
  cat_summary_list <- list()   # Store results for Categorical terms
  
  valid_model_names <- names(model_list)[!sapply(model_list, is.null)]
  if (length(valid_model_names) == 0) {
    warning("No valid models found in the provided list.")
    return(list(fp_powers = NULL, cat_selection = NULL)) # Return empty structure
  }
  
  # 1. Iterate through each valid model
  for (model_name in valid_model_names) {
    model <- model_list[[model_name]]
    
    # Check essential components exist
    if (is.null(model$fp_terms) || is.null(model$coefficients)) {
      warning(paste("Skipping model", model_name, "- missing fp_terms or coefficients."))
      next
    }
    
    fp_terms_table <- model$fp_terms
    coef_names <- names(model$coefficients)
    
    # Initialize empty dataframes *for this model* with correct column names
    fp_summary_for_model <- data.frame(model = character(0), term = character(0), DF = numeric(0), P1 = numeric(0), P2 = numeric(0), stringsAsFactors = FALSE)
    cat_summary_for_model <- data.frame(model = character(0), term = character(0), Selected = character(0), stringsAsFactors = FALSE)
    
    # 2. Process each term listed in fp_terms table
    for (term_name in rownames(fp_terms_table)) {
      term_info <- fp_terms_table[term_name, ]
      is_selected <- term_info$selected # TRUE/FALSE
      
      # --- Handle Continuous Variables (FP or Linear) ---
      # Check if the term name (e.g., "AgeDgc") is one of the expected continuous vars
      if (term_name %in% all_cont_vars) {
        p1 <- if (is_selected) term_info$power1 else NA
        # P2 is only relevant if df > 2 (i.e., FP2)
        p2 <- if (is_selected && term_info$df_final > 2) term_info$power2 else NA
        df <- if (is_selected) term_info$df_final else 0 # 0 DF if not selected (linear or removed)
        
        new_fp_row <- data.frame(model = model_name, term = term_name, DF = df, P1 = p1, P2 = p2, stringsAsFactors = FALSE)
        fp_summary_for_model <- rbind(fp_summary_for_model, new_fp_row)
        
        # --- Handle Categorical Variables ---
        # Assumes categorical terms appear like "factor(rtsurg)" or just "grade2_3" in rownames
      } else {
        # Assume anything not in all_cont_vars is categorical or a simple numeric term like grade2_3
        selection_status <- ifelse(is_selected, "+", "-")
        
        # Clean up term name (remove factor() wrapper if present)
        clean_term_name <- gsub("factor\\((.*?)\\)", "\\1", term_name)
        
        new_cat_row <- data.frame(model = model_name, term = clean_term_name, Selected = selection_status, stringsAsFactors = FALSE)
        cat_summary_for_model <- rbind(cat_summary_for_model, new_cat_row)
      }
    } # End loop through terms in fp_terms table
    
    # Store the completed dataframes for this model in the main lists
    fp_summary_list[[model_name]] <- fp_summary_for_model
    cat_summary_list[[model_name]] <- cat_summary_for_model
    
  } # End loop through models
  
  # 3. Combine and Reshape FP Summaries
  fp_results_wide <- NULL
  if (length(fp_summary_list) > 0) {
    fp_summary_long <- do.call(rbind, fp_summary_list)
    rownames(fp_summary_long) <- NULL # Reset row names after rbind
    
    if (!is.null(fp_summary_long) && nrow(fp_summary_long) > 0) {
      # Reshape using dcast
      fp_p1_wide <- dcast(fp_summary_long, model ~ term, value.var = "P1", fill = NA)
      fp_p2_wide <- dcast(fp_summary_long, model ~ term, value.var = "P2", fill = NA)
      fp_df_wide <- dcast(fp_summary_long, model ~ term, value.var = "DF", fill = 0) # Fill DF with 0 if term missing
      
      # Merge the reshaped dataframes
      fp_results_wide <- merge(fp_p1_wide, fp_p2_wide, by = "model", suffixes = c("_P1", "_P2"), all = TRUE)
      fp_results_wide <- merge(fp_results_wide, fp_df_wide, by = "model", all = TRUE) # Suffix for DF cols added below
      
      # Select and rename columns for final output
      final_fp_cols <- c("model")
      for (term in all_cont_vars) {
        col_p1 <- paste0(term, "_P1"); col_p2 <- paste0(term, "_P2"); col_df <- term; col_df_final <- paste0(term, "_DF")
        
        # Check if columns exist after merge before adding to selection
        if (col_p1 %in% names(fp_results_wide)) final_fp_cols <- c(final_fp_cols, col_p1)
        if (col_p2 %in% names(fp_results_wide)) final_fp_cols <- c(final_fp_cols, col_p2)
        if (col_df %in% names(fp_results_wide)) {
          names(fp_results_wide)[names(fp_results_wide) == col_df] <- col_df_final # Rename DF column
          final_fp_cols <- c(final_fp_cols, col_df_final)
        }
      }
      # Ensure only existing columns are selected
      final_fp_cols_exist <- intersect(final_fp_cols, names(fp_results_wide))
      fp_results_wide <- fp_results_wide[, final_fp_cols_exist, drop = FALSE]
      
    } else { warning("FP/Linear summary data frame is empty after combining.") }
  } else { warning("No FP/Linear summary results generated.") }
  
  
  # 4. Combine and Reshape Categorical Summaries
  cat_selection_wide <- NULL
  if (length(cat_summary_list) > 0) {
    cat_summary_long <- do.call(rbind, cat_summary_list)
    rownames(cat_summary_long) <- NULL # Reset row names
    
    if (!is.null(cat_summary_long) && nrow(cat_summary_long) > 0) {
      # Ensure Selection is only +/-
      cat_summary_long$Selected <- ifelse(cat_summary_long$Selected == "+", "+", "-")
      # Reshape using dcast
      cat_selection_wide <- dcast(cat_summary_long, model ~ term, value.var = "Selected", fill = "-") # Fill with '-' if term missing
      
      # Apply renaming rules (same as original script)
      current_names <- names(cat_selection_wide)
      new_names <- current_names
      new_names <- gsub("^rtsurg", "Surg_", new_names); new_names <- gsub("^racegrp", "Race_", new_names)
      new_names <- gsub("^ER_imp", "ER_", new_names); new_names <- gsub("^PR_imp", "PR_", new_names)
      new_names <- gsub("^married", "Married_", new_names); new_names <- gsub("^innerkwd", "InnerKwd_", new_names)
      new_names <- gsub("^Lateral", "Lateral_", new_names); new_names <- gsub("^grade2_3", "Grade_2or3", new_names); new_names <- gsub("^grade3", "Grade_3", new_names)
      # Handle levels (assuming they appear after the base name)
      new_names <- gsub("Positive$", "Pos", new_names); new_names <- gsub("bcrt$", "BCRT", new_names)
      new_names <- gsub("mart$", "MART", new_names); new_names <- gsub("ma$", "MA", new_names)
      new_names <- gsub("black$", "Black", new_names); new_names <- gsub("asian$", "Asian", new_names)
      new_names <- gsub("whisp$", "Hisp", new_names) # Adjust if needed
      # Add underscore before numeric levels if not already present (e.g., Lateral2 -> Lateral_2)
      new_names <- gsub("([a-zA-Z])(\\d+)$", "\\1_\\2", new_names)
      names(cat_selection_wide) <- new_names
      message("Applied renaming to categorical columns:")
      print(names(cat_selection_wide))
      
    } else { warning("Categorical summary data frame is empty after combining.") }
  } else { warning("No categorical summary results generated.") }
  
  
  # 5. Reorder rows based on the original list of models attempted
  expected_model_order <- names(model_list) # Use the names from the input list
  
  if (!is.null(fp_results_wide) && nrow(fp_results_wide) > 0) {
    # Ensure 'model' column is character for matching, then factor for ordering
    fp_results_wide$model <- as.character(fp_results_wide$model)
    # Match requires the column to exist
    if("model" %in% names(fp_results_wide)){
      fp_results_wide <- fp_results_wide[match(expected_model_order, fp_results_wide$model), ]
      # Remove rows that are all NA (might occur if a model failed completely but name was in list)
      fp_results_wide <- fp_results_wide[rowSums(is.na(fp_results_wide)) < ncol(fp_results_wide), ]
      fp_results_wide$model <- factor(fp_results_wide$model, levels = intersect(expected_model_order, fp_results_wide$model)) # Keep order
      rownames(fp_results_wide) <- NULL
    } else {
      warning("Column 'model' not found in fp_results_wide for reordering.")
    }
  }
  
  if (!is.null(cat_selection_wide) && nrow(cat_selection_wide) > 0) {
    cat_selection_wide$model <- as.character(cat_selection_wide$model)
    if("model" %in% names(cat_selection_wide)){
      cat_selection_wide <- cat_selection_wide[match(expected_model_order, cat_selection_wide$model), ]
      cat_selection_wide <- cat_selection_wide[rowSums(is.na(cat_selection_wide)) < ncol(cat_selection_wide), ]
      cat_selection_wide$model <- factor(cat_selection_wide$model, levels = intersect(expected_model_order, cat_selection_wide$model))
      rownames(cat_selection_wide) <- NULL
    } else {
      warning("Column 'model' not found in cat_selection_wide for reordering.")
    }
  }
  
  # 6. Return results
  return(list(fp_powers = fp_results_wide, cat_selection = cat_selection_wide))
}


# Extract and Save Results using the function
aic_results_mfp2 <- extract_mfp2_results(aic_models_mfp2)
bic_results_mfp2 <- extract_mfp2_results(bic_models_mfp2)

# Define output file names
aic_fp_file <- "mfp2_aic_fp_powers_summary_loaded.csv"
aic_cat_file <- "mfp2_aic_categorical_selection_summary_loaded.csv"
bic_fp_file <- "mfp2_bic_fp_powers_summary_loaded.csv"
bic_cat_file <- "mfp2_bic_categorical_selection_summary_loaded.csv"

# Save results if they exist
if (!is.null(aic_results_mfp2$fp_powers) && nrow(aic_results_mfp2$fp_powers) > 0) {
  write.csv(aic_results_mfp2$fp_powers, file = aic_fp_file, row.names = FALSE)
  message(paste("AIC FP powers summary saved to:", aic_fp_file))
} else { message("AIC FP powers results were empty or null, not saved.") }

if (!is.null(aic_results_mfp2$cat_selection) && nrow(aic_results_mfp2$cat_selection) > 0) {
  write.csv(aic_results_mfp2$cat_selection, file = aic_cat_file, row.names = FALSE)
  message(paste("AIC categorical selection summary saved to:", aic_cat_file))
} else { message("AIC categorical selection results were empty or null, not saved.") }

if (!is.null(bic_results_mfp2$fp_powers) && nrow(bic_results_mfp2$fp_powers) > 0) {
  write.csv(bic_results_mfp2$fp_powers, file = bic_fp_file, row.names = FALSE)
  message(paste("BIC FP powers summary saved to:", bic_fp_file))
} else { message("BIC FP powers results were empty or null, not saved.") }

if (!is.null(bic_results_mfp2$cat_selection) && nrow(bic_results_mfp2$cat_selection) > 0) {
  write.csv(bic_results_mfp2$cat_selection, file = bic_cat_file, row.names = FALSE)
  message(paste("BIC categorical selection summary saved to:", bic_cat_file))
} else { message("BIC categorical selection results were empty or null, not saved.") }

message("Script finished.")

# Optional: Clean up environment if needed
# rm(list = ls())
