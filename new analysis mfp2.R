# Script: mfp2_analysis_pipeline.R
# Replicates the original multi-subset mfp analysis using only the mfp2 package.

# --- 1. Load Libraries and Data ---
# install.packages(c("mfp2", "survival", "dplyr", "reshape2", "ggplot2", "corrplot")) # Run once if needed
library(mfp2)
library(survival)
library(dplyr)    # For glimpse, etc.
library(reshape2) # For dcast later
library(ggplot2)  # For plotting predictions
library(corrplot) # For correlation plot

# Set seed for reproducible sampling
set.seed(123)

# Load the main dataset
tryCatch({
  dat <- read.csv("data/dat_ida.csv", header = TRUE)
  message("Successfully loaded data/dat_ida.csv")
}, error = function(e) {
  stop("Error loading data/dat_ida.csv. Original error: ", e$message)
})

# Basic check
print(paste("Original data dimensions:", nrow(dat), "rows,", ncol(dat), "columns"))
required_cols <- c("Grade_imp", "futime", "Death", "AgeDgc", "tSizeMm", "EodNPos",
                   "YearDgc", "rtsurg", "racegrp", "Lateral", "ER_imp",
                   "PR_imp", "married", "innerkwd")
if (!all(required_cols %in% names(dat))) {
  stop("Missing required columns in the dataset. Needed: ", paste(required_cols, collapse=", "))
}

# --- 2. Preprocessing (Applied ONCE to the main dataframe) ---
# Although mfp2 uses factor() in formula, preprocessing ensures correct levels/types
message("Starting preprocessing...")
dat$grade2_3 <- ifelse(dat$Grade_imp > 1, 1, 0)
dat$grade3   <- ifelse(dat$Grade_imp == 3, 1, 0)

# Convert relevant columns to factors (essential for releveling and consistency)
factor_cols <- c("rtsurg", "racegrp", "Lateral", "ER_imp", "PR_imp", "married", "innerkwd")
for (col in factor_cols) {
  if (col %in% names(dat)) {
    dat[[col]] <- factor(dat[[col]])
  } else { warning(paste("Column", col, "not found for factor conversion.")) }
}

# Relevel racegrp (ensure 'white' exists as a level)
if ("racegrp" %in% names(dat) && "white" %in% levels(dat$racegrp)) {
  dat$racegrp <- relevel(dat$racegrp, ref = "white")
} else { warning("'racegrp' or level 'white' not found. Skipping releveling.") }
message("Preprocessing finished.")

# --- 3. Create Data Subsets and Store in a List ---
message("Creating data subsets...")
n_total <- nrow(dat)
sample_indices <- 1:n_total

# Define sample sizes and replications
sample_defs <- list(
  "Full" = list(size = n_total, reps = 1), # Add Full dataset explicitly
  "200K" = list(size = 200000, reps = 2),
  "100K" = list(size = 100000, reps = 2),
  "50K"  = list(size = 50000,  reps = 2),
  "25K"  = list(size = 25000,  reps = 3),
  "10K"  = list(size = 10000,  reps = 3),
  "5K"   = list(size = 5000,   reps = 3),
  "2K5"  = list(size = 2500,   reps = 3), # Use K for clarity
  "1K"   = list(size = 1000,   reps = 3)
)

# Store subsets in a list
dat_subsets <- list()

for (name in names(sample_defs)) {
  size <- sample_defs[[name]]$size
  reps <- sample_defs[[name]]$reps
  if (size > n_total) {
    warning(paste("Sample size", size, "is larger than total data", n_total, ". Skipping subset", name))
    next
  }
  if (name == "Full") {
    subset_name <- "dat_Full"
    dat_subsets[[subset_name]] <- dat
    message(paste("Added full dataset as:", subset_name))
  } else {
    for (i in 1:reps) {
      subset_name <- paste0("dat_", name, "_", i) # e.g., dat_200K_1
      sampled_rows <- sample(sample_indices, size, replace = FALSE)
      dat_subsets[[subset_name]] <- dat[sampled_rows, ]
      message(paste("Created subset:", subset_name, "with", size, "rows"))
    }
  }
}

# Optional: Save subsets for later use
# message("Saving the list of data subsets...")
saveRDS(dat_subsets, file = "dat_subsets_list_for_mfp2.rds")

# --- 4. Define Model Formula for mfp2 ---
# Using the syntax confirmed to work by the user
# Identify variables needed for NA checks
#mfp2_model_vars <- all.vars(mfp2_formula)

# --- 5. Run mfp2 Models (AIC-like and BIC) ---
message("Fitting mfp2 models for all subsets...")

# Define model names corresponding to subsets
# Ensure this matches the names generated in the dat_subsets list
model_names_map <- list(
  "M1" = "dat_Full",
  "M200K1" = "dat_200K_1", "M200K2" = "dat_200K_2",
  "M100K1" = "dat_100K_1", "M100K2" = "dat_100K_2",
  "M50K1" = "dat_50K_1", "M50K2" = "dat_50K_2",
  "M25K1" = "dat_25K_1", "M25K2" = "dat_25K_2", "M25K3" = "dat_25K_3",
  "M10K1" = "dat_10K_1", "M10K2" = "dat_10K_2", "M10K3" = "dat_10K_3",
  "M5K1" = "dat_5K_1", "M5K2" = "dat_5K_2", "M5K3" = "dat_5K_3",
  "M2K5_1" = "dat_2K5_1", "M2K5_2" = "dat_2K5_2", "M2K5_3" = "dat_2K5_3",
  "M1K1" = "dat_1K_1", "M1K2" = "dat_1K_2", "M1K3" = "dat_1K_3"
)

# Store models in lists
aic_models_mfp2 <- list()
bic_models_mfp2 <- list()

# --- AIC-like Models (select = 0.157) ---
message("\nFitting AIC-like models (select=0.157)...")
aic_select_level <- 0.157

for (model_name in names(model_names_map)) {
  data_name <- model_names_map[[model_name]]
  current_data <- dat_subsets[[data_name]]
  
  message(paste("Fitting", model_name, "(AIC) using data", data_name))
  aic_models_mfp2[[model_name]] <- tryCatch({
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
           factor(innerkwd) ,
         data = current_data, # Use complete data
         family = "cox",
         select = aic_select_level,
         verbose = FALSE)
  }, error = function(e) {
    warning(paste("Error fitting", model_name, "(AIC):", e$message)); return(NULL)
  })
  if (!is.null(aic_models_mfp2[[model_name]])) { message(paste("Finished", model_name, "(AIC)")) }
}


# --- BIC Models (criterion = "BIC") ---
message("\nFitting BIC models (criterion='BIC')...")

for (model_name in names(model_names_map)) {
  model_name_bic <- paste0(model_name, "b") # Add suffix for BIC model name
  data_name <- model_names_map[[model_name]]
  current_data <- dat_subsets[[data_name]]
  
  
  # Check events in complete data
  events <- sum(current_data$Death, na.rm = TRUE)
  slc <- round(pchisq(log(events), 1, lower.tail = FALSE), 4)
  
  message(paste("Fitting", model_name_bic, "using data", data_name))
  bic_models_mfp2[[model_name_bic]] <- tryCatch({
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
           factor(innerkwd) ,
         data = current_data, # Use complete data
         family = "cox",
         select = slc,
         verbose = FALSE)
  }, error = function(e) {
    warning(paste("Error fitting", model_name_bic, ":", e$message)); return(NULL)
  })
  if (!is.null(bic_models_mfp2[[model_name_bic]])) { message(paste("Finished", model_name_bic)) }
}
message("Finished all mfp2 model fitting.")

# --- 6. Extract and Save mfp2 Results Tables ---
message("\nGenerating results tables for mfp2 models...")

# Function to extract results from a list of mfp2 models
# --- Function to Extract Results from mfp2 Models ---

# --- Function to Extract Results from mfp2 Models (REVISED V4 - Fix Init) ---

extract_mfp2_results <- function(model_list) {
  fp_input_vars <- c("AgeDgc", "tSizeMm", "EodNPos")
  linear_vars <- c("YearDgc") # Initially linear, can become FP
  all_cont_vars <- c(fp_input_vars, linear_vars)
  
  fp_summary_list <- list()    # Store results for FP/Linear terms
  cat_summary_list <- list()   # Store results for Categorical terms
  
  valid_model_names <- names(model_list)[!sapply(model_list, is.null)]
  if (length(valid_model_names) == 0) {
    warning("No valid models found in the provided list.")
    return(NULL)
  }
  
  # 1. Iterate through each valid model
  for (model_name in valid_model_names) {
    model <- model_list[[model_name]]
    
    if (is.null(model$fp_terms) || is.null(model$coefficients)) {
      warning(paste("Skipping model", model_name, "- missing fp_terms or coefficients."))
      next
    }
    
    fp_terms_table <- model$fp_terms
    coef_names <- names(model$coefficients)
    
    # *** FIXED: Initialize empty dataframes with correct column names ***
    fp_summary_for_model <- data.frame(model = character(0), term = character(0), DF = numeric(0), P1 = numeric(0), P2 = numeric(0))
    cat_summary_for_model <- data.frame(model = character(0), term = character(0), Selected = character(0))
    # *** End Fix ***
    
    # 2. Process each term listed in fp_terms
    for (term_name in rownames(fp_terms_table)) {
      term_info <- fp_terms_table[term_name, ]
      is_selected <- term_info$selected
      
      # --- Handle Continuous Variables (FP or Linear) ---
      if (term_name %in% all_cont_vars) {
        p1 <- ifelse(is_selected, term_info$power1, NA)
        p2 <- ifelse(is_selected && term_info$df_final > 2, term_info$power2, NA)
        df <- ifelse(is_selected, term_info$df_final, 0)
        
        # Create a row to add
        new_fp_row <- data.frame(model = model_name, term = term_name, DF = df, P1 = p1, P2 = p2)
        # Append to the dataframe *for this model*
        fp_summary_for_model <- rbind(fp_summary_for_model, new_fp_row)
        
        # --- Handle Categorical Variables ---
      } else {
        selection_status <- ifelse(is_selected, "+", "-")
        clean_term_name <- gsub("factor\\((.*?)\\)(.*)", "\\1\\2", term_name)
        
        new_cat_row <- data.frame(model = model_name, term = clean_term_name, Selected = selection_status)
        # Append to the dataframe *for this model*
        cat_summary_for_model <- rbind(cat_summary_for_model, new_cat_row)
      }
    } # End loop through terms
    
    # Store the completed dataframes for this model in the lists
    fp_summary_list[[model_name]] <- fp_summary_for_model
    cat_summary_list[[model_name]] <- cat_summary_for_model
    
  } # End loop through models
  
  # 3. Combine and Reshape FP Summaries (Keep as before)
  fp_summary_long <- do.call(rbind, fp_summary_list)
  fp_results_wide <- NULL
  if (!is.null(fp_summary_long) && nrow(fp_summary_long) > 0) {
    fp_p1_wide <- dcast(fp_summary_long, model ~ term, value.var = "P1", fill = NA)
    fp_p2_wide <- dcast(fp_summary_long, model ~ term, value.var = "P2", fill = NA)
    fp_df_wide <- dcast(fp_summary_long, model ~ term, value.var = "DF", fill = 0)
    fp_results_wide <- merge(fp_p1_wide, fp_p2_wide, by = "model", suffixes = c("_P1", "_P2"), all = TRUE)
    fp_results_wide <- merge(fp_results_wide, fp_df_wide, by = "model", all = TRUE)
    
    final_fp_cols <- c("model")
    for (term in all_cont_vars) {
      col_p1 <- paste0(term, "_P1"); col_p2 <- paste0(term, "_P2"); col_df <- term; col_df_final <- paste0(term, "_DF")
      if (col_p1 %in% names(fp_results_wide)) final_fp_cols <- c(final_fp_cols, col_p1)
      if (col_p2 %in% names(fp_results_wide)) final_fp_cols <- c(final_fp_cols, col_p2)
      if (col_df %in% names(fp_results_wide)) {
        names(fp_results_wide)[names(fp_results_wide) == col_df] <- col_df_final
        final_fp_cols <- c(final_fp_cols, col_df_final)
      }
    }
    final_fp_cols_exist <- final_fp_cols[final_fp_cols %in% names(fp_results_wide)]
    fp_results_wide <- fp_results_wide[, final_fp_cols_exist, drop = FALSE]
  } else { warning("No FP/Linear summary results generated.") }
  
  
  # 4. Combine and Reshape Categorical Summaries (Keep as before)
  cat_summary_long <- do.call(rbind, cat_summary_list)
  cat_selection_wide <- NULL
  if (!is.null(cat_summary_long) && nrow(cat_summary_long) > 0) {
    cat_summary_long$Selected <- ifelse(cat_summary_long$Selected == "+", "+", "-") # Ensure only +/-
    cat_selection_wide <- dcast(cat_summary_long, model ~ term, value.var = "Selected", fill = "-")
    
    # Renaming (Keep as before)
    current_names <- names(cat_selection_wide)
    new_names <- current_names
    new_names <- gsub("^rtsurg", "Surg_", new_names); new_names <- gsub("^racegrp", "Race_", new_names)
    new_names <- gsub("^ER_imp", "ER_", new_names); new_names <- gsub("^PR_imp", "PR_", new_names)
    new_names <- gsub("^married", "Married_", new_names); new_names <- gsub("^innerkwd", "InnerKwd_", new_names)
    new_names <- gsub("^Lateral", "Lateral_", new_names); new_names <- gsub("^grade2_3", "Grade_2or3", new_names); new_names <- gsub("^grade3", "Grade_3", new_names)
    new_names <- gsub("Positive$", "Pos", new_names); new_names <- gsub("bcrt$", "BCRT", new_names)
    new_names <- gsub("mart$", "MART", new_names); new_names <- gsub("ma$", "MA", new_names)
    new_names <- gsub("black$", "Black", new_names); new_names <- gsub("asian$", "Asian", new_names)
    new_names <- gsub("whisp$", "Hisp", new_names)
    new_names <- gsub("([a-zA-Z])(\\d+)$", "\\1_\\2", new_names)
    names(cat_selection_wide) <- new_names
    message("Applied renaming to categorical columns:")
    print(names(cat_selection_wide))
    
  } else { warning("No categorical summary results generated.") }
  
  
  # 5. Reorder rows (Keep as before)
  expected_model_order <- names(model_list)
  if (!is.null(fp_results_wide)) {
    fp_results_wide$model <- as.character(fp_results_wide$model)
    fp_results_wide <- fp_results_wide[match(expected_model_order, fp_results_wide$model), ]
    fp_results_wide$model <- factor(fp_results_wide$model, levels = expected_model_order)
    rownames(fp_results_wide) <- NULL
  }
  if (!is.null(cat_selection_wide)) {
    cat_selection_wide$model <- as.character(cat_selection_wide$model)
    cat_selection_wide <- cat_selection_wide[match(expected_model_order, cat_selection_wide$model), ]
    cat_selection_wide$model <- factor(cat_selection_wide$model, levels = expected_model_order)
    rownames(cat_selection_wide) <- NULL
  }
  
  # 6. Return results (Keep as before)
  return(list(fp_powers = fp_results_wide, cat_selection = cat_selection_wide))
}
# Extract and Save Results
aic_results_mfp2 <- extract_mfp2_results(aic_models_mfp2)
bic_results_mfp2 <- extract_mfp2_results(bic_models_mfp2)

if (!is.null(aic_results_mfp2$fp_powers)) write.csv(aic_results_mfp2$fp_powers, file = "mfp2_aic_fp_powers_summary.csv", row.names = FALSE)
if (!is.null(aic_results_mfp2$cat_selection)) write.csv(aic_results_mfp2$cat_selection, file = "mfp2_aic_categorical_selection_summary.csv", row.names = FALSE)
if (!is.null(bic_results_mfp2$fp_powers)) write.csv(bic_results_mfp2$fp_powers, file = "mfp2_bic_fp_powers_summary.csv", row.names = FALSE)
if (!is.null(bic_results_mfp2$cat_selection)) write.csv(bic_results_mfp2$cat_selection, file = "mfp2_bic_categorical_selection_summary.csv", row.names = FALSE)
message("mfp2 results tables saved.")


