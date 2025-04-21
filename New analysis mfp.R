# --- 1. Load Libraries and Data ---
# install.packages(c("ggplot2", "mfp", "reshape2", "survival", "corrplot")) # Run once if needed
library(ggplot2)
library(mfp)
library(reshape2)
library(survival) # Needed for Surv()
library(corrplot) # Needed for corrplot()

# Set seed for reproducible sampling
set.seed(123)

# Load the main dataset
tryCatch({
  dat <- read.csv("data/dat_ida.csv", header = TRUE)
  message("Successfully loaded data/dat_ida.csv")
}, error = function(e) {
  stop("Error loading data/dat_ida.csv. Please ensure the file path is correct and the file exists. Original error: ", e$message)
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
message("Starting preprocessing...")

# Create grade dummy variables
if ("Grade_imp" %in% names(dat)) {
  dat$grade2_3 <- ifelse(dat$Grade_imp > 1, 1, 0)
  dat$grade3 <- ifelse(dat$Grade_imp == 3, 1, 0)
  message("Created grade dummy variables.")
} else {
  warning("Column 'Grade_imp' not found. Skipping grade variable creation.")
}

# Convert relevant columns to factors
factor_cols <- c("rtsurg", "racegrp", "Lateral", "ER_imp", "PR_imp", "married", "innerkwd")
for (col in factor_cols) {
  if (col %in% names(dat)) {
    dat[[col]] <- factor(dat[[col]])
    message(paste("Converted", col, "to factor."))
  } else {
    warning(paste("Column", col, "not found. Skipping factor conversion."))
  }
}

# Relevel racegrp (ensure 'white' exists as a level)
if ("racegrp" %in% names(dat) && "white" %in% levels(dat$racegrp)) {
  dat$racegrp <- relevel(dat$racegrp, ref = "white")
  message("Releveled 'racegrp' with 'white' as reference.")
} else if ("racegrp" %in% names(dat)) {
  warning("'white' not found in levels of 'racegrp'. Skipping releveling.")
  print(levels(dat$racegrp))
} else {
  warning("Column 'racegrp' not found. Skipping releveling.")
}

message("Preprocessing finished.")

# --- 3. Create Data Subsets (Sampling done ONCE) ---
message("Creating data subsets...")
n_total <- nrow(dat)
sample_indices <- 1:n_total

# Define sample sizes and replications
sample_defs <- list(
  "200K" = list(size = 200000, reps = 2),
  "100K" = list(size = 100000, reps = 2),
  "50K" = list(size = 50000, reps = 2),
  "25K" = list(size = 25000, reps = 3),
  "10K" = list(size = 10000, reps = 3),
  "5K" = list(size = 5000, reps = 3),
  "2K5" = list(size = 2500, reps = 3), # Changed 2.5K to 2K5
  "1K" = list(size = 1000, reps = 3)
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
  for (i in 1:reps) {
    subset_name <- paste0("dat_", name, "_", i)
    sampled_rows <- sample(sample_indices, size, replace = FALSE) # Sample without replacement
    dat_subsets[[subset_name]] <- dat[sampled_rows, ]
    message(paste("Created subset:", subset_name, "with", size, "rows"))
  }
}
# Add the full dataset for M1 model
dat_subsets[["dat_Full"]] <- dat

# Clean up intermediate objects if needed
# rm(sampled_rows, sample_indices, n_total, size, reps, i, name, subset_name, col, factor_cols, required_cols)

# Store subsets in a list
dat_subsets <- list()

for (name in names(sample_defs)) {
  size <- sample_defs[[name]]$size
  reps <- sample_defs[[name]]$reps
  if (size > n_total) {
    warning(paste("Sample size", size, "is larger than total data", n_total, ". Skipping subset", name))
    next
  }
  for (i in 1:reps) {
    # Ensure valid R object names for list elements
    subset_name <- paste0("dat_", gsub("\\.", "K", name), "_", i) # Replaces '.' with 'K'
    subset_name <- gsub("-", "_", subset_name) # Replace any hyphens if used
    sampled_rows <- sample(sample_indices, size, replace = FALSE) # Sample without replacement
    dat_subsets[[subset_name]] <- dat[sampled_rows, ]
    message(paste("Created subset:", subset_name, "with", size, "rows"))
  }
}

# --- ADD THIS SECTION TO SAVE THE LIST ---
message("Saving the list of data subsets...")
tryCatch({
  # Choose a filename
  saveRDS(dat_subsets, file = "dat_subsets_list.rds")
  message("Data subsets list successfully saved to dat_subsets_list.rds")
}, error = function(e) {
  warning("Error saving data subsets list: ", e$message)
  # Consider adding more robust error handling if needed
})

# --- 4. Define Model Formula ---
# Using the structure from the user's "working" example (assuming variables are factors)
# Note: factor() is kept for Lateral as per original code, check if Lateral is numeric 0/1
# If Lateral is already a factor after preprocessing, remove factor()
# If ER_imp, PR_imp etc are NOT factors after preprocessing, add factor() back
# For simplicity, we assume preprocessing step correctly made them factors
mfp_formula <- Surv(futime, Death) ~
  fp(AgeDgc, df = 4) +
  fp(tSizeMm, df = 4) +
  fp(EodNPos, df = 4) +
  YearDgc +          # Assumed linear unless selected otherwise by fp() implicitly
  grade2_3 + grade3 + # Already numeric 0/1
  rtsurg +           # Assumed factor from preprocessing
  racegrp +          # Assumed factor with correct ref level
  Lateral +          # Assumed factor from preprocessing
  ER_imp +           # Assumed factor
  PR_imp +           # Assumed factor
  married +          # Assumed factor
  innerkwd           # Assumed factor

# --- 5. Run MFP Models (Using Loops) ---
message("Fitting MFP models...")

# Define model names corresponding to subsets (adjust if subset names change)
model_names_map <- list(
  "M1" = "dat_Full",
  "M200K1" = "dat_200K_1", "M200K2" = "dat_200K_2",
  "M100K1" = "dat_100K_1", "M100K2" = "dat_100K_2",
  "M50K1" = "dat_50K_1", "M50K2" = "dat_50K_2",
  "M25K1" = "dat_25K_1", "M25K2" = "dat_25K_2", "M25K3" = "dat_25K_3",
  "M10K1" = "dat_10K_1", "M10K2" = "dat_10K_2", "M10K3" = "dat_10K_3",
  "M5K1" = "dat_5K_1", "M5K2" = "dat_5K_2", "M5K3" = "dat_5K_3",
  "M2K5_1" = "dat_2K5_1", "M2K5_2" = "dat_2K5_2", "M2K5_3" = "dat_2K5_3", # Changed names
  "M1K1" = "dat_1K_1", "M1K2" = "dat_1K_2", "M1K3" = "dat_1K_3"
)

# Store models in lists
aic_models <- list()
bic_models <- list()

# --- 5a. AIC Models (select = 0.157) ---
message("Fitting AIC-based models (select=0.157)...")
aic_select_level <- 0.157



for (model_name in names(model_names_map)) {
  data_name <- model_names_map[[model_name]]
  if (!data_name %in% names(dat_subsets)) {
    warning(paste("Data subset", data_name, "not found for model", model_name, ". Skipping."))
    next
  }
  current_data <- dat_subsets[[data_name]]
  message(paste("Fitting", model_name, "using data", data_name, "(", nrow(current_data), "rows ) with AIC criterion..."))
  
  # Use tryCatch to handle potential errors during fitting
  aic_models[[model_name]] <- tryCatch({
    mfp(formula = mfp_formula,
        data = current_data,
        family = cox,
        select = aic_select_level)
  }, error = function(e) {
    warning(paste("Error fitting", model_name, ":", e$message))
    return(NULL) # Return NULL if model fails
  })
  if (!is.null(aic_models[[model_name]])) {
    message(paste("Finished fitting", model_name))
  }
}


# --- 5b. BIC Models (select based on events) ---
message("Fitting BIC-based models (select=pchisq(log(events),1))...")

for (model_name in names(model_names_map)) {
  model_name_bic <- paste0(model_name, "b") # e.g., M1b, M200K1b
  data_name <- model_names_map[[model_name]]
  if (!data_name %in% names(dat_subsets)) {
    # Warning already given in AIC loop
    next
  }
  current_data <- dat_subsets[[data_name]]
  
  # Calculate BIC-based selection level (slc)
  events <- sum(current_data$Death, na.rm = TRUE)
  if (events <= 0) {
    warning(paste("No events found in", data_name, "for model", model_name_bic, ". Skipping."))
    next
  }
  slc <- round(pchisq(log(events), 1, lower.tail = FALSE), 4)
  message(paste("Fitting", model_name_bic, "using data", data_name, "(", nrow(current_data), "rows, ", events, " events) with BIC criterion (select=", slc, ")..."))
  
  # Use tryCatch
  bic_models[[model_name_bic]] <- tryCatch({
    mfp(formula = mfp_formula,
        data = current_data,
        family = cox,
        select = slc)
  }, error = function(e) {
    warning(paste("Error fitting", model_name_bic, ":", e$message))
    return(NULL)
  })
  if (!is.null(bic_models[[model_name_bic]])) {
    message(paste("Finished fitting", model_name_bic))
  }
}

message("Finished all model fitting.")

# --- 6. Create Results Tables ---
message("Generating results tables...")

extract_mfp_results <- function(model_list, model_type_suffix = "") {
  
  fp_vars <- c("AgeDgc", "tSizeMm", "EodNPos", "YearDgc") # Continuous variables potentially having FP terms
  
  fp_powers_list <- list()
  cat_selection_list <- list()
  
  valid_model_names <- names(model_list)[!sapply(model_list, is.null)]
  if (length(valid_model_names) == 0) {
    warning("No valid models found in the provided list.")
    return(list(fp_powers = NULL, cat_selection = NULL))
  }
  
  # Dynamically get categorical variable names from the first valid model's fptable
  first_valid_model <- model_list[[valid_model_names[1]]]
  cat_vars <- character(0) # Initialize empty vector
  
  if (!is.null(first_valid_model$fptable)) {
    fptable_rownames <- rownames(first_valid_model$fptable)
    # Identify potential categorical vars by excluding fp_vars and FP interaction terms (like .1, .2)
    potential_cat_vars <- fptable_rownames[!fptable_rownames %in% fp_vars & !grepl("\\.[1-9]$", fptable_rownames)]
    cat_vars <- potential_cat_vars
    message("Dynamically identified potential categorical terms: ", paste(cat_vars, collapse=", "))
  } else {
    warning("Could not access fptable of first valid model to get categorical variables dynamically. Cannot process categorical selection.")
  }
  
  
  for (model_name in valid_model_names) {
    model <- model_list[[model_name]]
    if (is.null(model$fptable) || is.null(model$coefficients)) {
      warning(paste("Skipping results extraction for", model_name, "as fptable or coefficients are missing."))
      next
    }
    fptable <- model$fptable
    coefs <- model$coefficients
    coef_names <- names(coefs)
    
    # Extract FP powers (same as before)
    fp_idx <- rownames(fptable) %in% fp_vars
    if(any(fp_idx)) {
      fp_nms <- rownames(fptable)[fp_idx]
      # Ensure power columns exist, provide NA otherwise
      power1_col <- if("power1" %in% colnames(fptable)) fptable[fp_idx, "power1"] else rep(NA, length(fp_nms))
      power2_col <- if("power2" %in% colnames(fptable)) fptable[fp_idx, "power2"] else rep(NA, length(fp_nms))
      fp_powers_list[[model_name]] <- data.frame(model = model_name, vars = fp_nms,
                                                 power1 = power1_col, power2 = power2_col)
    } else {
      message(paste("No FP variables found in fptable for model", model_name))
    }
    
    # --- Revised Categorical Selection ---
    # Check if dynamically identified cat_vars exist in *this* model's fptable
    current_cat_vars_idx <- rownames(fptable) %in% cat_vars
    if(any(current_cat_vars_idx)) {
      current_cat_nms <- rownames(fptable)[current_cat_vars_idx]
      # Determine selection by checking presence in coefficients
      # A term from fptable might correspond to one *or more* coefficients (e.g., multi-level factors)
      # Simple check: is the exact name from fptable present in coefficients?
      # More robust: handle factor levels e.g. `rtsurg` in fptable becomes `rtsurgbcrt`, `rtsurgma` etc. in coefs
      # We'll use the exact names found in fptable (which often already include the level)
      cat_sel <- ifelse(current_cat_nms %in% coef_names, "+", "-")
      
      cat_selection_list[[model_name]] <- data.frame(model = model_name, vars = current_cat_nms, selection = cat_sel)
    } else {
      # This case should be less common if cat_vars was derived correctly
      message(paste("None of the dynamically identified categorical variables found in fptable for model", model_name))
    }
  }
  
  # Combine results into data frames
  fp_powers_df <- do.call(rbind, fp_powers_list)
  cat_selection_df <- do.call(rbind, cat_selection_list)
  
  # Reshape FP powers (same as before)
  fp_results_wide <- NULL
  if (!is.null(fp_powers_df) && nrow(fp_powers_df) > 0) {
    # Use fill=NA for models where a variable might be missing
    fp_p1_wide <- dcast(fp_powers_df, model ~ vars, value.var = "power1", fill = NA)
    fp_p2_wide <- dcast(fp_powers_df, model ~ vars, value.var = "power2", fill = NA)
    
    all_fp_vars <- unique(fp_powers_df$vars)
    fp_results_wide <- data.frame(model = fp_p1_wide$model)
    for(v in all_fp_vars) {
      # Check if columns exist before accessing
      if (v %in% names(fp_p1_wide)) fp_results_wide[[paste0(v, "_P1")]] <- fp_p1_wide[[v]]
      if (v %in% names(fp_p2_wide)) fp_results_wide[[paste0(v, "_P2")]] <- fp_p2_wide[[v]]
    }
    # Specific handling for YearDgc if often linear
    if("YearDgc_P1" %in% names(fp_results_wide)) names(fp_results_wide)[names(fp_results_wide)=="YearDgc_P1"] <- "YearDgc_SelectedP1" # Indicate it's selected P1=1 often
    if("YearDgc_P2" %in% names(fp_results_wide)) fp_results_wide$YearDgc_P2 <- NULL # Remove P2 if usually NA
  } else {
    warning("No FP power results to reshape.")
  }
  
  
  # Reshape categorical selections
  cat_results_wide <- NULL
  if (!is.null(cat_selection_df) && nrow(cat_selection_df) > 0) {
    # Use fill="-" for models where a variable might be missing (implies not selected)
    cat_results_wide <- dcast(cat_selection_df, model ~ vars, value.var = "selection", fill = "-")
    
    # --- Refined Renaming based on observed names ---
    current_names <- names(cat_results_wide)
    new_names <- current_names
    
    # Apply specific replacements based on the "Dynamically identified" list:
    # grade3, married1, rtsurgbcrt, rtsurgma, rtsurgmart, racegrpasian, racegrpblack,
    # racegrpwhisp, grade2_3, PR_impPositive, ER_impPositive, innerkwd1, Lateral2
    new_names <- gsub("^grade3$", "Grade_3", new_names)
    new_names <- gsub("^married1$", "Married", new_names) # Assuming 1 is the non-ref level shown
    new_names <- gsub("^rtsurgbcrt$", "Surg_BCRT", new_names)
    new_names <- gsub("^rtsurgma$", "Surg_MA", new_names)
    new_names <- gsub("^rtsurgmart$", "Surg_MART", new_names)
    new_names <- gsub("^racegrpasian$", "Race_Asian", new_names)
    new_names <- gsub("^racegrpblack$", "Race_Black", new_names)
    new_names <- gsub("^racegrpwhisp$", "Race_Hispanic", new_names) # Assuming 'whisp' means White Hispanic
    new_names <- gsub("^grade2_3$", "Grade_2or3", new_names)
    new_names <- gsub("^PR_impPositive$", "PR_Positive", new_names)
    new_names <- gsub("^ER_impPositive$", "ER_Positive", new_names)
    new_names <- gsub("^innerkwd1$", "InnerKwd", new_names) # Assuming 1 is the non-ref level shown
    new_names <- gsub("^Lateral2$", "Lateral_2", new_names) # Assuming 2 is the non-ref level shown
    
    # Check if any names failed to rename (for debugging)
    if(any(new_names == current_names & new_names != "model")){
      failed_rename <- current_names[new_names == current_names & new_names != "model"]
      warning("Could not rename the following columns: ", paste(failed_rename, collapse=", "))
    }
    
    names(cat_results_wide) <- new_names
    message("Applied refined renaming to categorical columns. Please verify names:")
    print(names(cat_results_wide))
    
  } else {
    warning("No categorical selection results to reshape.")
  }
  
  # Reorder rows based on original model order (approximate, adapt as needed)
  # Create the full list of expected model names including the suffix
  expected_model_order <- names(model_list) # Use names from the input list directly
  
  if (!is.null(fp_results_wide)) {
    # Match expects exact names; ensure model names in fp_results_wide match expected_model_order
    fp_results_wide <- fp_results_wide[match(expected_model_order, fp_results_wide$model), ]
    rownames(fp_results_wide) <- NULL # Reset row numbers
  }
  if (!is.null(cat_results_wide)) {
    cat_results_wide <- cat_results_wide[match(expected_model_order, cat_results_wide$model), ]
    rownames(cat_results_wide) <- NULL # Reset row numbers
  }
  
  
  return(list(fp_powers = fp_results_wide, cat_selection = cat_results_wide))
}



# Extract for AIC models
aic_results <- extract_mfp_results(aic_models)
aic_fp_powers <- aic_results$fp_powers
aic_cat_selection <- aic_results$cat_selection


# Extract for BIC models
bic_results <- extract_mfp_results(bic_models, model_type_suffix = "b") # Suffix 'b' might not be needed if model names in list already have it
bic_fp_powers <- bic_results$fp_powers
bic_cat_selection <- bic_results$cat_selection



# Rename columns more descriptively if needed (example based on original j1/j2 structure)
# This requires inspecting the actual column names produced by dcast above
# Example reordering (adjust based on actual names and desired order):
# if (!is.null(aic_fp_powers)) {
#   # aic_fp_powers <- aic_fp_powers[, c("model", "AgeDgc_P1", "AgeDgc_P2", ...)]
# }
# if (!is.null(aic_cat_selection)) {
#    # aic_cat_selection <- aic_cat_selection[, c("model", "Grade2or3", "Grade3", ...)]
#}

# Save results with descriptive names
if (!is.null(aic_fp_powers)) write.csv(aic_fp_powers, file = "mfp_aic_fp_powers_summary.csv", row.names = FALSE)
if (!is.null(aic_cat_selection)) write.csv(aic_cat_selection, file = "mfp_aic_categorical_selection_summary.csv", row.names = FALSE)
if (!is.null(bic_fp_powers)) write.csv(bic_fp_powers, file = "mfp_bic_fp_powers_summary.csv", row.names = FALSE)
if (!is.null(bic_cat_selection)) write.csv(bic_cat_selection, file = "mfp_bic_categorical_selection_summary.csv", row.names = FALSE)

message("Results tables saved.")

# --- 7. Correlation Plot ---
message("Generating correlation plot...")
# Select columns for correlation (ensure they exist and are numeric/convertible)
cor_cols_options <- c("AgeDgc", "tSizeMm", "EodNPos", "YearDgc", "Grade_imp",
                      "rtsurg", "racegrp", "Lateral", "ER_imp", "PR_imp",
                      "married", "innerkwd")
cor_cols_exist <- intersect(cor_cols_options, names(dat))

cor_dat_list <- list()
for (col in cor_cols_exist) {
  if (is.numeric(dat[[col]])) {
    cor_dat_list[[col]] <- dat[[col]]
  } else if (is.factor(dat[[col]])) {
    # Convert factor to numeric - use with caution, depends on encoding
    cor_dat_list[[col]] <- as.numeric(dat[[col]])
    message(paste("Converting factor", col, "to numeric for correlation."))
  } else {
    warning(paste("Column", col, "is not numeric or factor. Skipping for correlation."))
  }
}
cor_dat <- as.data.frame(cor_dat_list)
cor_dat <- na.omit(cor_dat) # Handle missing values

if (ncol(cor_dat) > 1) {
  cor_table <- round(cor(cor_dat, method = "spearman"), 2)
  # Plotting - consider saving to a file
  png("correlation_plot_spearman.png", width=800, height=800)
  corrplot(cor_table, method = "number", type="lower", tl.col="black", tl.srt=45)
  dev.off()
  message("Correlation plot saved to correlation_plot_spearman.png")
  # print(cor_table) # Optionally print the table
} else {
  warning("Could not compute correlation matrix (less than 2 valid columns).")
}


# --- 8. Plot Functional Forms (Example for AgeDgc, requires careful adaptation) ---
message("Generating functional form plots (example for AgeDgc)...")

# Helper function to get FP prediction for one variable
# NOTE: This is complex because predict.mfp doesn't easily isolate one variable's effect.
# We need to manually parse the transformation and apply coefficients.
# This is prone to errors and specific to how mfp stores results.
# Simplified version: Extracts formula terms and coefficients by name.

get_fp_prediction <- function(model, var_name, var_values) {
  if (is.null(model) || is.null(model$coefficients) || is.null(model$trafo)) {
    warning(paste("Model invalid or missing components for", var_name))
    return(rep(NA, length(var_values)))
  }
  
  # Find transformation string for the variable
  trafo_string <- model$trafo[var_name, 1]
  if (is.na(trafo_string) || nchar(trafo_string) == 0) {
    message(paste(var_name, "not found or not transformed in model", deparse(substitute(model))))
    # Check if it's included linearly (check coefficients directly)
    if(var_name %in% names(model$coefficients)) {
      coef_val <- model$coefficients[var_name]
      return(coef_val * var_values)
    } else {
      return(rep(0, length(var_values))) # Assume zero contribution if not found
    }
  }
  
  # Find corresponding coefficients (e.g., AgeDgc.1, AgeDgc.2)
  coef_names <- names(model$coefficients)
  var_coef_indices <- grep(paste0("^", var_name, "\\.[1-9]"), coef_names) # Matches Var.1, Var.2 etc.
  var_coefs <- model$coefficients[var_coef_indices]
  
  if (length(var_coefs) == 0) {
    message(paste("Coefficients for", var_name, "not found in model", deparse(substitute(model))))
    # Check for linear term again
    if(var_name %in% names(model$coefficients)) {
      coef_val <- model$coefficients[var_name]
      return(coef_val * var_values)
    } else {
      return(rep(0, length(var_values)))
    }
  }
  
  # Create a temporary environment for evaluation
  eval_env <- new.env()
  assign(var_name, var_values, envir = eval_env)
  
  # Evaluate the transformation string parts
  # This requires parsing the string like "I((AgeDgc/100)^0.5)+I((AgeDgc/100)^1)"
  # Simplified: assumes '+' separates terms and 'I()' wraps them.
  # Robust parsing is needed for complex cases.
  terms <- strsplit(trafo_string, "\\+")[[1]]
  prediction <- rep(0, length(var_values))
  
  if (length(terms) != length(var_coefs)) {
    warning(paste("Mismatch between terms in trafo string and coefficients for", var_name, "in model", deparse(substitute(model))))
    # Attempt to match based on names if possible (e.g. if some coefs are zero)
    # Fallback: return NA or zero
    return(rep(NA, length(var_values)))
  }
  
  for (i in 1:length(terms)) {
    term_expr <- try(parse(text = terms[i]), silent = TRUE)
    if (!inherits(term_expr, "try-error")) {
      term_value <- try(eval(term_expr, envir = eval_env), silent = TRUE)
      if (!inherits(term_value, "try-error")) {
        prediction <- prediction + var_coefs[i] * term_value
      } else {
        warning(paste("Could not evaluate term:", terms[i]))
        return(rep(NA, length(var_values)))
      }
    } else {
      warning(paste("Could not parse term:", terms[i]))
      return(rep(NA, length(var_values)))
    }
  }
  return(prediction)
}


# Example Plot for AgeDgc (AIC models)
png("functional_forms_AgeDgc_AIC.png", width = 1000, height = 600)
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1) # Reset layout
age_range <- seq(min(dat$AgeDgc, na.rm=TRUE), max(dat$AgeDgc, na.rm=TRUE), length.out = 100)
reference_age <- median(dat$AgeDgc, na.rm=TRUE) # Or choose a specific age like 18 if needed

# Get prediction for reference model (M1)
ref_model_name <- "M1"
ref_pred_age <- NULL
if (ref_model_name %in% names(aic_models) && !is.null(aic_models[[ref_model_name]])) {
  ref_pred_age <- get_fp_prediction(aic_models[[ref_model_name]], "AgeDgc", age_range)
  ref_pred_ref_val <- get_fp_prediction(aic_models[[ref_model_name]], "AgeDgc", reference_age)
  if (!is.null(ref_pred_age) && !is.null(ref_pred_ref_val) && !any(is.na(ref_pred_age)) && !any(is.na(ref_pred_ref_val))) {
    plot(age_range, ref_pred_age - ref_pred_ref_val, type = 'l', lwd = 4, col = "red",
         ylim = c(-2, 2), # Adjust ylim based on observed range
         main = "Functional Form for Age (AIC Models, Centered)",
         xlab = "Age (years)", ylab = "Log Hazard Ratio (Centered)", axes = TRUE)
  } else {
    plot(1, type="n", xlim=range(age_range), ylim=c(-2,2), main="Functional Form Error", xlab="Age", ylab="Log HR")
    text(mean(age_range), 0, "Could not generate reference plot for AgeDgc")
    ref_pred_age <- NULL # Ensure it's null if plot failed
  }
} else {
  plot(1, type="n", xlim=range(age_range), ylim=c(-2,2), main="Functional Form Error", xlab="Age", ylab="Log HR")
  text(mean(age_range), 0, "Reference Model M1 (AIC) not found or invalid.")
}


# Overlay predictions for other models
for (model_name in names(aic_models)) {
  if (model_name == ref_model_name) next # Skip reference model
  model <- aic_models[[model_name]]
  if (is.null(model)) next
  
  pred_age <- get_fp_prediction(model, "AgeDgc", age_range)
  pred_ref_val <- get_fp_prediction(model, "AgeDgc", reference_age)
  
  if (!is.null(pred_age) && !is.null(pred_ref_val) && !any(is.na(pred_age)) && !any(is.na(pred_ref_val))) {
    lines(age_range, pred_age - pred_ref_val, col = "grey", lwd = 1)
  } else {
    message(paste("Could not get valid AgeDgc prediction for model", model_name))
  }
}

# Redraw reference line on top
if (!is.null(ref_pred_age) && !is.null(ref_pred_ref_val)) {
  lines(age_range, ref_pred_age - ref_pred_ref_val, lwd = 4, col = "red")
}
legend("topleft", c("Reference (M1)", "Sample Models"), col = c("red", "grey"), lwd = c(4, 1), bty = "n")
dev.off()
message("Functional form plot for AgeDgc (AIC) saved.")

# Repeat for BIC models, and other variables (tSizeMm, EodNPos) if desired,
# remembering to adapt the variable name and potentially the x/y ranges.

# --- 9. Compare Predictions (Using predict()) ---
message("Comparing model predictions...")

# Create a new validation dataset (if needed, or use an existing subset)
# Using a subset already created, e.g., dat_10K_1
validation_data_name <- "dat_10K_1"
if (validation_data_name %in% names(dat_subsets)) {
  new_dat <- dat_subsets[[validation_data_name]]
  message(paste("Using subset", validation_data_name, "for prediction comparison."))
} else {
  warning("Validation subset not found. Skipping prediction comparison.")
  new_dat <- NULL
}


if (!is.null(new_dat)) {
  # Store predictions
  predictions_lp <- list() # lp = linear predictor
  
  # Get prediction from reference model (e.g., M1b - BIC full data model)
  ref_model_bic_name <- "M1b"
  if (ref_model_bic_name %in% names(bic_models) && !is.null(bic_models[[ref_model_bic_name]])) {
    predictions_lp[[ref_model_bic_name]] <- tryCatch({
      predict(bic_models[[ref_model_bic_name]], newdata = new_dat, type = "lp")
    }, error = function(e) {
      warning(paste("Could not predict for reference model", ref_model_bic_name, ":", e$message))
      return(NULL)
    })
  } else {
    warning(paste("Reference BIC model", ref_model_bic_name, "not found or invalid for prediction."))
  }
  
  # Get predictions from selected sample models (e.g., M10K1b, M5K2b)
  sample_model_names_pred <- c("M10K1b", "M5K2b") # Add more as needed
  for (model_name in sample_model_names_pred) {
    if (model_name %in% names(bic_models) && !is.null(bic_models[[model_name]])) {
      predictions_lp[[model_name]] <- tryCatch({
        predict(bic_models[[model_name]], newdata = new_dat, type = "lp")
      }, error = function(e) {
        warning(paste("Could not predict for sample model", model_name, ":", e$message))
        return(NULL)
      })
    } else {
      warning(paste("Sample model", model_name, "not found or invalid for prediction."))
    }
  }
  
  # Create comparison plots (Scatter and Bland-Altman)
  ref_pred <- predictions_lp[[ref_model_bic_name]]
  
  if(!is.null(ref_pred)) {
    for (model_name in sample_model_names_pred) {
      sample_pred <- predictions_lp[[model_name]]
      if (!is.null(sample_pred)) {
        plot_title_scatter <- paste("Prediction Comparison:", ref_model_bic_name, "vs", model_name)
        plot_title_ba <- paste("Bland-Altman:", ref_model_bic_name, "vs", model_name)
        file_suffix <- paste0(ref_model_bic_name, "_vs_", model_name)
        
        # Scatter Plot
        png(paste0("prediction_scatter_", file_suffix, ".png"))
        plot(ref_pred, sample_pred,
             xlab = paste("Reference Model (", ref_model_bic_name, ") Linear Predictor"),
             ylab = paste("Sample Model (", model_name, ") Linear Predictor"),
             main = plot_title_scatter, pch = 19, col = rgb(0,0,0,0.3))
        abline(0, 1, col = "red", lwd = 2) # Line of identity
        dev.off()
        
        # Bland-Altman Plot
        model_mean <- (ref_pred + sample_pred) / 2
        model_diff <- sample_pred - ref_pred
        ba_data <- data.frame(mean = model_mean, diff = model_diff)
        
        ba_plot <- ggplot(ba_data, aes(x = mean, y = diff)) +
          geom_point(alpha = 0.5) +
          geom_hline(yintercept = mean(model_diff, na.rm=TRUE), color = "blue", linetype = "dashed") +
          geom_hline(yintercept = mean(model_diff, na.rm=TRUE) + 1.96 * sd(model_diff, na.rm=TRUE), color = "red", linetype = "dashed") +
          geom_hline(yintercept = mean(model_diff, na.rm=TRUE) - 1.96 * sd(model_diff, na.rm=TRUE), color = "red", linetype = "dashed") +
          labs(title = plot_title_ba,
               x = "Mean of Linear Predictors",
               y = "Difference of Linear Predictors (Sample - Reference)") +
          theme_bw()
        
        ggsave(paste0("prediction_bland_altman_", file_suffix, ".png"), plot = ba_plot)
        
        message(paste("Prediction comparison plots saved for", file_suffix))
      }
    }
  } else {
    message("Reference prediction missing, cannot create comparison plots.")
  }
}

message("Script finished.")