
# ============================================================
# MCD Pipeline: Metadata-driven comparison definition + export
# ============================================================

# Suggested packages:
# install.packages(c("openxlsx"))

`%||%` <- function(x, y) if (is.null(x)) y else x


# ------------------------------------------------------------
# Helper: flatten named list for export
# ------------------------------------------------------------
flatten_settings <- function(x, prefix = NULL) {
  out <- list()
  
  recurse <- function(obj, nm = NULL) {
    if (is.list(obj) && !is.data.frame(obj)) {
      nms <- names(obj)
      if (is.null(nms)) nms <- seq_along(obj)
      for (i in seq_along(obj)) {
        next_nm <- if (is.null(nm)) as.character(nms[i]) else paste(nm, nms[i], sep = ".")
        recurse(obj[[i]], next_nm)
      }
    } else {
      out[[nm]] <<- paste(obj, collapse = "; ")
    }
  }
  
  recurse(x, prefix)
  
  data.frame(
    setting = names(out),
    value = unlist(out, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# Helper: interactive choice
# ------------------------------------------------------------
choose_comparison_mode <- function() {
  
  
  show_header("Select MCD analysis mode")
  cat("1. global\n")
  cat("2. treatment_vs_control\n")
  cat("3. pairwise_by_treatment\n")
  cat("4. custom\n\n")
  
  repeat{
  
    if (interactive()) {
      choice <- trimws(readline("Enter choice [1-4]: "))
    } else {
      cat("Enter choice [1-4]: ")
      choice <- trimws(readLines("stdin", n=1))
    }
    
    
    mode <- switch(
      choice,
      "1" = "global",
      "2" = "treatment_vs_control",
      "3" = "pairwise_by_treatment",
      "4" = "custom"
    )
    
    if (!is.null(mode)) return(mode)
      
    cat("Invalid choice. Please enter 1, 2, 3, or 4.")
  }

}

# ------------------------------------------------------------
# Helper: prompt for grouping variables
# ------------------------------------------------------------
choose_grouping_vars <- function(metadata) {
  cat("\nAvailable metadata columns:\n")
  cat(paste0(" - ", names(metadata), collapse = "\n"), "\n\n")
  
  default = c("cell_line", "timepoint")  
  prompt <- paste0(
    "Enter grouping variables separated by commas",
     paste0(" [default: ", paste(default, collapse = ", "), "]"),
    ": ", "\n If default, Press ENTER \n"
  )
  
  repeat{
    
    if (interactive()) {
      x <- trimws(readline(prompt))
    } else {
      cat(prompt)
      x <- trimws(readLines("stdin", n=1))
    }
    
    if (nchar(trimws(x)) == 0) {
      vars <- default
    } else {
      vars <- trimws(strsplit(x, ",")[[1]])
    }
    
    vars <- vars[vars %in% names(metadata)]
    
    if (!length(vars) == 0) return(vars)
    
    cat("No valid grouping variables selected. \n")
    cat("Select columns available in metadata as grouping variable. \n")
  }
}

# choose treatment var
choose_treatment_vars <- function(metadata) {
  
  cat("\nAvailable metadata columns:\n")
  cat(paste0(" - ", names(metadata), collapse = "\n"), "\n\n")
  
  default = 'treatment'
  prompt <- paste0(
    "Choose treatment column.",
    paste0(" [default: ", default, "]"),
    ": ", "\nIf default, Press ENTER: \n")
  
  repeat{
    
    if (interactive()) {
      x <- trimws(readline(prompt))
    } else {
      cat(prompt)
      x <- trimws(readLines("stdin", n=1))
    }
    
    if (nchar(trimws(x)) == 0) {
      vars <- default
    } else {
      vars <- trimws(strsplit(x, ",")[[1]])
    }
    
    vars <- vars[vars %in% names(metadata)]
    
    if (!length(vars) == 0) return(vars)
    
    cat("No valid column selected. \n")
    cat("Select columns available in metadata as treatment variable. \n")
  }
}
# ------------------------------------------------------------
# Helper: prompt for control label
# ------------------------------------------------------------
choose_control_label <- function(metadata, treatment_col) {
  
  if (!treatment_col %in% names(metadata)) {
    stop(sprintf("Metadata does not contain '%s' column.", treatment_col))
  }
  
  vals <- sort(unique(as.character(metadata[[treatment_col]])))
  cat("\nAvailable treatment labels:\n")
  cat(paste0(" - ", vals, collapse = "\n"), "\n\n")
  
  repeat{
    
    if (interactive()) {
      ctl <- trimws(readline("Enter control label (e.g. DMSO): "))
    } else {
      cat("Enter control label (e.g. DMSO): ")
      ctl <- trimws(readLines("stdin", n=1))
    }
    
    
    if (ctl %in% vals) return(ctl)
    
    cat(sprintf("Control label '%s' not found in metadata$%s.", 
                ctl, treatment_col))
  }
}


# ------------------------------------------------------------
# Helper: validate final_input object
# ------------------------------------------------------------
validate_final_input <- function(final_input) {
  req <- c("input_type", "score_matrix", "metadata")
  miss <- setdiff(req, names(final_input))
  if (length(miss) > 0) {
    stop("final_input is missing required fields: ", paste(miss, collapse = ", "))
  }
  
  if (!is.matrix(final_input$score_matrix) && !is.data.frame(final_input$score_matrix)) {
    stop("final_input$score_matrix must be a matrix or data.frame.")
  }
  
  score_matrix <- as.matrix(final_input$score_matrix)
  
  if (!is.numeric(score_matrix)) {
    suppressWarnings(storage.mode(score_matrix) <- "numeric")
  }
  
  if (!all(colnames(score_matrix) %in% final_input$metadata$sample)) {
    stop("All score_matrix column names must be present in metadata$sample.")
  }
  
  if (!all(final_input$metadata$sample %in% colnames(score_matrix))) {
    warning("Some metadata samples are not present in score_matrix columns.")
  }
  
  invisible(TRUE)
}

# ------------------------------------------------------------
# Define comparisons
# ------------------------------------------------------------

define_comparisons <- function(metadata,
                               mode,
                               grouping_vars,
                               treatment_col,
                               sample_col = "sample",
                               control_label) {
  if (!sample_col %in% names(metadata)) {
    stop(sprintf("Metadata must contain '%s'.", sample_col))
  }
  
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  
  if (mode == "global") {
    cmp <- data.frame(
      comparison_id = "CMP001",
      comparison_name = "global_all_samples",
      mode = "global",
      group_id = "all_samples",
      group_test = "all",
      group_ref = NA_character_,
      n_test = nrow(metadata),
      n_ref = NA_integer_,
      sample_test = paste(metadata[[sample_col]], collapse = ";"),
      sample_ref = NA_character_,
      valid = TRUE,
      notes = "",
      stringsAsFactors = FALSE
    )
    
    return(cmp)
  }
  
  if (is.null(grouping_vars)) {
    stop("grouping_vars must be supplied for non-global modes.")
  }
  
  missing_group_vars <- setdiff(grouping_vars, names(metadata))
  if (length(missing_group_vars) > 0) {
    stop("Missing grouping variables in metadata: ", paste(missing_group_vars, collapse = ", "))
  }
  
  if (!treatment_col %in% names(metadata)) {
    stop(sprintf("Metadata must contain '%s' for mode '%s'.", treatment_col, mode))
  }
  
  if (length(grouping_vars) == 0) {
    metadata$.__group_id__. <- rep("all_samples", nrow(metadata))
  } else {
    metadata$.__group_id__. <- apply(
      metadata[, grouping_vars, drop = FALSE], 1, 
      function(z) paste(z, collapse = " | ")
      )
  }
  
  res <- list()
  cmp_counter <- 1L
  
  if (mode == "treatment_vs_control") {
    if (is.null(control_label)) {
      stop("control_label must be supplied for treatment_vs_control mode.")
    }
    
    split_groups <- split(metadata, metadata$.__group_id__.)
    
    for (gid in names(split_groups)) {
      gdf <- split_groups[[gid]]
      trts <- sort(unique(as.character(gdf[[treatment_col]])))
      tests <- setdiff(trts, control_label)
      
      for (tx in tests) {
        ref_df <- gdf[gdf[[treatment_col]] == control_label, , drop = FALSE]
        test_df <- gdf[gdf[[treatment_col]] == tx, , drop = FALSE]
        
        valid <- TRUE
        notes <- character()
        
        if (nrow(ref_df) == 0) {
          valid <- FALSE
          notes <- c(notes, "No control samples found")
        }
        if (nrow(test_df) == 0) {
          valid <- FALSE
          notes <- c(notes, "No treatment samples found")
        }
        
        cmp_row <- data.frame(
          comparison_id = sprintf("CMP%03d", cmp_counter),
          comparison_name = paste(c(gsub(" \\| ", "_", gid), paste0(tx, "_vs_", control_label)), collapse = "_"),
          mode = mode,
          group_id = gid,
          group_test = tx,
          group_ref = control_label,
          n_test = nrow(test_df),
          n_ref = nrow(ref_df),
          sample_test = paste(test_df[[sample_col]], collapse = ";"),
          sample_ref = paste(ref_df[[sample_col]], collapse = ";"),
          valid = valid,
          notes = paste(notes, collapse = "; "),
          stringsAsFactors = FALSE
        )
        
        for (v in grouping_vars) {
          cmp_row[[v]] <- unique(gdf[[v]])[1]
        }
        
        res[[length(res) + 1L]] <- cmp_row
        cmp_counter <- cmp_counter + 1L
      }
    }
  }
  
  if (mode == "pairwise_by_treatment") {
    split_groups <- split(metadata, metadata$.__group_id__.)
    
    for (gid in names(split_groups)) {
      gdf <- split_groups[[gid]]
      trts <- sort(unique(as.character(gdf[[treatment_col]])))
      
      if (length(trts) < 2) next
      
      pairs <- utils::combn(trts, 2, simplify = FALSE)
      
      for (pp in pairs) {
        ref_lab <- pp[1]
        test_lab <- pp[2]
        
        ref_df <- gdf[gdf[[treatment_col]] == ref_lab, , drop = FALSE]
        test_df <- gdf[gdf[[treatment_col]] == test_lab, , drop = FALSE]
        
        valid <- TRUE
        notes <- character()
        
        if (nrow(ref_df) == 0 || nrow(test_df) == 0) {
          valid <- FALSE
          notes <- c(notes, "One or both groups missing")
        }
        
        cmp_row <- data.frame(
          comparison_id = sprintf("CMP%03d", cmp_counter),
          comparison_name = paste(c(gsub(" \\| ", "_", gid), paste0(test_lab, "_vs_", ref_lab)), collapse = "_"),
          mode = mode,
          group_id = gid,
          group_test = test_lab,
          group_ref = ref_lab,
          n_test = nrow(test_df),
          n_ref = nrow(ref_df),
          sample_test = paste(test_df[[sample_col]], collapse = ";"),
          sample_ref = paste(ref_df[[sample_col]], collapse = ";"),
          valid = valid,
          notes = paste(notes, collapse = "; "),
          stringsAsFactors = FALSE
        )
        
        for (v in grouping_vars) {
          cmp_row[[v]] <- unique(gdf[[v]])[1]
        }
        
        res[[length(res) + 1L]] <- cmp_row
        cmp_counter <- cmp_counter + 1L
      }
    }
  }
  
  if (mode == "custom") {
    stop("Custom mode is not yet implemented in this skeleton.")
  }
  
  if (length(res) == 0) {
    out <- data.frame(
      comparison_id = character(),
      comparison_name = character(),
      mode = character(),
      group_id = character(),
      group_test = character(),
      group_ref = character(),
      n_test = integer(),
      n_ref = integer(),
      sample_test = character(),
      sample_ref = character(),
      valid = logical(),
      notes = character(),
      stringsAsFactors = FALSE
    )
    return(out)
  }
  
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}


# ------------------------------------------------------------
# Validate comparisons
# ------------------------------------------------------------
validate_comparisons <- function(comparison_table,
                                 min_n_test = 1L,
                                 min_n_ref = 1L) {
  if (nrow(comparison_table) == 0) {
    warning("No comparisons were defined.")
    return(comparison_table)
  }
  
  comparison_table$valid <- as.logical(comparison_table$valid)
  comparison_table$notes <- as.character(comparison_table$notes)
  
  bad_test <- !is.na(comparison_table$n_test) & comparison_table$n_test < min_n_test
  bad_ref  <- !is.na(comparison_table$n_ref) & comparison_table$n_ref < min_n_ref
  
  comparison_table$valid[bad_test | bad_ref] <- FALSE
  
  comparison_table$notes[bad_test] <- paste(
    trimws(comparison_table$notes[bad_test]),
    sprintf("n_test < %d", min_n_test)
  )
  comparison_table$notes[bad_ref] <- paste(
    trimws(comparison_table$notes[bad_ref]),
    sprintf("n_ref < %d", min_n_ref)
  )
  
  comparison_table$notes <- trimws(comparison_table$notes)
  
  comparison_table
}


# ------------------------------------------------------------
# Build comparison object
# ------------------------------------------------------------
build_comparison_object <- function(comparison_row, metadata, score_matrix, sample_col = "sample") {
  test_samples <- if (!is.na(comparison_row$sample_test) && nzchar(comparison_row$sample_test)) {
    trimws(strsplit(comparison_row$sample_test, ";", fixed = TRUE)[[1]])
  } else {
    character()
  }
  
  ref_samples <- if (!is.na(comparison_row$sample_ref) && nzchar(comparison_row$sample_ref)) {
    trimws(strsplit(comparison_row$sample_ref, ";", fixed = TRUE)[[1]])
  } else {
    character()
  }
  
  keep_samples <- unique(c(test_samples, ref_samples))
  keep_samples <- keep_samples[keep_samples %in% colnames(score_matrix)]
  
  metadata_subset <- metadata[metadata[[sample_col]] %in% keep_samples, , drop = FALSE]
  metadata_subset <- metadata_subset[match(keep_samples, metadata_subset[[sample_col]]), , drop = FALSE]
  
  score_subset <- score_matrix[, keep_samples, drop = FALSE]
  
  list(
    comparison_id = comparison_row$comparison_id,
    comparison_name = comparison_row$comparison_name,
    mode = comparison_row$mode,
    group_id = comparison_row$group_id,
    group_test = comparison_row$group_test,
    group_ref = comparison_row$group_ref,
    test_samples = test_samples,
    ref_samples = ref_samples,
    metadata_subset = metadata_subset,
    score_subset = score_subset,
    comparison_row = comparison_row
  )
}


# ------------------------------------------------------------
# Build analysis matrix for MCD
# ------------------------------------------------------------
build_analysis_matrix <- function(score_matrix,
                                  test_samples,
                                  ref_samples = NULL,
                                  mode = c("global", 
                                           "treatment_vs_control", 
                                           "pairwise_by_treatment"),
                                  aggregate_fun = c("mean", "median")) {
  mode <- match.arg(mode)
  aggregate_fun <- match.arg(aggregate_fun)
  
  agg <- switch(
    aggregate_fun,
    mean = function(x) rowMeans(x, na.rm = TRUE),
    median = function(x) apply(x, 1, stats::median, na.rm = TRUE)
  )
  
  if (mode == "global") {
    return(score_matrix)
  }
  
  test_mat <- score_matrix[, test_samples, drop = FALSE]
  ref_mat  <- score_matrix[, ref_samples, drop = FALSE]
  
  mean_test <- agg(test_mat)
  mean_ref  <- agg(ref_mat)
  delta_mean <- mean_test - mean_ref
  
  out <- cbind(
    mean_test = mean_test,
    mean_ref = mean_ref,
    delta_mean = delta_mean
  )
  
  rownames(out) <- rownames(score_matrix)
  out
}



build_group_analysis_matrix <- function(comparison_subset,
                                        score_matrix,
                                        aggregate_fun = c("mean", "median"),
                                        mode = c("global", 
                                                 "treatment_vs_control", 
                                                 "pairwise_by_treatment")) {
  aggregate_fun <- match.arg(aggregate_fun)
  mode <- match.arg(mode)
  
  agg <- switch(
    aggregate_fun,
    mean = function(x) rowMeans(x, na.rm = TRUE),
    median = function(x) apply(x, 1, stats::median, na.rm = TRUE)
  )
  
  # ----------------------------------------------------------
  # Global mode:
  # For each group, feed the raw per-sample matrix into MCD.
  # If comparison_subset contains one row describing the group,
  # use sample_test as the sample set.
  # ----------------------------------------------------------
  if (mode == "global") {
    if (nrow(comparison_subset) != 1) {
      stop("Global mode expects one row per group in comparison_subset.")
    }
    
    samples <- trimws(strsplit(comparison_subset$sample_test[1], ";", fixed = TRUE)[[1]])
    samples <- samples[samples %in% colnames(score_matrix)]
    
    mat <- score_matrix[, samples, drop = FALSE]
    return(mat)
  }
  
  # ----------------------------------------------------------
  # Non-global modes:
  # Build one column per comparison inside the group.
  # Each column is a delta:
  #   aggregate(test) - aggregate(ref)
  # ----------------------------------------------------------
  out_list <- vector("list", nrow(comparison_subset))
  out_names <- character(nrow(comparison_subset))
  
  for (i in seq_len(nrow(comparison_subset))) {
    cmp <- comparison_subset[i, , drop = FALSE]
    
    test_samples <- trimws(strsplit(cmp$sample_test, ";", fixed = TRUE)[[1]])
    ref_samples  <- trimws(strsplit(cmp$sample_ref,  ";", fixed = TRUE)[[1]])
    
    test_samples <- test_samples[test_samples %in% colnames(score_matrix)]
    ref_samples  <- ref_samples[ref_samples %in% colnames(score_matrix)]
    
    if (length(test_samples) == 0 || length(ref_samples) == 0) {
      delta <- rep(NA_real_, nrow(score_matrix))
    } else {
      test_mat <- score_matrix[, test_samples, drop = FALSE]
      ref_mat  <- score_matrix[, ref_samples, drop = FALSE]
      
      delta <- agg(test_mat) - agg(ref_mat)
    }
    
    out_list[[i]] <- delta
    out_names[i] <- cmp$comparison_name
  }
  
  out <- do.call(cbind, out_list)
  colnames(out) <- make.unique(out_names)
  rownames(out) <- rownames(score_matrix)
  
  out
}


# ------------------------------------------------------------
# Placeholder MCD engine
# Replace this with your actual compute.MCD wrapper
# ------------------------------------------------------------


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# using the grenander ebp
# BiocManager::install("HybridMTest")
# library("HybridMTest")

library(fdrtool)

grenander =  function(F, type=c("decreasing", "increasing"))  # from fdrtool package
  {
    if( !any(class(F) == "ecdf") ) stop("ecdf object required as input!")
    type = match.arg(type)
    if (type == "decreasing")
    {
      # find least concave majorant of ECDF
      ll = fdrtool::gcmlcm(environment(F)$x, environment(F)$y, type="lcm")
    }
    else
    {
      # find greatest convex minorant of ECDF
      l = length(environment(F)$y)
      ll = fdrtool::gcmlcm(environment(F)$x, c(0,environment(F)$y[-l]), type="gcm")
    }
    f.knots = ll$slope.knots
    f.knots = c(f.knots, f.knots[length(f.knots)])
    g = list(F=F,
             x.knots=ll$x.knots,
             F.knots=ll$y.knots,
             f.knots=f.knots)
    class(g) = "grenander"
    return(g)
  }

grenander.ebp =function(p)     # Compute the grenander.ebp from a vector of p-values
  {
    na=is.na(p)
    p.edf=ecdf(p[!na])
    gren.res=grenander(p.edf)
    gren.pdf=approx(gren.res$x.knots,gren.res$f.knots,xout=p)$y
    gren.ebp=min(gren.res$f.knots)/gren.pdf
    ebp.null=pval.pdf=rep(NA,length(p))
    ebp.null[!na]=gren.ebp
    pval.pdf[!na]=gren.pdf   
    return(cbind.data.frame(pval = p,pval.pdf=pval.pdf,ebp.null=ebp.null))
  }

# Langaas, M., Lindqvist, B., 2005. Estimating the proportion of true null hypotheses, 
# with application to DNA microarray data. J.R. Statist. Soc. B67, part4, 555-572. Strimmer, K. 2008. A
# unified approach to false discovery rate estimation. BMC Bioinformatics 9: 303. Strimmer, K.
# 2008. fdrtool: a versatile R package for estimating local and tail area-based false discovery rates.
# Bioinformatics 24: 1461-1462. Grenander, U. (1956). On the theory of mortality measurement, 
# Part II. Skan. Aktuarietidskr, 39, 125–153.

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 


run_single_mcd <- function(analysis_matrix,
                           genes = rownames(analysis_matrix),
                           comparison_info = NULL, 
                           B = 5) {
  
  if (is.null(genes)) {
    genes <- paste0("gene_", seq_len(nrow(analysis_matrix)))
  }

  # ----------------------------------------------------------
  # This is intentionally a placeholder.
  # Replace with your actual MCD implementation.
  
  X =as.matrix(analysis_matrix)
  mcd.res=robustbase::covMcd(X)                                                  # apply MCD
  mhd=mahalanobis(X,mcd.res$center,mcd.res$cov)                                  # compute Mahalanobis distances
  
  k=ncol(X)                                                       # dimension of X
  p.chisq = pchisq(mhd, k, lower.tail=F)
  
  # median scaling to chi-square distribution
  qqr=sort(mhd)/qchisq(1:length(mhd)/(length(mhd)+1),k)     # quantile ratios
  mdn.ratio=median(qqr)                                           # median ratio
  mr.mhd=mhd/mdn.ratio                                            # median ratio scaled Mahalanobis distance
  p.mdn.scl=pchisq(mr.mhd,k,lower.tail=F)                         # p-values for median ratio scaled Mahalanobis distances
  
  gren.ebp = grenander.ebp(p.mdn.scl)
  
  # simulation based
  
  sim.list = replicate(B, {                                                     # loop over sim reps
    idx = sample.int(m, size = m, replace = TRUE, prob = 1-gren.ebp$ebp.null)     # bootstrap X
    sim.X = X[idx, , drop = FALSE]                                              # apply mcd to the simulated data
    
    # null_idx = which(gren.ebp$ebp.null > 0.8)
    # sim.X = X[sample(null_idx, size = m, replace = TRUE), ]
    
    sim.mcd = robustbase::covMcd(sim.X, raw=T)                                                     # mahalanobis distances from the mcd on simulated data
    mahalanobis(sim.X, sim.mcd$center, sim.mcd$cov)
  }, simplify = FALSE)
  
  sim.dist.all = unlist(sim.list, use.names = FALSE)
  
  p.sim = vapply(mhd, function(d) mean(sim.dist.all >= d), numeric(1))          # get the empirical tail probabilities
  
  # Expected output:
  # - one row per gene
  # - at minimum: gene, mahalanobis, pval, fdr, outlier
  # ----------------------------------------------------------

  res <- data.frame(
    gene = genes,
    mahalanobis = mhd,
    p.mdn = p.mdn.scl,
    p.sim = p.sim,
    rank_mahalanobis = rank(-mhd, ties.method = "min"),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(comparison_info)) {
    res$comparison_id <- comparison_info$comparison_id
    res$comparison_name <- comparison_info$comparison_name
    res$mode <- comparison_info$mode
    res$group_test <- comparison_info$group_test
    res$group_ref <- comparison_info$group_ref
  }
  
  res
}


# ------------------------------------------------------------
# Run MCD across comparisons
# ------------------------------------------------------------
run_mcd_per_group <- function(comparison_table,
                                       metadata,
                                       score_matrix,
                                       gene_data = NULL,
                                       aggregate_fun = "mean") {
  per_group <- list()
  qc_list <- list()

  
  if (nrow(comparison_table) == 0) {
    return(list(
      per_group = per_group,
      per_comparison = list(),
      qc_metrics = data.frame(),
      pca_coordinates = NULL
    ))
  }
    
  cmp <- comparison_table[i, , drop = FALSE]
    
  if (!isTRUE(cmp$valid)) {
    
    qc_metrics <- data.frame(
      comparison_id = character(),
      comparison_name = character(),
      group_id = character(),
      n_samples_total = integer(),
      n_features = integer(),
      n_genes = integer(),
      n_missing = integer(),
      n_zero_var_genes = integer(),
      matrix_rank = integer(),
      alpha = numeric(),
      status = character(),
      warning = character(),
      stringsAsFactors = FALSE
    )
    
    return(list(
      per_group = per_group,
      per_comparison = list(),
      qc_metrics = qc_metrics,
      pca_coordinates = NULL
    ))
  }
  
  # ----------------------------------------------------------
  # Determine grouping key
  # If group_id exists in comparison_table, use it directly.
  # Otherwise treat everything as one group.
  # ----------------------------------------------------------
  
  valid_cmp <- comparison_table[comparison_table$valid, , drop = FALSE]
  
  
  if ("group_id" %in% names(valid_cmp)) {
    split_groups <- split(valid_cmp, valid_cmp$group_id)
  } else {
    split_groups <- list(all_samples = valid_cmp)
  }
  
  for (gid in names(split_groups)) {
    cmp_df <- split_groups[[gid]]
    
    mode_i <- unique(cmp_df$mode)
    if (length(mode_i) != 1) {
      warning(sprintf("Group '%s' has multiple modes; using the first.", gid))
      mode_i <- mode_i[1]
    }
    
    analysis_matrix <- build_group_analysis_matrix(
      comparison_subset = cmp_df,
      score_matrix = score_matrix,
      aggregate_fun = aggregate_fun,
      mode = mode_i
    )
    
    genes <- rownames(analysis_matrix)
    if (is.null(genes) && !is.null(gene_data) && "gene" %in% names(gene_data)) {
      genes <- gene_data$gene
    }
    if (is.null(genes)) {
      genes <- paste0("gene_", seq_len(nrow(analysis_matrix)))
    }
    
    # --------------------------------------------------------
    # Run one MCD per group
    # --------------------------------------------------------
    mcd_res <- run_single_mcd(
      analysis_matrix = analysis_matrix,
      genes = genes,
      comparison_info = data.frame(
        comparison_id = paste(cmp_df$comparison_id, collapse = ";"),
        comparison_name = gid,
        mode = mode_i,
        group_test = if ("group_test" %in% names(cmp_df)) paste(unique(cmp_df$group_test), collapse = ";") else NA_character_,
        group_ref  = if ("group_ref"  %in% names(cmp_df)) paste(unique(cmp_df$group_ref),  collapse = ";") else NA_character_,
        stringsAsFactors = FALSE
      )
    )
    
    
    # Tag as group-level output
    mcd_res$group_id <- gid
    mcd_res$n_group_comparisons <- nrow(cmp_df)
    mcd_res$feature_names <- paste(colnames(analysis_matrix), collapse = "; ")
    
    # carry grouping variables if present
    if (!is.null(grouping_vars)) {
      for (v in grouping_vars) {
        if (v %in% names(cmp_df)) {
          vals <- unique(cmp_df[[v]])
          mcd_res[[v]] <- if (length(vals) == 1) vals else paste(vals, collapse = ";")
        }
      }
    }
    
    per_group[[gid]] <- mcd_res
    
    # --------------------------------------------------------
    # QC
    # --------------------------------------------------------
    zero_var <- tryCatch(
      apply(analysis_matrix, 1, stats::var, na.rm = TRUE),
      error = function(e) rep(NA_real_, nrow(analysis_matrix))
    )
    
    qc_row <- data.frame(
      group_id = gid,
      mode = mode_i,
      n_group_comparisons = nrow(cmp_df),
      n_samples_total = sum(unique(c(cmp_df$n_test, cmp_df$n_ref)), na.rm = TRUE),
      n_features = ncol(analysis_matrix),
      n_genes = nrow(analysis_matrix),
      n_missing = sum(is.na(analysis_matrix)),
      n_zero_var_genes = sum(zero_var == 0, na.rm = TRUE),
      matrix_rank = tryCatch(qr(analysis_matrix)$rank, error = function(e) NA_integer_),
      alpha = alpha,
      status = "success",
      warning = "",
      stringsAsFactors = FALSE
    )
    
    if (!is.null(grouping_vars)) {
      for (v in grouping_vars) {
        if (v %in% names(cmp_df)) {
          vals <- unique(cmp_df[[v]])
          qc_row[[v]] <- if (length(vals) == 1) vals else paste(vals, collapse = ";")
        }
      }
    }
    
    qc_list[[length(qc_list) + 1L]] <- qc_row
    
  }
    
  qc_metrics <- if (length(qc_list) > 0) do.call(rbind, qc_list) else data.frame()

  # keep backward compatibility name if you want
  list(
    per_group = per_group,
    qc_metrics = qc_metrics
  )

}

# ------------------------------------------------------------
# Combine all gene-level results
# ------------------------------------------------------------
combine_mcd_results <- function(per_group) {
  if (length(per_group) == 0) return(data.frame())
  out <- do.call(rbind, per_group)
  rownames(out) <- NULL
  out
}


# ------------------------------------------------------------
# Build outlier hits table
# ------------------------------------------------------------
build_outlier_hits <- function(combined_gene_results, fdr_cutoff = 0.05) {
  if (nrow(combined_gene_results) == 0) return(combined_gene_results)
  
  keep <- !is.na(combined_gene_results$fdr) & combined_gene_results$fdr <= fdr_cutoff
  out <- combined_gene_results[keep, , drop = FALSE]
  rownames(out) <- NULL
  out
}


# ------------------------------------------------------------
# Summarize by comparison
# ------------------------------------------------------------
summarize_by_comparison <- function(combined_gene_results) {
  if (nrow(combined_gene_results) == 0) return(data.frame())
  
  split_res <- split(combined_gene_results, combined_gene_results$comparison_id)
  
  out <- lapply(split_res, function(df) {
    first_row <- df[1, , drop = FALSE]
    
    extra_cols <- intersect(
      c("comparison_name", "mode", "group_test", "group_ref", "n_test", "n_ref",
        "cell_line", "timepoint", "treatment", "control"),
      names(df)
    )
    
    ans <- data.frame(
      comparison_id = first_row$comparison_id,
      n_genes = nrow(df),
      n_outliers_fdr_005 = sum(df$fdr <= 0.05, na.rm = TRUE),
      n_outliers_fdr_010 = sum(df$fdr <= 0.10, na.rm = TRUE),
      median_mahalanobis = stats::median(df$mahalanobis, na.rm = TRUE),
      max_mahalanobis = max(df$mahalanobis, na.rm = TRUE),
      min_fdr = min(df$fdr, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    
    for (nm in extra_cols) {
      ans[[nm]] <- first_row[[nm]]
    }
    
    ans
  })
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  
  front <- intersect(
    c("comparison_id", "comparison_name", "mode", "cell_line", "timepoint",
      "group_test", "group_ref", "n_test", "n_ref"),
    names(out)
  )
  
  out[, c(front, setdiff(names(out), front)), drop = FALSE]
}


# ------------------------------------------------------------
# Summarize by gene
# ------------------------------------------------------------
summarize_by_gene <- function(combined_gene_results) {
  if (nrow(combined_gene_results) == 0) return(data.frame())
  
  split_res <- split(combined_gene_results, combined_gene_results$gene)
  
  out <- lapply(split_res, function(df) {
    data.frame(
      gene = df$gene[1],
      n_comparisons_tested = nrow(df),
      n_outlier_fdr_005 = sum(df$fdr <= 0.05, na.rm = TRUE),
      n_outlier_fdr_010 = sum(df$fdr <= 0.10, na.rm = TRUE),
      min_fdr = min(df$fdr, na.rm = TRUE),
      best_rank = min(df$rank_mahalanobis, na.rm = TRUE),
      mean_mahalanobis = mean(df$mahalanobis, na.rm = TRUE),
      max_mahalanobis = max(df$mahalanobis, na.rm = TRUE),
      mean_delta_mean = if ("delta_mean" %in% names(df)) mean(df$delta_mean, na.rm = TRUE) else NA_real_,
      max_abs_delta_mean = if ("delta_mean" %in% names(df)) max(abs(df$delta_mean), na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  
  out[order(out$min_fdr, out$best_rank), , drop = FALSE]
}


# ------------------------------------------------------------
# README sheet contents
# ------------------------------------------------------------
build_readme_sheet <- function() {
  data.frame(
    sheet = c(
      "metadata",
      "comparison_table",
      "analysis_settings",
      "combined_gene_results",
      "outlier_hits",
      "summary_by_comparison",
      "summary_by_gene",
      "qc_metrics",
      "pca_coordinates"
    ),
    description = c(
      "Sample-level metadata used in the analysis.",
      "Defined comparisons, included samples, validity checks, and grouping information.",
      "Pipeline settings used for this run.",
      "Full gene-level MCD results across all comparisons.",
      "Subset of combined_gene_results filtered by FDR threshold.",
      "Summary statistics for each comparison.",
      "Gene-level recurrence summary across comparisons.",
      "Quality control and diagnostic metrics for each comparison.",
      "Optional PCA coordinates for genes within each comparison."
    ),
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# Excel export
# ------------------------------------------------------------
write_mcd_excel <- function(mcd_pipeline, out_file) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required to write Excel output.")
  }
  
  wb <- openxlsx::createWorkbook()
  
  add_sheet <- function(name, x) {
    if (is.null(x)) return(invisible(NULL))
    if (!is.data.frame(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
    if (nrow(x) == 0 && ncol(x) == 0) x <- data.frame(note = "No data")
    
    openxlsx::addWorksheet(wb, name)
    openxlsx::writeDataTable(
      wb,
      sheet = name,
      x = x,
      withFilter = TRUE,
      tableStyle = "TableStyleMedium2"
    )
    openxlsx::freezePane(wb, sheet = name, firstRow = TRUE)
    openxlsx::setColWidths(wb, sheet = name, cols = seq_len(ncol(x)), widths = "auto")
  }
  
  add_sheet("metadata", mcd_pipeline$metadata)
  add_sheet("comparison_table", mcd_pipeline$comparison_table)
  add_sheet("analysis_settings", flatten_settings(mcd_pipeline$analysis_settings))
  add_sheet("combined_gene_results", mcd_pipeline$results$combined_gene_results)
  add_sheet("outlier_hits", mcd_pipeline$results$outlier_hits)
  add_sheet("summary_by_comparison", mcd_pipeline$results$summary_by_comparison)
  add_sheet("summary_by_gene", mcd_pipeline$results$summary_by_gene)
  add_sheet("qc_metrics", mcd_pipeline$results$qc_metrics)
  
  if (!is.null(mcd_pipeline$results$pca_coordinates)) {
    add_sheet("pca_coordinates", mcd_pipeline$results$pca_coordinates)
  }
  
  add_sheet("README", build_readme_sheet())
  
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
  invisible(out_file)
}


# ------------------------------------------------------------
# Main wrapper
# ------------------------------------------------------------
run_mcd_pipeline <- function(final_input,
                             out_dir = getwd(),
                             alpha = 0.75,
                             p_adjust_method = "fdr",
                             min_n_test = 1L,
                             min_n_ref = 1L,
                             save_excel = TRUE,
                             save_rds = TRUE,
                             excel_file = "mcd_analysis.xlsx",
                             rds_file = "mcd_pipeline.rds") {
  validate_final_input(final_input)
  
  aggregate_fun = final_input$preprocessing$aggregation_fun
  score_matrix <- as.matrix(final_input$score_matrix)
  storage.mode(score_matrix) <- "numeric"
  
  metadata <- final_input$metadata
  
  idx <- match(colnames(score_matrix), metadata$sample)
  if (any(is.na(idx))) {
    stop("Could not align metadata to score_matrix columns.")
  }
  metadata = metadata[idx, , drop = FALSE]
  
  comparison_mode <- choose_comparison_mode()
  
  if (comparison_mode != "global") {
    grouping_vars <- choose_grouping_vars(metadata)
    treatment_col = choose_treatment_vars(metadata)
  }
  
  if (comparison_mode == "treatment_vs_control") {
    treatment_col = choose_treatment_vars(metadata)
    control_label <- choose_control_label(metadata, treatment_col)
  }
  
  cat("\nDefining comparisons...\n")
  
  comparison_table <- define_comparisons(
    metadata = metadata,
    mode = comparison_mode,
    grouping_vars = grouping_vars,
    treatment_col = treatment_col,
    sample_col = 'sample',
    control_label = control_label
  )
  
  # comparison_table <- validate_comparisons(
  #   comparison_table = comparison_table,
  #   min_n_test = min_n_test,
  #   min_n_ref = min_n_ref
  # )

  cat(sprintf("Defined %d comparison(s).\n", nrow(comparison_table)))
  cat(sprintf("Valid comparisons: %d\n", sum(comparison_table$valid, na.rm = TRUE)))
  
  run_res <- run_mcd_per_group(
    comparison_table = comparison_table,
    metadata = metadata,
    score_matrix = score_matrix,
    gene_data = final_input$gene_data %||% NULL,
    aggregate_fun = aggregate_fun,
    sample_col = 'sample'
  )
  
  combined_gene_results <- combine_mcd_results(run_res$per_group)
  outlier_hits <- build_outlier_hits(combined_gene_results, fdr_cutoff = 0.05)
  summary_cmp <- summarize_by_comparison(combined_gene_results)
  summary_gene <- summarize_by_gene(combined_gene_results)
  
  mcd_pipeline <- list(
    call = match.call(),
    timestamp = Sys.time(),
    input_type = final_input$input_type %||% NA_character_,
    selected_metric = final_input$selected_metric %||% NA_character_,
    source_info = final_input$source_info %||% NA,
    
    metadata = metadata,
    score_matrix = score_matrix,
    gene_data = final_input$gene_data %||% NULL,
    
    analysis_settings = list(
      comparison_mode = comparison_mode,
      grouping_vars = grouping_vars,
      treatment_col = treatment_col,
      sample_col = sample_col,
      control_label = control_label,
      aggregate_fun = aggregate_fun,
      mcd_alpha = alpha,
      p_adjust_method = p_adjust_method,
      min_n_test = min_n_test,
      min_n_ref = min_n_ref,
      save_excel = save_excel,
      save_rds = save_rds
    ),
    
    comparisons = comparison_table$comparison_id,
    comparison_table = comparison_table,
    
    results = list(
      per_comparison = run_res$per_comparison,
      combined_gene_results = combined_gene_results,
      summary_by_comparison = summary_cmp,
      summary_by_gene = summary_gene,
      outlier_hits = outlier_hits,
      pca_coordinates = run_res$pca_coordinates,
      qc_metrics = run_res$qc_metrics
    ),
    
    export_paths = list(
      excel = NULL,
      rds = NULL
    )
  )
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (isTRUE(save_excel)) {
    excel_path <- file.path(out_dir, excel_file)
    write_mcd_excel(mcd_pipeline, excel_path)
    mcd_pipeline$export_paths$excel <- excel_path
  }
  
  if (isTRUE(save_rds)) {
    rds_path <- file.path(out_dir, rds_file)
    saveRDS(mcd_pipeline, rds_path)
    mcd_pipeline$export_paths$rds <- rds_path
  }
  
  cat("\nMCD pipeline complete.\n")
  if (!is.null(mcd_pipeline$export_paths$excel)) {
    cat("Excel output: ", mcd_pipeline$export_paths$excel, "\n", sep = "")
  }
  if (!is.null(mcd_pipeline$export_paths$rds)) {
    cat("RDS output: ", mcd_pipeline$export_paths$rds, "\n", sep = "")
  }
  
  return(mcd_pipeline)
}