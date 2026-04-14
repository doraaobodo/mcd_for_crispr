###### HELPERS #######

check_required_packages = function(pkgs) {
  missing = pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  
  if (length(missing) == 0) {
    return(invisible(TRUE))
  }
  
  show_header("Missing Required Packages")
  cat("The following packages are required but not installed:\n")
  for (pkg in missing) cat(" - ", pkg, "\n", sep = "")
  
  cat("\nInstall them in R with:\n")
  cat(sprintf('install.packages(c(%s))\n',
              paste(sprintf('"%s"', missing), collapse = ", ")))
  cat("\n")
  
  pause("Press [Enter] to exit...")
  quit(status = 1)
}


get_script_path = function() {
  args = commandArgs(trailingOnly = FALSE)
  file_arg = "--file="
  idx = grep(file_arg, args)
  
  if (length(idx) > 0) {
    return(normalizePath(sub(file_arg, "", args[idx][1]), winslash = "/", mustWork = TRUE))
  }
  
  NULL
}

get_script_dir = function() {
  script_path = get_script_path()
  if (!is.null(script_path)) {
    return(dirname(script_path))
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

# generate a pause for user input in launcher
pause = function(prompt_text = "Press Enter to continue...") {
  if (interactive()) {
    return(invisible(readline(prompt = prompt_text)))
  } else {
    cat(prompt_text)
    return(invisible(readLines("stdin", n = 1)))
  }
}

say_line = function(char = "=", n = 60) {
  cat(paste(rep(char, n), collapse = ""), "\n", sep = "")
}

show_header = function(title) {
  cat("\n")
  say_line("=")
  cat(title, "\n")
  say_line("=")
}

ask_yes_no = function(prompt, default = TRUE) {
  suffix = if (default) " [Y/n]: " else " [y/N]: "
  
  repeat {
    
    if (interactive()) {
      ans = trimws(tolower(readline(paste0(prompt, suffix))))
    } else {
      cat(paste0(prompt, suffix))
      ans = trimws(tolower(readLines("stdin", n=1)))
    }
    
    if (ans == "") return(default)
    if (ans %in% c("y", "yes")) return(TRUE)
    if (ans %in% c("n", "no")) return(FALSE)
    
    cat("Please enter y or n.\n")
  }
}


choose_input_type = function() {
  show_header("Select Input Type")
  cat("What type of file are you providing?\n\n")
  cat("  1) Raw counts matrix\n")
  cat("  2) MAGeCK RRA gene summary results\n")
  cat("  3) MAGeCK MLE results\n\n")
  
  #we want to repeat until we get the right input
  repeat {
    
    if (interactive()) {
      ans = trimws(readline("Enter choice number: "))
    } else {
      cat("Enter choice number: ")
      ans = trimws(readLines("stdin", n=1))
    }

    if (ans == "1") return("raw_counts")
    if (ans == "2") return("mageck_rra")
    if (ans == "3") return("mageck_mle")
    
    cat("Invalid choice. Please enter 1, 2, or 3.\n")
  }
}


choose_input_source_mode = function() {
  show_header("Select Input Source")
  cat("How would you like to provide the input?\n\n")
  cat("  1) Single file (e.g., counts data, MAGECK MLE results for one screen)\n")
  cat("  2) Multiple files (e.g., MAGeCK RRA/MLE gene summary results for multiple screens)\n")
  cat("  3) Folder containing files (e.g., folder containing only results files \n with the same structure \n\n")
  
  repeat {
    if (interactive()) {
      ans = trimws(readline("Enter choice number: "))
    } else {
      cat("Enter choice number: ")
      ans = trimws(readLines("stdin", n=1))
    }
    
    if (ans == "1") return("single")
    if (ans == "2") return("multiple")
    if (ans == "3") return("directory")
    cat("Invalid choice. Please enter 1, 2, or 3.\n")
  }
}

choose_input_files = function(
    prompt_text = "Enter full paths separated by commas: ",
    caption = "Select input files"
) {
  repeat {
    selected = character(0)
    
    if (.Platform$OS.type == "windows") {
      selected = tryCatch(
        utils::choose.files(caption = caption, multi = TRUE),
        error = function(e) character(0)
      )
    }
    
    if (length(selected) > 0 && all(file.exists(selected))) {
      return(normalizePath(selected, winslash = "/", mustWork = TRUE))
    }
    
    cat("File chooser unavailable or no files selected.\n")
    if (interactive()) {
      ans = trimws(readline(prompt_text))
    } else {
      cat(prompt_text)
      ans = trimws(readLines("stdin", n=1))
    }
    parts = trimws(strsplit(ans, ",", fixed = TRUE)[[1]])
    parts = path.expand(parts[nzchar(parts)])
    
    if (length(parts) > 0 && all(file.exists(parts))) {
      return(normalizePath(parts, winslash = "/", mustWork = TRUE))
    }
    
    cat("One or more files were not found. Please try again.\n")
  }
}

collect_input_paths = function() {
  mode = choose_input_source_mode()
  
  if (mode == "single") {
    paths = choose_input_files(
      prompt_text = "Enter full path to input file: ",
      caption = "Select input file"
    )
    return(list(mode = mode, paths = paths))
  }
  
  if (mode == "multiple") {
    paths = choose_input_files(
      prompt_text = "Enter full paths to input files, separated by commas: ",
      caption = "Select input files"
    )
    return(list(mode = mode, paths = paths))
  }
  
  dir_path = choose_output_dir(
    prompt_text = "Enter folder containing input files: ",
    caption = "Select input folder",
    create = FALSE
  )
  
  files = list.files(dir_path, full.names = TRUE)
  files = files[file.info(files)$isdir %in% FALSE]
  
  if (length(files) == 0) {
    stop("No files were found in the selected directory.")
  }
  
  list(mode = mode, paths = normalizePath(files, winslash = "/", mustWork = TRUE))
}


load_counts_data = function(file_path) {
  ext = tolower(tools::file_ext(file_path))
  
  dat = switch(
    ext,
    "csv" = read.csv(file_path, stringsAsFactors = FALSE, check.names = FALSE),
    "tsv" = read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE),
    "txt" = read.delim(file_path, stringsAsFactors = FALSE, check.names = FALSE),
    stop("Unsupported file type: ", ext, ". Please use csv, tsv, or txt.")
  )
  
  dat
}

parse_one_mageck_rra_file = function(path, metric_col = "neg|lfc", gene_col = "id") {
  dat = load_counts_data(path)
  
  required_cols = c(gene_col, metric_col)
  missing_cols = setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop(
      basename(path),
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out = dat[, c(gene_col, metric_col), drop = FALSE]
  colnames(out) = c("gene", "score")
  out$gene = as.character(out$gene)
  out
}

parse_mageck_rra_set = function(paths, metric_col = "neg|lfc") {
  score_list = vector("list", length(paths))
  sample_names = character(length(paths))
  
  for (i in seq_along(paths)) {
    path = paths[i]
    one = parse_one_mageck_rra_file(path, metric_col = metric_col)
    
    sample_label = tools::file_path_sans_ext(basename(path))
    sample_label = gsub('.gene_summary|_gene_summary', '', sample_label)
    sample_label = gsub('.rra|_rra', '', sample_label)
    
    colnames(one)[2] = sample_label
    
    score_list[[i]] = one
    sample_names[i] = sample_label
  }
  
  merged = Reduce(
    function(x, y) merge(x, y, by = "gene", all = TRUE),
    score_list
  )
  
  rownames(merged) = merged$gene
  score_matrix = merged[, sample_names, drop = FALSE]
  
  list(
    input_type = "mageck_rra",
    gene_data = data.frame(gene = merged$gene, stringsAsFactors = FALSE),
    score_matrix = score_matrix,
    sample_names = colnames(score_matrix),
    source_info = data.frame(file = paths, label = sample_names, stringsAsFactors = FALSE),
    selected_metric = metric_col
  )
}

parse_one_mageck_mle_file = function(path) {
  dat = load_counts_data(path)
  
  required_cols = "Gene"
  missing_cols = setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop(
      basename(path),
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  z_cols = grep("\\|z$", colnames(dat), value = TRUE)
  if (length(z_cols) == 0) {
    stop(basename(path), " has no columns ending in '|z'.")
  }
  
  out = dat[, c("Gene", z_cols), drop = FALSE]
  colnames(out)[1] = "gene"
  out$gene = as.character(out$gene)
  
  file_stub = tools::file_path_sans_ext(basename(path))
  new_names = sub("\\|z$", "", z_cols)
  colnames(out)[-1] = new_names
  
  out
}

parse_mageck_mle_set = function(paths) {
  parsed = lapply(paths, parse_one_mageck_mle_file)
  
  merged = Reduce(
    function(x, y) merge(x, y, by = "gene", all = TRUE),
    parsed
  )
  
  rownames(merged) = merged$gene
  score_matrix = merged[, setdiff(colnames(merged), "gene"), drop = FALSE]
  
  list(
    input_type = "mageck_mle",
    gene_data = data.frame(gene = merged$gene, stringsAsFactors = FALSE),
    score_matrix = score_matrix,
    sample_names = colnames(score_matrix),
    source_info = data.frame(file = paths, stringsAsFactors = FALSE),
    selected_metric = "z"
  )
}

show_vector_choices = function(x) {
  for (i in seq_along(x)) {
    cat(sprintf("  %d) %s\n", i, x[i]))
  }
}

ask_single_column = function(prompt, choices, allow_skip = FALSE) {
  repeat {
    cat("\n", prompt, "\n", sep = "")
    show_vector_choices(choices)
    
    if (allow_skip) {
      cat("  0) Skip\n")
    }
    
    if (interactive()) {
      ans = trimws(readline("Enter choice number: "))
    } else {
      cat("Enter choice number: ")
      ans = trimws(readLines("stdin", n=1))
    }
    
    if (allow_skip && ans == "0") {
      return(NA_character_)
    }
    
    if (ans %in% as.character(seq_along(choices))) {
      return(choices[as.integer(ans)])
    }
    
    cat("Invalid choice. Please enter one of the listed numbers.\n")
  }
}

ask_multiple_columns = function(prompt, choices) {
  repeat {
    cat("\n", prompt, "\n", sep = "")
    show_vector_choices(choices)
    cat("Enter one or more numbers separated by commas.\n")
    cat("Example: 3,4,5\n")
    
    if (interactive()) {
      ans = trimws(readline("Selection: "))
    } else {
      cat("Selection: ")
      ans = trimws(readLines("stdin", n=1))
    }
    
    if (!nzchar(ans)) {
      cat("Please select at least one column.\n")
      next
    }
    
    idx = trimws(strsplit(ans, ",", fixed = TRUE)[[1]])
    
    if (all(idx %in% as.character(seq_along(choices)))) {
      idx_num = unique(as.integer(idx))
      return(choices[idx_num])
    }
    
    cat("Invalid selection. Please enter valid column numbers separated by commas.\n")
  }
}

guess_counts_columns = function(cnts) {
  nms = colnames(cnts)
  low = tolower(nms)
  
  is_numeric_col = vapply(cnts, is.numeric, logical(1))
  
  guide_candidates = nms[
    grepl("guide|grna|sgRNA|sgrna|barcode|construct", nms, ignore.case = TRUE)
  ]
  
  gene_candidates = nms[
    grepl("gene|symbol|hgnc|target", nms, ignore.case = TRUE)
  ]
  
  sample_candidates = nms[is_numeric_col]
  
  list(
    guide_candidates = unique(guide_candidates),
    gene_candidates = unique(gene_candidates),
    sample_candidates = unique(sample_candidates)
  )
}

preview_column_guesses = function(guesses) {
  cat("\nDetected likely columns:\n")
  
  cat("\nGuide column candidates:\n")
  if (length(guesses$guide_candidates) > 0) {
    cat(paste0(" - ", guesses$guide_candidates), sep = "\n")
  } else {
    cat(" - None detected\n")
  }
  
  cat("\nGene column candidates:\n")
  if (length(guesses$gene_candidates) > 0) {
    cat(paste0(" - ", guesses$gene_candidates), sep = "\n")
  } else {
    cat(" - None detected\n")
  }
  
  cat("\nSample/count column candidates:\n")
  if (length(guesses$sample_candidates) > 0) {
    cat(paste0(" - ", guesses$sample_candidates), sep = "\n")
  } else {
    cat(" - None detected\n")
  }
  
  cat("\n")
}

validate_counts_mapping = function(cnts, mapping) {
  required_names = c("guide_col", "gene_col", "sample_cols")
  missing_names = setdiff(required_names, names(mapping))
  
  if (length(missing_names) > 0) {
    stop("Mapping object is missing required entries: ",
         paste(missing_names, collapse = ", "))
  }
  
  all_cols = colnames(cnts)
  
  if (!mapping$guide_col %in% all_cols) {
    stop("Mapped guide column not found in counts data: ", mapping$guide_col)
  }
  
  if (!mapping$gene_col %in% all_cols) {
    stop("Mapped gene column not found in counts data: ", mapping$gene_col)
  }
  
  if (length(mapping$sample_cols) == 0) {
    stop("No sample/count columns were selected.")
  }
  
  if (!all(mapping$sample_cols %in% all_cols)) {
    bad = mapping$sample_cols[!mapping$sample_cols %in% all_cols]
    stop("Mapped sample/count columns not found in counts data: ",
         paste(bad, collapse = ", "))
  }
  
  if (mapping$guide_col == mapping$gene_col) {
    stop("Guide column and gene column cannot be the same.")
  }
  
  overlap = intersect(c(mapping$guide_col, mapping$gene_col), mapping$sample_cols)
  if (length(overlap) > 0) {
    stop("Annotation columns cannot also be sample/count columns: ",
         paste(overlap, collapse = ", "))
  }
  
  non_numeric_sample_cols = mapping$sample_cols[
    !vapply(cnts[mapping$sample_cols], is.numeric, logical(1))
  ]
  
  if (length(non_numeric_sample_cols) > 0) {
    stop("These selected sample/count columns are not numeric: ",
         paste(non_numeric_sample_cols, collapse = ", "))
  }
  
  invisible(TRUE)
}

extract_mapped_counts = function(cnts, mapping) {
  annotation = cnts[, c(mapping$guide_col, mapping$gene_col), drop = FALSE]
  counts_mat = cnts[, mapping$sample_cols, drop = FALSE]
  rownames(counts_mat) = cnts[[mapping$guide_col]]

  
  list(
    annotation = annotation,
    counts = counts_mat
  )
}

map_counts_columns = function(cnts) {
  show_header("Map Counts Columns")
  
  all_cols = colnames(cnts)
  guesses = guess_counts_columns(cnts)
  
  cat("The counts table has the following columns:\n\n")
  show_vector_choices(all_cols)
  cat("\n")
  
  preview_column_guesses(guesses)
  
  guide_default = if (length(guesses$guide_candidates) >= 1) guesses$guide_candidates[1] else NA_character_
  gene_default  = if (length(guesses$gene_candidates) >= 1)  guesses$gene_candidates[1]  else NA_character_
  sample_default = guesses$sample_candidates
  
  if (length(sample_default) == 0) {
    cat("No numeric sample/count columns were automatically detected.\n")
    cat("You will need to select them manually.\n")
  }
  
  use_guess = FALSE
  if (!is.na(guide_default) || !is.na(gene_default) || length(sample_default) > 0) {
    cat("Proposed mapping:\n")
    cat("  Guide column: ", ifelse(is.na(guide_default), "<none>", guide_default), "\n", sep = "")
    cat("  Gene column:  ", ifelse(is.na(gene_default), "<none>", gene_default), "\n", sep = "")
    cat("  Sample columns:\n")
    if (length(sample_default) > 0) {
      cat(paste0("   - ", sample_default), sep = "\n")
    } else {
      cat("   - <none>\n")
    }
    cat("\n")
    
    use_guess = ask_yes_no("Use this proposed mapping?", default = TRUE)
  }
  
  if (use_guess) {
    guide_col = guide_default
    gene_col = gene_default
    sample_cols = sample_default
  } else {
    guide_col = ask_single_column(
      prompt = "Select the guide column:",
      choices = all_cols,
      allow_skip = FALSE
    )
    
    gene_col = ask_single_column(
      prompt = "Select the gene column:",
      choices = all_cols,
      allow_skip = FALSE
    )
    
    remaining_cols = setdiff(all_cols, unique(c(guide_col, gene_col)))
    
    sample_cols = ask_multiple_columns(
      prompt = "Select the sample/count columns:",
      choices = remaining_cols
    )
  }
  
  mapping = list(
    guide_col = guide_col,
    gene_col = gene_col,
    sample_cols = sample_cols
  )
  
  validate_counts_mapping(cnts, mapping)
  
  show_header("Counts Column Mapping")
  cat("Guide column: ", mapping$guide_col, "\n", sep = "")
  cat("Gene column:  ", mapping$gene_col, "\n", sep = "")
  cat("Sample columns:\n")
  cat(paste0(" - ", mapping$sample_cols), sep = "\n")
  cat("\n")
  
  if (!ask_yes_no("Is this mapping correct?", default = TRUE)) {
    cat("Let's try again.\n")
    return(map_counts_columns(cnts))
  }
  
  mapping
}

validate_sample_metadata_alignment = function(sample_names, meta) {
  required_col = "sample"
  
  if (!required_col %in% colnames(meta)) {
    stop("Metadata must contain a 'sample' column.")
  }
  
  meta_samples = as.character(meta$sample)
  
  missing_in_meta = setdiff(sample_names, meta_samples)
  missing_in_counts = setdiff(meta_samples, sample_names)
  
  if (length(missing_in_meta) > 0) {
    stop(
      "These count samples are missing from metadata: ",
      paste(missing_in_meta, collapse = ", ")
    )
  }
  
  if (length(missing_in_counts) > 0) {
    stop(
      "These metadata samples are not present in counts data: ",
      paste(missing_in_counts, collapse = ", ")
    )
  }
  
  invisible(TRUE)
}

align_metadata_to_samples = function(sample_names, meta) {
  meta = as.data.frame(meta, stringsAsFactors = FALSE)
  idx = match(sample_names, meta$sample)
  
  if (any(is.na(idx))) {
    stop("Could not align metadata to sample names.")
  }
  
  meta[idx, , drop = FALSE]
}

aggregate_raw_counts_to_gene = function(parsed_raw, aggregation_fun = c("sum", "mean", "median", "none")) {
  aggregation_fun = match.arg(aggregation_fun)
  
  if (is.null(parsed_raw$annotation) || is.null(parsed_raw$counts_matrix)) {
    stop("parsed_raw must contain 'annotation' and 'counts_matrix'.")
  }
  
  mapping = parsed_raw$mapping
  gene_col = mapping$gene_col
  ann = parsed_raw$annotation
  cnt = parsed_raw$counts_matrix
  
  if (ncol(ann) < 2) {
    stop("Annotation must contain at least guide and gene columns.")
  }
  
  genes = as.character(parsed_raw$annotation[[gene_col]])
  
  if (any(is.na(genes)) || any(trimws(genes) == "")) {
    stop("Gene column contains missing or empty values.")
  }
  
  # TODO: split by finding gene name in sgrna col
  split_idx = split(seq_len(nrow(cnt)), genes)
  
  agg_list = lapply(split_idx, function(i) {
    block = cnt[i, , drop = FALSE]
    
    if (aggregation_fun == "sum") {
      colSums(block, na.rm = TRUE)
    } else if (aggregation_fun == "mean") {
      colMeans(block, na.rm = TRUE)
    } else if (aggregation_fun == "median"){
      apply(block, 2, median, na.rm = TRUE)
    } else{
      block
    }
  })
  
  gene_mat = do.call(rbind, agg_list)
  gene_mat = as.data.frame(gene_mat, check.names = FALSE)
  gene_mat$gene = rownames(gene_mat)
  rownames(gene_mat) = NULL
  
  gene_data = gene_mat["gene"]
  score_matrix = gene_mat[, setdiff(colnames(gene_mat), "gene"), drop = FALSE]
  
  rownames(score_matrix) = gene_data$gene
  
  list(
    gene_data = gene_data,
    score_matrix = score_matrix
  )
}

normalize_gene_counts = function(score_matrix, method = c("log2", "cpm_log2", "none")) {
  method = match.arg(method)
  
  mat = as.matrix(score_matrix)
  storage.mode(mat) = "numeric"
  
  if (method == "none") {
    return(as.data.frame(mat, check.names = FALSE))
  }
  
  if (method == "log2") {
    mat = log2(mat + 1)
    return(as.data.frame(mat, check.names = FALSE))
  }
  
  lib_sizes = colSums(mat, na.rm = TRUE)
  
  if (any(lib_sizes == 0)) {
    stop("One or more sample columns have total count zero; cannot compute CPM normalization.")
  }
  
  mat = sweep(mat, 2, lib_sizes / 1e6, FUN = "/")
  mat = log2(mat + 1)
  
  as.data.frame(mat, check.names = FALSE)
}

standardize_columns_to_z = function(score_matrix) {
  mat = as.matrix(score_matrix)
  storage.mode(mat) = "numeric"
  
  zmat = apply(mat, 2, function(x) {
    s = stats::sd(x, na.rm = TRUE)
    m = mean(x, na.rm = TRUE)
    
    if (is.na(s) || s == 0) {
      rep(0, length(x))
    } else {
      (x - m) / s
    }
  })
  
  zmat = as.matrix(zmat)
  rownames(zmat) = rownames(mat)
  
  as.data.frame(zmat, check.names = FALSE)
}

build_gene_score_matrix_from_raw_counts = function(
    parsed_raw,
    meta,
    aggregation_fun = c("sum", "mean", "median", "none"),
    normalize_method = c("cpm_log2", "log2", "none"),
    make_z = TRUE
) {
  aggregation_fun = match.arg(aggregation_fun)
  normalize_method = match.arg(normalize_method)
  
  if (is.null(parsed_raw$sample_names)) {
    stop("parsed_raw must contain 'sample_names'.")
  }
  
  validate_sample_metadata_alignment(parsed_raw$sample_names, meta)
  meta_aligned = align_metadata_to_samples(parsed_raw$sample_names, meta)
  
  aggregated = aggregate_raw_counts_to_gene(
    parsed_raw = parsed_raw,
    aggregation_fun = aggregation_fun
  )
  
  norm_scores = normalize_gene_counts(
    score_matrix = aggregated$score_matrix,
    method = normalize_method
  )
  
  rownames(norm_scores) = aggregated$gene_data$gene
  
  if (make_z) {
    final_scores = standardize_columns_to_z(norm_scores)
  } else {
    final_scores = norm_scores
  }
  
  rownames(final_scores) = aggregated$gene_data$gene
  
  list(
    input_type = "raw_counts",
    gene_data = aggregated$gene_data,
    score_matrix = final_scores,
    sample_names = colnames(final_scores),
    source_info = data.frame(
      sample = parsed_raw$sample_names,
      stringsAsFactors = FALSE
    ),
    selected_metric = if (make_z) "gene_level_z_from_raw_counts" else "normalized_gene_counts",
    metadata = meta_aligned,
    preprocessing = list(
      aggregation_fun = aggregation_fun,
      normalize_method = normalize_method,
      make_z = make_z
    )
  )
}

choose_raw_count_processing_options = function() {
  show_header("Raw Count Processing Options")
  
  cat("Choose how guide-level raw counts should be converted to gene-level scores.\n\n")
  
  cat("Aggregation method:\n")
  cat("  1) Sum guide counts within each gene\n")
  cat("  2) Mean guide counts within each gene\n")
  cat("  3) Median guide counts within each gene\n")
  cat("  4) No aggregation\n\n")
  
  repeat {
    
    if (interactive()) {
      ans = trimws(readline("Choose aggregation method [1]: "))
    } else {
      cat("Choose aggregation method [1]: ")
      ans = trimws(readLines("stdin", n=1))
    }
    
    if (ans == "" || ans == "1") {
      aggregation_fun = "sum"
      break
    }
    if (ans == "2") {
      aggregation_fun = "mean"
      break
    }
    if (ans == "3") {
      aggregation_fun = "median"
      break
    }
    if (ans == "4") {
      aggregation_fun = "none"
      break
    }
    
    cat("Invalid choice. Please enter 1, 2, 3, or 4.\n")
  }
  
  cat("\nNormalization method:\n")
  cat("  1) CPM + log2(x + 1)\n")
  cat("  2) log2(x + 1)\n")
  cat("  3) None\n\n")
  
  repeat {
    
    if (interactive()) {
      ans = trimws(readline("Choose normalization method [1]: "))
    } else {
      cat("Choose normalization method [1]: ")
      ans = trimws(readLines("stdin", n=1))
    }
    
    if (ans == "" || ans == "1") {
      normalize_method = "cpm_log2"
      break
    }
    if (ans == "2") {
      normalize_method = "log2"
      break
    }
    if (ans == "3") {
      normalize_method = "none"
      break
    }
    cat("Invalid choice. Please enter 1, 2, or 3.\n")
  }
  
  make_z = ask_yes_no(
    "Convert normalized gene-level values to Z-scores across genes within each sample?",
    default = TRUE
  )
  
  list(
    aggregation_fun = aggregation_fun,
    normalize_method = normalize_method,
    make_z = make_z
  )
}

parse_raw_counts_set = function(path) {
  dat = load_counts_data(path)
  mapping = map_counts_columns(dat)
  mapped = extract_mapped_counts(dat, mapping)
  
  list(
    input_type = "raw_counts",
    raw_table = dat,
    annotation = mapped$annotation,
    counts_matrix = mapped$counts,
    sample_names = colnames(mapped$counts),
    mapping = mapping
  )
}

parse_input_set = function(input_type, paths) {
  if (input_type == "raw_counts") {
    if (length(paths) != 1) {
      stop("Raw counts input currently expects exactly one file.")
    }
    return(parse_raw_counts_set(paths))
  }
  
  if (input_type == "mageck_rra") {
    return(parse_mageck_rra_set(paths, metric_col = "neg|lfc"))
  }
  
  if (input_type == "mageck_mle") {
    return(parse_mageck_mle_set(paths))
  }
  
  stop("Unsupported input type: ", input_type)
}

show_standardized_matrix_summary = function(parsed_input) {
  show_header("Standardized Matrix Summary")
  
  if (!is.null(parsed_input$score_matrix)) {
    cat("Rows (genes): ", nrow(parsed_input$score_matrix), "\n", sep = "")
    cat("Columns (samples/comparisons): ", ncol(parsed_input$score_matrix), "\n\n", sep = "")
    cat("Column names:\n")
    cat(paste0(" - ", colnames(parsed_input$score_matrix)), sep = "\n")
    cat("\n")
    return(invisible(NULL))
  }
  
  if (!is.null(parsed_input$counts_matrix)) {
    cat("Raw counts parsed successfully.\n")
    cat("Rows: ", nrow(parsed_input$counts_matrix), "\n", sep = "")
    cat("Columns: ", ncol(parsed_input$counts_matrix), "\n\n", sep = "")
    cat("Sample names:\n")
    cat(paste0(" - ", colnames(parsed_input$counts_matrix)), sep = "\n")
    cat("\n")
  }
}

show_final_score_matrix_summary = function(final_input) {
  show_header("Final Gene Score Matrix Summary")
  
  cat("Input type: ", final_input$input_type, "\n", sep = "")
  cat("Selected metric: ", final_input$selected_metric, "\n", sep = "")
  cat("Genes: ", nrow(final_input$score_matrix), "\n", sep = "")
  cat("Samples/comparisons: ", ncol(final_input$score_matrix), "\n\n", sep = "")
  
  cat("Column names:\n")
  cat(paste0(" - ", colnames(final_input$score_matrix)), sep = "\n")
  cat("\n")
}

split_sample_tokens = function(sample_name) {
  tokens = unlist(strsplit(sample_name, "[-_.]+", perl = TRUE))
  tokens = tokens[nzchar(tokens)]
  tokens
}

clean_token_core = function(token) {
  sub("\\d+$", "", token)
}

detect_treatment_type = function(treatment) {
  if (is.na(treatment) || !nzchar(treatment)) {
    return(NA_character_)
  }
  
  trt = toupper(treatment)
  if (trt %in% c("DMSO", "CTRL", "CONTROL", "UNTREATED", "VEHICLE", "PBS")) {
    return("control")
  }
  
  "treatment"
}

parse_sample_name = function(sample_name) {
  tokens = split_sample_tokens(sample_name)
  tok_upper = toupper(tokens)
  
  out = list(
    sample = sample_name,
    cell_line = NA_character_,
    treatment = NA_character_,
    treatment_type = NA_character_,
    rep = NA_integer_,
    timepoint = NA_character_
  )
  
  if (length(tokens) == 0) {
    return(out)
  }
  
  time_idx = grep("^(T|TP|D|DAY|TIME|TIMEPOINT).*$", tok_upper)
  if (length(time_idx) > 0) {
    out$timepoint = tokens[time_idx[1]]
    time_idx = time_idx[1]
  }
  
  rep_idx = grep("^(REP|R|REPLICATE).*$", tok_upper)
  if (length(rep_idx) > 0) {
    out$rep = suppressWarnings(
      as.integer(sub("^(REP|R|REPLICATE)", "", tok_upper[rep_idx[1]]))
    )
    rep_idx=rep_idx[1]
  }
  
  if (length(tokens) >= 2) {
    out$cell_line = tokens[1]
  } else {
    out$cell_line = tokens[1]
  }
  
  remaining_idx = seq_along(tokens)
  remaining_idx = setdiff(remaining_idx, c(time_idx, rep_idx))
  
  if (length(tokens) >= 2) {
    remaining_idx = setdiff(remaining_idx, 1)
  }
  
  if (length(remaining_idx) > 0) {
    trt_tok = tokens[remaining_idx[length(remaining_idx)]]
    trt_core = clean_token_core(trt_tok)
    out$treatment = toupper(trt_core)
    
    if (is.na(out$rep) && grepl("\\d+$", trt_tok)) {
      out$rep = suppressWarnings(as.integer(sub("^.*?(\\d+)$", "\\1", trt_tok)))
    }
  }
  
  out$treatment_type = detect_treatment_type(out$treatment)
  
  out
}

autodetect_metadata_from_samples = function(sample_names) {
  parsed = lapply(sample_names, parse_sample_name)
  
  meta = data.frame(
    sample = vapply(parsed, function(x) x$sample, character(1)),
    cell_line = vapply(parsed, function(x) x$cell_line, character(1)),
    treatment = vapply(parsed, function(x) x$treatment, character(1)),
    treatment_type = vapply(parsed, function(x) x$treatment_type, character(1)),
    rep = vapply(parsed, function(x) {
      if (is.na(x$rep)) NA_integer_ else as.integer(x$rep)
    }, integer(1)),
    timepoint = vapply(parsed, function(x) x$timepoint, character(1)),
    stringsAsFactors = FALSE
  )
  
  meta[] = lapply(meta, function(col) {
    if (is.character(col)) {
      col[trimws(col) == ""] = NA_character_
    }
    col
  })
  
  meta
}

edit_metadata_interactively = function(meta) {
  show_header("Metadata Preview")
  cat("Metadata inferred from sample names looks like:\n\n")
  print(meta)
  cat("\nRequired columns: sample, cell_line, treatment, treatment_type\n")
  cat("You may now edit the metadata table.\n\n")
  
  pause("Press Enter to open the metadata editor...")
  
  edited = utils::edit(meta)
  
  if (is.null(edited)) {
    stop("Metadata editing was cancelled.")
  }
  
  as.data.frame(edited, stringsAsFactors = FALSE)
}

validate_metadata = function(meta) {
  required_cols = c("sample", "cell_line", "treatment", "treatment_type", "rep", "timepoint")
  missing_cols = setdiff(required_cols, colnames(meta))
  
  if (length(missing_cols) > 0) {
    stop(
      "Metadata is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  required_nonempty = c("sample", "cell_line", "treatment", "treatment_type")
  
  empty_required = vapply(
    required_nonempty,
    function(col) any(is.na(meta[[col]]) | trimws(as.character(meta[[col]])) == ""),
    logical(1)
  )
  
  if (any(empty_required)) {
    stop(
      "Metadata has empty values in required columns: ",
      paste(required_nonempty[empty_required], collapse = ", ")
    )
  }
  
  invisible(TRUE)
}
choose_output_dir = function(
    prompt_text = "Enter full path to output directory: ",
    caption = "Select output directory",
    create = TRUE
) {
  repeat {
    selected = NULL
    
    if (.Platform$OS.type == "windows") {
      selected = tryCatch(
        utils::choose.dir(caption = caption),
        error = function(e) NA_character_
      )
    }
    
    if (!is.null(selected) && !is.na(selected) && nzchar(selected) && dir.exists(selected)) {
      return(normalizePath(selected, winslash = "/", mustWork = TRUE))
    }
    
    cat("Directory chooser unavailable or no folder selected.\n")
    
    if (interactive()) {
      path = trimws(readline(prompt_text))
    } else {
      cat(prompt_text)
      path = trimws(readLines("stdin", n=1))
    }
        
    path = path.expand(path)
    
    if (dir.exists(path)) {
      return(normalizePath(path, winslash = "/", mustWork = TRUE))
    }
    
    if (create) {
      if (ask_yes_no("Directory does not exist. Create it?", default = TRUE)) {
        ok = dir.create(path, recursive = TRUE, showWarnings = FALSE)
        if (ok && dir.exists(path)) {
          return(normalizePath(path, winslash = "/", mustWork = TRUE))
        }
      }
    }
    
    cat("Directory not found or could not be created. Please try again.\n")
  }
}
###### CLI #############

source("run_mcd.R")

main = function() {
  
  check_required_packages(c("tools", "utils", "robustbase", "fdrtool"))
  library(tools)
  library(utils)
  library(robustbase)
  library(fdrtool)
  
  
  # script_dir = get_script_dir()
  
  show_header("CRISPR Drug Screen MCD UI")
  cat("Welcome. ")
  cat("This tool will guide you through loading counts data,\n")
  cat("reviewing metadata, and running an MCD analysis on CRISPR data.\n\n")
  cat("R version: ", R.version.string, "\n", sep = "")
  cat("OS: ", Sys.info()[["sysname"]], "\n", sep = "")
  
  input_type = choose_input_type()
  input_src = collect_input_paths()

  parsed_input = tryCatch(
    parse_input_set(input_type, input_src$paths),
    error = function(e) {
      show_header("Input Parsing Error")
      cat(conditionMessage(e), "\n")
      pause("Press Enter to exit...")
      quit(status = 1)
    }
  )

  
  if (input_type != "raw_counts") {
    show_standardized_matrix_summary(parsed_input)
  }

  meta = autodetect_metadata_from_samples(parsed_input$sample_names)
  
  show_header("Autodetect Metadata")
  cat("Metadata inferred from sample or condition names:\n\n")
  print(meta)
  cat("\n")
  
  if (ask_yes_no("Would you like to edit the metadata table now?", default = TRUE)) {
    meta = tryCatch(
      edit_metadata_interactively(meta),
      error = function(e) {
        cat("Error during metadata editing:\n")
        cat(conditionMessage(e), "\n")
        pause("Press Enter to exit...")
        quit(status = 1)
      }
    )
  }
  
  tryCatch(
    validate_metadata(meta),
    error = function(e) {
      show_header("Metadata Validation Error")
      cat(conditionMessage(e), "\n")
      pause("Press Enter to exit...")
      quit(status = 1)
    }
  )
  
  if (input_type == "raw_counts") {
  
    raw_opts = choose_raw_count_processing_options()
  
    final_input = tryCatch(
      build_gene_score_matrix_from_raw_counts(
        parsed_raw = parsed_input,
        meta = meta,
        aggregation_fun = raw_opts$aggregation_fun,
        normalize_method = raw_opts$normalize_method,
        make_z = raw_opts$make_z
      ),
      error = function(e) {
        show_header("Raw Count Processing Error")
        cat(conditionMessage(e), "\n")
        pause("Press Enter to exit...")
        quit(status = 1)
      }
    )
  
  } else {
  
    final_input = parsed_input
  
    final_input$metadata = tryCatch(
      align_metadata_to_samples(final_input$sample_names, meta),
      error = function(e) {
        show_header("Metadata Alignment Error")
        cat(conditionMessage(e), "\n")
        pause("Press Enter to exit...")
        quit(status = 1)
      }
    )
  }
  
  show_final_score_matrix_summary(final_input)
  
  out_dir = choose_output_dir(
    prompt_text = "Enter full path to output directory: ",
    caption = "Select output directory",
    create = TRUE
  )

  show_header("Review Analysis Settings")
  cat("Input type: ", final_input$input_type, "\n\n", sep = "")

  cat("Input path(s): \n")
  cat(paste0(" - ", input_src$paths), sep = "\n")
  cat("\n\n")

  cat("Output directory:\n", out_dir, "\n\n", sep = "")

  cat("Samples / conditions:\n")
  cat(paste0(" - ", final_input$sample_names), sep = "\n")
  cat("\n")

  if (!ask_yes_no("Proceed with MCD analysis?", default = TRUE)) {
    cat("Operation cancelled by user.\n")
    pause("Press Enter to exit...")
    quit(status = 0)
  }

  show_header("Running MCD Analysis")
  # cat("This is where your MCD workflow will run.\n")
  cat("Rows (genes): ", nrow(final_input$score_matrix), "\n", sep = "")
  cat("Conditions: ", ncol(final_input$score_matrix), "\n\n", sep = "")

  res = run_mcd_pipeline(
    final_input = final_input,
    out_dir = out_dir)
  
  # add excel file
  
  
  
  write.csv(final_input$metadata, file.path(out_dir, "edited_metadata.csv"), row.names = FALSE)
  write.csv(final_input$score_matrix, file.path(out_dir, "gene_score_matrix.csv"), row.names = TRUE)
  write.csv(res$out_df, file.path(out_dir, "mcd_res.csv"), row.names = TRUE)
  
  
  # 
  # show_header("Complete")
  # cat("Analysis setup completed successfully.\n")
  # cat("Saved files to:\n")
  # cat(" - ", file.path(out_dir, "edited_metadata.csv"), "\n", sep = "")
  # cat(" - ", file.path(out_dir, "gene_score_matrix.csv"), "\n", sep = "")
  # cat("\n")
  # 
  # pause("Press Enter to exit...")
}

main()