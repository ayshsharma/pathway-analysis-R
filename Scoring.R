# =========================================================
# 1) Helper: normalize KEGG compound IDs
# =========================================================
normalize_compound_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^cpd:", "", x, ignore.case = TRUE)
  x <- sub("^compound:", "", x, ignore.case = TRUE)
  x <- toupper(x)
  x <- x[!is.na(x) & x != ""]
  unique(x)
}


# =========================================================
# 2) Helper: normalize KEGG pathway IDs
#    Examples normalized to "hsa00400":
#    - "hsa00400"
#    - "map00400"
#    - "path:hsa00400"
#    - "00400"
# =========================================================
normalize_pathway_ids <- function(x, organism_code) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^path:", "", x, ignore.case = TRUE)
  x <- tolower(x)
  
  # Extract trailing 5 digits if present
  digits <- sub(".*?(\\d{5})$", "\\1", x)
  
  # Keep only valid 5-digit pathway codes
  valid <- grepl("^\\d{5}$", digits)
  digits[!valid] <- NA_character_
  
  out <- ifelse(is.na(digits), NA_character_, paste0(tolower(organism_code), digits))
  out
}


# =========================================================
# 3) Helper: build normalized pathway -> compound map
#    from kegg_pathway_nested_maps[[organism_code]]
# =========================================================
build_normalized_organism_map <- function(kegg_pathway_nested_maps, organism_code) {
  if (!organism_code %in% names(kegg_pathway_nested_maps)) {
    return(NULL)
  }
  
  organism_map_raw <- kegg_pathway_nested_maps[[organism_code]]
  raw_names <- names(organism_map_raw)
  norm_names <- normalize_pathway_ids(raw_names, organism_code = organism_code)
  
  out <- list()
  
  for (i in seq_along(organism_map_raw)) {
    norm_id <- norm_names[i]
    if (is.na(norm_id) || norm_id == "") next
    
    entry <- organism_map_raw[[i]]
    
    # Try to extract compound_ids safely
    pathway_compounds <- character(0)
    
    if (is.list(entry) && "compound_ids" %in% names(entry)) {
      pathway_compounds <- entry$compound_ids
    } else if (is.character(entry)) {
      # fallback if entry itself is already a character vector
      pathway_compounds <- entry
    }
    
    pathway_compounds <- normalize_compound_ids(pathway_compounds)
    
    if (is.null(out[[norm_id]])) {
      out[[norm_id]] <- pathway_compounds
    } else {
      out[[norm_id]] <- unique(c(out[[norm_id]], pathway_compounds))
    }
  }
  
  out
}


# =========================================================
# 4) Helper: safely pull the compound table to score
#    Tries compounds_key first, then optional fallbacks
# =========================================================
get_compound_table <- function(validation_main_db, pathway, frac_key, compounds_key,
                               fallback_keys = c("validation", "10", "sample")) {
  
  if (is.null(validation_main_db[[pathway]]) ||
      is.null(validation_main_db[[pathway]][[frac_key]])) {
    return(NULL)
  }
  
  frac_entry <- validation_main_db[[pathway]][[frac_key]]
  
  # Try requested key first
  if (!is.null(compounds_key) &&
      compounds_key %in% names(frac_entry) &&
      is.data.frame(frac_entry[[compounds_key]]) &&
      "kegg_id" %in% names(frac_entry[[compounds_key]])) {
    return(frac_entry[[compounds_key]])
  }
  
  # Try fallback keys
  for (k in fallback_keys) {
    if (k %in% names(frac_entry) &&
        is.data.frame(frac_entry[[k]]) &&
        "kegg_id" %in% names(frac_entry[[k]])) {
      return(frac_entry[[k]])
    }
  }
  
  NULL
}


# =========================================================
# 5) Score one pathway + one fraction
# =========================================================
score_validation_compounds_corrected <- function(
    pathway,
    frac,
    validation_main_db,
    all_validation_results,
    organism_code,
    kegg_pathway_nested_maps,
    compounds_key = "validation",   # can also be "10" or "sample"
    truth_frac = "1",
    truth_key = "sample",
    hits_col = "Hits",
    total_col = "Total",
    verbose = FALSE, 
    degree = 1
) {
  
  frac_key <- as.character(frac)
  truth_frac_key <- as.character(truth_frac)
  
  # -----------------------------------------------------
  # 1) Check results entry exists
  # -----------------------------------------------------
  if (is.null(all_validation_results[[pathway]]) ||
      is.null(all_validation_results[[pathway]][[frac_key]])) {
    return(list(
      status = "missing_results_entry",
      reason = "Missing pathway/fraction in all_validation_results",
      scores = NULL
    ))
  }
  
  result_entry <- all_validation_results[[pathway]][[frac_key]]
  
  # Respect analysis status
  if (is.null(result_entry$status) || result_entry$status != "ok") {
    return(list(
      status = ifelse(is.null(result_entry$status), "not_ok", result_entry$status),
      reason = if (!is.null(result_entry$reason)) result_entry$reason else "No valid pathway analysis available",
      scores = NULL
    ))
  }
  
  if (is.null(result_entry$pathway_anal)) {
    return(list(
      status = "missing_pathway_anal",
      reason = "status == 'ok' but pathway_anal is missing",
      scores = NULL
    ))
  }
  
  # -----------------------------------------------------
  # 2) Build normalized KEGG pathway membership map
  # -----------------------------------------------------
  organism_map <- build_normalized_organism_map(
    kegg_pathway_nested_maps = kegg_pathway_nested_maps,
    organism_code = organism_code
  )
  
  if (is.null(organism_map)) {
    return(list(
      status = "missing_organism_map",
      reason = paste("No KEGG map found for organism:", organism_code),
      scores = NULL
    ))
  }
  
  # -----------------------------------------------------
  # 3) Pull compound table to score
  # -----------------------------------------------------
  compound_table <- get_compound_table(
    validation_main_db = validation_main_db,
    pathway = pathway,
    frac_key = frac_key,
    compounds_key = compounds_key
  )
  
  if (is.null(compound_table)) {
    return(list(
      status = "missing_compounds",
      reason = paste("Could not find compound table for pathway =", pathway,
                     "frac =", frac_key, "using compounds_key =", compounds_key),
      scores = NULL
    ))
  }
  
  compounds_norm <- normalize_compound_ids(compound_table$kegg_id)
  
  if (length(compounds_norm) == 0) {
    return(list(
      status = "empty_compounds",
      reason = "No compounds available after normalization",
      scores = data.frame(
        compound_id = character(0),
        compound_score = numeric(0),
        matched_pathway_count = integer(0),
        is_true_source_pathway_compound = logical(0),
        matched_pathways = character(0),
        stringsAsFactors = FALSE
      )
    ))
  }
  
  # Creating a mass lookup table
  
  # Build compound -> mass lookup safely (preserve row alignment first)
  compound_ids_rowwise <- normalize_compound_ids_rowwise(compound_table$kegg_id)
  compound_masses_rowwise <- suppressWarnings(as.numeric(compound_table$mol_mass))
  
  mass_lookup_df <- data.frame(
    compound_id = compound_ids_rowwise,
    mol_mass = compound_masses_rowwise,
    stringsAsFactors = FALSE
  )
  
  # Remove blank/NA IDs
  mass_lookup_df <- mass_lookup_df[
    !is.na(mass_lookup_df$compound_id) &
      mass_lookup_df$compound_id != "",
    ,
    drop = FALSE
  ]
  
  # For each compound_id, keep the first finite mass if available
  mass_lookup_df <- mass_lookup_df[order(!is.finite(mass_lookup_df$mol_mass)), , drop = FALSE]
  mass_lookup_df <- mass_lookup_df[!duplicated(mass_lookup_df$compound_id), , drop = FALSE]
  
  mass_lookup <- mass_lookup_df$mol_mass
  names(mass_lookup) <- mass_lookup_df$compound_id
  
  
  # -----------------------------------------------------
  # 4) Pull truth set from frac = 1 sample
  # -----------------------------------------------------
  truth_ids <- character(0)
  
  if (!is.null(validation_main_db[[pathway]]) &&
      !is.null(validation_main_db[[pathway]][[truth_frac_key]]) &&
      !is.null(validation_main_db[[pathway]][[truth_frac_key]][[truth_key]]) &&
      is.data.frame(validation_main_db[[pathway]][[truth_frac_key]][[truth_key]]) &&
      "kegg_id" %in% names(validation_main_db[[pathway]][[truth_frac_key]][[truth_key]])) {
    
    truth_ids <- normalize_compound_ids(
      validation_main_db[[pathway]][[truth_frac_key]][[truth_key]]$kegg_id
    )
  }
  
  # -----------------------------------------------------
  # 5) Read pathway_anal correctly (matrix OR data.frame)
  # -----------------------------------------------------
  pathway_results_raw <- result_entry$pathway_anal
  
  if (!(is.matrix(pathway_results_raw) || is.data.frame(pathway_results_raw))) {
    return(list(
      status = "bad_pathway_table",
      reason = "pathway_anal is neither a matrix nor a data.frame",
      scores = NULL
    ))
  }
  
  # Convert to data.frame for easier named-column handling
  pathway_results <- as.data.frame(pathway_results_raw, stringsAsFactors = FALSE)
  
  # Preserve row names (these are pathway IDs in your object)
  pathway_ids_raw <- rownames(pathway_results_raw)
  
  if (is.null(pathway_ids_raw) || length(pathway_ids_raw) != nrow(pathway_results)) {
    return(list(
      status = "bad_pathway_table",
      reason = "Pathway IDs are not available in rownames(pathway_anal)",
      scores = NULL
    ))
  }
  
  if (!hits_col %in% colnames(pathway_results)) {
    return(list(
      status = "bad_pathway_table",
      reason = paste("Missing column:", hits_col),
      scores = NULL
    ))
  }
  
  if (!total_col %in% colnames(pathway_results)) {
    return(list(
      status = "bad_pathway_table",
      reason = paste("Missing column:", total_col),
      scores = NULL
    ))
  }
  
  pathway_ids_norm <- normalize_pathway_ids(pathway_ids_raw, organism_code = organism_code)
  
  # -----------------------------------------------------
  # 6) Initialize score table
  # -----------------------------------------------------
  matched_pathways_list <- vector("list", length(compounds_norm))
  
  score_df <- data.frame(
    compound_id = compounds_norm,
    compound_score = 0,
    matched_pathway_count = 0L,
    is_true_source_pathway_compound = compounds_norm %in% truth_ids,
    stringsAsFactors = FALSE
  )
  
  n_result_rows <- nrow(pathway_results)
  n_rows_used <- 0L
  n_total_compound_matches <- 0L
  
  # -----------------------------------------------------
  # 7) Scoring loop
  #    score += Hits / Total for each matched pathway
  # -----------------------------------------------------
  for (j in seq_len(nrow(pathway_results))) {
    
    pathway_id_j <- pathway_ids_norm[j]
    hits_j <- suppressWarnings(as.numeric(pathway_results[j, hits_col]))
    total_j <- suppressWarnings(as.numeric(pathway_results[j, total_col]))
    
    if (is.na(pathway_id_j) || pathway_id_j == "") next
    if (is.na(hits_j) || is.na(total_j) || total_j == 0) next
    if (!pathway_id_j %in% names(organism_map)) next
    
    pathway_compounds <- organism_map[[pathway_id_j]]
    pathway_compounds <- normalize_compound_ids(pathway_compounds)
    
    matched_compounds <- intersect(score_df$compound_id, pathway_compounds)
    if (length(matched_compounds) == 0) next
    
    contribution <- (hits_j / total_j)^degree
    
    # Skip pathological cases
    if (!is.finite(contribution) || length(contribution) != 1) next
    
    # Group matched compounds by near-equal mass within this pathway
    mass_clusters <- cluster_compounds_by_mass(
      compounds = matched_compounds,
      mass_lookup = mass_lookup,
      tolerance_ppm = 10
    )
    
    if (length(mass_clusters) == 0) next
    
    for (cluster in mass_clusters) {
      
      cluster <- unique(as.character(cluster))
      cluster <- cluster[!is.na(cluster) & cluster != ""]
      
      if (length(cluster) == 0) next
      
      idx <- which(score_df$compound_id %in% cluster)
      
      # If none of these compounds are present in score_df, skip safely
      if (length(idx) == 0) next
      
      # # split one pathway contribution across redundant candidates
      # cluster_share <- contribution / length(cluster)
      
      
      # if you want NO contribution if redundant compound present
      cluster_share <- contribution

      if (length(cluster) > 1)
        cluster_share = 0
      
      # Extra safety
      if (!is.finite(cluster_share) || length(cluster_share) != 1) next
      
      score_df$compound_score[idx] <- score_df$compound_score[idx] + cluster_share
      score_df$matched_pathway_count[idx] <- score_df$matched_pathway_count[idx] + 1L
      
      for (ii in idx) {
        matched_pathways_list[[ii]] <- unique(c(matched_pathways_list[[ii]], pathway_id_j))
      }
    }
    
    n_rows_used <- n_rows_used + 1L
    n_total_compound_matches <- n_total_compound_matches + length(matched_compounds)
    
    for (ii in idx) {
      matched_pathways_list[[ii]] <- unique(c(matched_pathways_list[[ii]], pathway_id_j))
    }
    
    n_rows_used <- n_rows_used + 1L
    n_total_compound_matches <- n_total_compound_matches + length(matched_compounds)
  }
  
  # Collapse matched pathways into a readable string
  score_df$matched_pathways <- vapply(
    matched_pathways_list,
    function(x) {
      if (length(x) == 0) "" else paste(sort(unique(x)), collapse = ";")
    },
    character(1)
  )
  
  # Rank compounds
  score_df <- score_df[
    order(-score_df$compound_score,
          -score_df$matched_pathway_count,
          score_df$compound_id),
    ,
    drop = FALSE
  ]
  
  rownames(score_df) <- NULL
  
  # -----------------------------------------------------
  # 8) Optional diagnostics
  # -----------------------------------------------------
  if (verbose) {
    message("Pathway: ", pathway, " | frac: ", frac_key)
    message("Compounds scored: ", nrow(score_df))
    message("Pathway result rows: ", n_result_rows)
    message("Pathway rows used after normalization: ", n_rows_used)
    message("Total compound-pathway matches: ", n_total_compound_matches)
    message("Compounds with non-zero score: ", sum(score_df$compound_score > 0))
    message("True source-pathway compounds found in score table: ",
            sum(score_df$is_true_source_pathway_compound))
  }
  
  list(
    status = "ok",
    reason = NULL,
    n_input_compounds = length(compounds_norm),
    n_rows_in_pathway_results = n_result_rows,
    n_rows_used = n_rows_used,
    n_total_compound_matches = n_total_compound_matches,
    scores = score_df
  )
}


# =========================================================
# 6) Wrapper: score all pathway/fraction combinations
# =========================================================
score_all_validation_compounds <- function(
    validation_main_db,
    all_validation_results,
    organism_code,
    kegg_pathway_nested_maps,
    compounds_key = "validation",   # can also be "10" or "sample"
    truth_frac = "1",
    truth_key = "sample",
    hits_col = "Hits",
    total_col = "Total",
    verbose = FALSE, 
    degree = 1
) {
  
  out <- list()
  
  pathways <- intersect(names(validation_main_db), names(all_validation_results))
  pathways <- pathways[!is.na(pathways) & pathways != ""]
  
  for (pathway in pathways) {
    if (is.null(validation_main_db[[pathway]]) ||
        is.null(all_validation_results[[pathway]]) ||
        !is.list(validation_main_db[[pathway]]) ||
        !is.list(all_validation_results[[pathway]])) {
      next
    }
    
    out[[pathway]] <- list()
    
    frac_keys <- intersect(
      names(validation_main_db[[pathway]]),
      names(all_validation_results[[pathway]])
    )
    frac_keys <- frac_keys[!is.na(frac_keys) & frac_keys != ""]
    
    for (frac_key in frac_keys) {
      out[[pathway]][[frac_key]] <- score_validation_compounds_corrected(
        pathway = pathway,
        frac = frac_key,
        validation_main_db = validation_main_db,
        all_validation_results = all_validation_results,
        organism_code = organism_code,
        kegg_pathway_nested_maps = kegg_pathway_nested_maps,
        compounds_key = compounds_key,
        truth_frac = truth_frac,
        truth_key = truth_key,
        hits_col = hits_col,
        total_col = total_col,
        verbose = verbose,
        degree = degree
      )
    }
  }
  
  out
}




normalize_compound_ids_rowwise <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^cpd:", "", x, ignore.case = TRUE)
  x <- sub("^compound:", "", x, ignore.case = TRUE)
  x <- toupper(x)
  x
}










evaluate_compound_rankings <- function(
    compound_scores_db,
    validation_main_db = NULL,
    all_validation_results = NULL,
    top_n = 10,
    k_mode = c("fixed", "sample_size", "analyzed_size")
) {
  
  k_mode <- match.arg(k_mode)
  
  summary_rows <- list()
  rank_rows <- list()
  
  for (pathway in names(compound_scores_db)) {
    if (is.null(compound_scores_db[[pathway]]) || !is.list(compound_scores_db[[pathway]])) next
    
    for (frac_key in names(compound_scores_db[[pathway]])) {
      entry <- compound_scores_db[[pathway]][[frac_key]]
      
      if (is.null(entry) || !is.list(entry)) next
      if (is.null(entry$status) || entry$status != "ok") next
      if (is.null(entry$scores) || !is.data.frame(entry$scores)) next
      
      df <- entry$scores
      
      needed_cols <- c("compound_id", "compound_score", "is_true_source_pathway_compound")
      if (!all(needed_cols %in% names(df))) next
      if (nrow(df) == 0) next
      
      # Ensure sorted in descending score order
      extra_sort_col <- if ("matched_pathway_count" %in% names(df)) df$matched_pathway_count else rep(0, nrow(df))
      df <- df[order(-df$compound_score, -extra_sort_col, df$compound_id), , drop = FALSE]
      rownames(df) <- NULL
      
      # Add ranks
      df$rank <- seq_len(nrow(df))
      n_total <- nrow(df)
      truth <- as.logical(df$is_true_source_pathway_compound)
      n_true <- sum(truth, na.rm = TRUE)
      
      if (n_true == 0) next
      
      # -----------------------------------------
      # Choose k dynamically
      # -----------------------------------------
      if (k_mode == "fixed") {
        k <- min(top_n, n_total)
        
      } else if (k_mode == "sample_size") {
        if (is.null(validation_main_db) ||
            is.null(validation_main_db[[pathway]]) ||
            is.null(validation_main_db[[pathway]][[frac_key]]) ||
            is.null(validation_main_db[[pathway]][[frac_key]][["sample"]]) ||
            !is.data.frame(validation_main_db[[pathway]][[frac_key]][["sample"]])) {
          next
        }
        
        k <- nrow(validation_main_db[[pathway]][[frac_key]][["sample"]])
        k <- min(k, n_total)
        
      } else if (k_mode == "analyzed_size") {
        if (is.null(all_validation_results) ||
            is.null(all_validation_results[[pathway]]) ||
            is.null(all_validation_results[[pathway]][[frac_key]]) ||
            is.null(all_validation_results[[pathway]][[frac_key]]$n_compounds)) {
          next
        }
        
        k <- as.integer(all_validation_results[[pathway]][[frac_key]]$n_compounds)
        k <- min(k, n_total)
      }
      
      if (is.na(k) || k < 1) next
      
      topk_truth <- truth[seq_len(k)]
      
      precision_at_k <- mean(topk_truth)
      recall_at_k <- sum(topk_truth) / n_true
      top1_is_true <- truth[1]
      
      # Random baseline = prevalence of true compounds in this ranked set
      random_baseline <- n_true / n_total
      enrichment_over_random <- if (random_baseline > 0) precision_at_k / random_baseline else NA_real_
      
      # Average precision over full ranking
      rel <- as.integer(truth)
      precision_at_i <- cumsum(rel) / seq_along(rel)
      average_precision <- sum(precision_at_i * rel) / sum(rel)
      
      # Ranks of true compounds
      true_ranks <- df$rank[truth]
      mean_rank_true <- mean(true_ranks)
      median_rank_true <- median(true_ranks)
      best_rank_true <- min(true_ranks)
      
      # Rank percentiles
      if (n_total == 1) {
        rank_percentile <- 0
      } else {
        rank_percentile <- (df$rank - 1) / (n_total - 1)
      }
      
      mean_rank_pct_true <- mean(rank_percentile[truth])
      median_rank_pct_true <- median(rank_percentile[truth])
      
      # Score separation
      mean_score_true <- mean(df$compound_score[truth])
      mean_score_false <- mean(df$compound_score[!truth])
      score_gap <- mean_score_true - mean_score_false
      
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        pathway = pathway,
        frac = as.numeric(frac_key),
        frac_key = frac_key,
        n_ranked = n_total,
        n_true = n_true,
        k = k,
        k_mode = k_mode,
        top1_is_true = top1_is_true,
        precision_at_k = precision_at_k,
        recall_at_k = recall_at_k,
        average_precision = average_precision,
        random_baseline = random_baseline,
        enrichment_over_random = enrichment_over_random,
        best_rank_true = best_rank_true,
        mean_rank_true = mean_rank_true,
        median_rank_true = median_rank_true,
        mean_rank_pct_true = mean_rank_pct_true,
        median_rank_pct_true = median_rank_pct_true,
        mean_score_true = mean_score_true,
        mean_score_false = mean_score_false,
        score_gap = score_gap,
        stringsAsFactors = FALSE
      )
      
      rank_rows[[length(rank_rows) + 1]] <- data.frame(
        pathway = pathway,
        frac = as.numeric(frac_key),
        frac_key = frac_key,
        compound_id = df$compound_id,
        compound_score = df$compound_score,
        rank = df$rank,
        rank_percentile = rank_percentile,
        is_true = truth,
        k = k,
        k_mode = k_mode,
        in_top_k = df$rank <= k,
        stringsAsFactors = FALSE
      )
    }
  }
  
  summary_df <- if (length(summary_rows) > 0) do.call(rbind, summary_rows) else data.frame()
  rank_df <- if (length(rank_rows) > 0) do.call(rbind, rank_rows) else data.frame()
  
  list(
    summary = summary_df,
    ranks = rank_df
  )
}




cluster_compounds_by_mass <- function(compounds, mass_lookup, tolerance_ppm = 10) {
  
  compounds <- unique(as.character(compounds))
  if (length(compounds) == 0) return(list())
  
  masses <- mass_lookup[compounds]
  
  # Keep valid masses
  valid <- !is.na(masses) & is.finite(masses)
  valid_compounds <- compounds[valid]
  valid_masses <- as.numeric(masses[valid])
  
  invalid_compounds <- compounds[!valid]
  
  # No valid masses -> every invalid one becomes its own singleton cluster
  if (length(valid_compounds) == 0) {
    return(as.list(invalid_compounds))
  }
  
  # Pairwise ppm distance using average mass in denominator
  ppm_dist <- outer(
    valid_masses,
    valid_masses,
    function(a, b) abs(a - b) / ((a + b) / 2) * 1e6
  )
  
  # Adjacency matrix: TRUE if masses are within tolerance
  adj <- ppm_dist <= tolerance_ppm
  diag(adj) <- TRUE
  
  # Find connected components with DFS
  n <- length(valid_compounds)
  visited <- rep(FALSE, n)
  clusters <- list()
  
  for (i in seq_len(n)) {
    if (visited[i]) next
    
    stack <- i
    component <- integer(0)
    
    while (length(stack) > 0) {
      node <- stack[[1]]
      stack <- stack[-1]
      
      if (visited[node]) next
      visited[node] <- TRUE
      component <- c(component, node)
      
      neighbors <- which(adj[node, ] & !visited)
      if (length(neighbors) > 0) {
        stack <- c(neighbors, stack)
      }
    }
    
    clusters[[length(clusters) + 1]] <- valid_compounds[component]
  }
  
  # Invalid masses are treated as singleton clusters
  if (length(invalid_compounds) > 0) {
    for (cmp in invalid_compounds) {
      clusters[[length(clusters) + 1]] <- cmp
    }
  }
  
  clusters
}





# Running the functions

compound_scores_squared_corrected_db <- score_all_validation_compounds(
  validation_main_db = validation_main_db,
  all_validation_results = all_validation_results,
  organism_code = "hsa",
  kegg_pathway_nested_maps = kegg_pathway_nested_maps,
  compounds_key = "validation",
  truth_frac = "1",
  truth_key = "sample",
  verbose = TRUE,
  degree = 3
)



# Check if it worked

eval_results_squared_corrected <- evaluate_compound_rankings(
  compound_scores_db = compound_scores_cubed_corrected_db,
  validation_main_db = validation_main_db,
  k_mode = "sample_size"
)




# Plotting



library(ggplot2)
library(dplyr)


# Comparision of corrected precision@k to uncorrected and random baseline
ggplot() +
  # geom_point(aes(y = precision_at_k), alpha = 0.5, color = "steelblue") +
  # geom_smooth(data = subset(eval_results_squared$summary, frac < 1), aes(x = frac, y = precision_at_k), se = FALSE, color = "orange") +
  geom_smooth(data = subset(eval_results_cubed_split$summary, frac < 1), aes(x = frac, y = precision_at_k), se = FALSE, color = "red") +
  geom_smooth(data = subset(eval_results_cubed_corrected$summary, frac < 1), aes(x = frac, y = precision_at_k), se = FALSE, color = "blue") +
  # geom_smooth(data = subset(eval_results$summary, frac < 1), aes(x = frac, y = precision_at_k), se = FALSE, color = "red") +
  # geom_smooth(data = subset(eval_results_split$summary, frac < 1), aes(x = frac, y = precision_at_k), se = FALSE, color = "green") +
  # geom_smooth(data = subset(eval_results_corrected$summary, frac < 1), aes(x = frac, y = precision_at_k), se = FALSE, color = "steelblue") +
  geom_smooth(data = subset(eval_results_corrected$summary, frac < 1), aes(x = frac, y = random_baseline), se = FALSE, color = "black", linetype = "dashed") +
  labs(
    title = "Precision@k vs random baseline (Corrected vs. Split vs. Uncorrected)",
    x = "Sampling fraction",
    y = "Value"
  )



# Degree = 1


# Precision@k across sampling fractions
ggplot(
  subset(eval_results_split$summary, frac < 1),
  aes(x = frac, y = precision_at_k)
) +
  # geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE) +
  labs(
    title = "Precision@k across sampling fractions (excluding frac = 1) (Split)",
    x = "Sampling fraction",
    y = "Precision@k"
  )



# library(dplyr)
# library(ggplot2)

# Mean precision@k across sampling fractions
plot_df <- eval_results_cubed_corrected$summary %>%
  filter(frac < 1) %>%
  group_by(frac) %>%
  summarise(
    mean_precision = mean(precision_at_k, na.rm = TRUE),
    sd_precision = sd(precision_at_k, na.rm = TRUE),
    n = n(),
    se_precision = sd_precision / sqrt(n),
    .groups = "drop"
  )

ggplot(plot_df, aes(x = frac, y = mean_precision)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = mean_precision - se_precision,
      ymax = mean_precision + se_precision
    ),
    width = 0.005
  ) +
  geom_errorbar(
    aes(
      ymin = mean_precision - 2 * se_precision,
      ymax = mean_precision + 2 * se_precision
    ),
    width = 0.002
  ) +
  labs(
    title = "Mean Precision@k across sampling fractions (excluding frac = 1) (Split)",
    x = "Sampling fraction",
    y = "Mean Precision@k"
  ) +
  theme_minimal()






# Precision@k across sampling fractions
ggplot(
  subset(eval_results_split$summary, frac < 1),
  aes(x = frac, y = precision_at_k)
) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE) +
  scale_x_continuous(breaks = seq(0, 0.2, by = 0.02)) +
  labs(
    title = "Precision@k across sampling fractions (excluding frac = 1) (Split)",
    x = "Sampling fraction",
    y = "Precision@k"
  ) +
  theme_minimal()







# Comparision of precision@k to random baseline
ggplot(plot_df, aes(x = frac)) +
  
  # geom_point(aes(y = precision_at_k), alpha = 0.5, color = "steelblue") +
  
  
  # geom_line(aes(y = mean_precision), color = "green") +
  
  geom_smooth(aes(y = mean_precision), se = FALSE, color = "green") +
  
  geom_smooth(aes(y = mean_precision + sd_precision),
              se = FALSE, linetype = "dashed", color = "blue") +
  
  geom_smooth(aes(y = mean_precision - sd_precision),
              se = FALSE, linetype = "dashed", color = "blue") +
  
  geom_smooth(aes(y = mean_precision + 2*sd_precision),
              se = FALSE, linetype = "dashed", color = "red") +
  
  geom_smooth(aes(y = mean_precision - 2*sd_precision),
              se = FALSE, linetype = "dashed", color = "red") +
  
  geom_smooth(data = subset(eval_results_corrected$summary, frac < 1), 
              aes(x = frac, y = random_baseline), se = FALSE, color = "black", linetype = "dashed") +
  
  ylim(0, 1) +
  
  xlim(0, 0.2) +
  
  labs(
    title = "Precision@k vs random baseline (Degree = 3) (Corrected)",
    x = "Sampling fraction",
    y = "Value"
  )



# Percentile ranking of true vs false compounds
ggplot(eval_results_split$ranks, aes(x = is_true, y = rank_percentile)) +
  geom_boxplot() +
  scale_y_reverse() +
  labs(
    title = "Rank percentile of true vs false compounds (Corrected)",
    x = "True source-pathway compound",
    y = "Rank percentile (top = better)"
  )





# Score distribution of true vs false compounds
ggplot(subset(eval_results_split$ranks, frac < 1), aes(x = is_true, y = compound_score)) +
  geom_boxplot() +
  labs(
    title = "Compound scores for true vs false compounds (Corrected)",
    x = "True source-pathway compound",
    y = "Compound score"
  )



# 

ggplot(subset(eval_results_split$summary, frac < 1), aes(x = frac, y = average_precision)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE) +
  labs(
    title = "Average precision across sampling fractions (Corrected)",
    x = "Sampling fraction",
    y = "Average precision"
  )



plot_df <- eval_results_split$summary %>%
  filter(frac < 1) %>%
  group_by(frac) %>%
  summarise(
    mean_precision = mean(average_precision, na.rm = TRUE),
    sd_precision = sd(average_precision, na.rm = TRUE),
    n = n(),
    se_precision = sd_precision / sqrt(n),
    .groups = "drop"
  )

ggplot(plot_df, aes(x = frac, y = mean_precision)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = mean_precision - se_precision,
      ymax = mean_precision + se_precision
    ),
    width = 0.005
  ) +
  labs(
    title = "Average precision across sampling fractions (Degree = 3) (Corrected)",
    x = "Sampling fraction",
    y = "Average precision"
  ) +
  theme_minimal()








# Degree = 2

# Precision@k across sampling fractions
ggplot(
  subset(eval_results_split$summary, frac < 1),
  aes(x = frac, y = precision_at_k)
) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE) +
  labs(
    title = "Precision@k across sampling fractions (excluding frac = 1) (Corrected)",
    x = "Sampling fraction",
    y = "Precision@k"
  )



# library(dplyr)
# library(ggplot2)

# Mean precision@k across sampling fractions
plot_df <- eval_results_cubed_corrected$summary %>%
  filter(frac < 1) %>%
  group_by(frac) %>%
  summarise(
    mean_precision = mean(precision_at_k, na.rm = TRUE),
    sd_precision = sd(precision_at_k, na.rm = TRUE),
    n = n(),
    se_precision = sd_precision / sqrt(n),
    .groups = "drop"
  )

ggplot(plot_df, aes(x = frac, y = mean_precision)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = mean_precision - se_precision,
      ymax = mean_precision + se_precision
    ),
    width = 0.005
  ) +
  geom_errorbar(
    aes(
      ymin = mean_precision - 2 * se_precision,
      ymax = mean_precision + 2 * se_precision
    ),
    width = 0.002
  ) +
  labs(
    title = "Mean Precision@k across sampling fractions (excluding frac = 1) (Split)",
    x = "Sampling fraction",
    y = "Mean Precision@k"
  ) +
  theme_minimal()






# Precision@k across sampling fractions
ggplot(
  subset(eval_results_split$summary, frac < 1),
  aes(x = frac, y = precision_at_k)
) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE) +
  scale_x_continuous(breaks = seq(0, 0.2, by = 0.02)) +
  labs(
    title = "Precision@k across sampling fractions (excluding frac = 1) (Split)",
    x = "Sampling fraction",
    y = "Precision@k"
  ) +
  theme_minimal()







# Comparision of precision@k to random baseline
ggplot(plot_df, aes(x = frac)) +
  
  # geom_point(aes(y = precision_at_k), alpha = 0.5, color = "steelblue") +
  
  
  # geom_line(aes(y = mean_precision), color = "green") +
  
  geom_smooth(aes(y = mean_precision), se = FALSE, color = "green") +
  
  geom_smooth(aes(y = mean_precision + sd_precision),
              se = FALSE, linetype = "dashed", color = "blue") +
  
  geom_smooth(aes(y = mean_precision - sd_precision),
              se = FALSE, linetype = "dashed", color = "blue") +
  
  geom_smooth(aes(y = mean_precision + 2*sd_precision),
              se = FALSE, linetype = "dashed", color = "red") +
  
  geom_smooth(aes(y = mean_precision - 2*sd_precision),
              se = FALSE, linetype = "dashed", color = "red") +
  
  geom_smooth(data = subset(eval_results_corrected$summary, frac < 1), 
              aes(x = frac, y = random_baseline), se = FALSE, color = "black", linetype = "dashed") +
  
  ylim(0, 1) +
  
  xlim(0, 0.2) +
  
  labs(
    title = "Precision@k vs random baseline (Degree = 3) (Corrected)",
    x = "Sampling fraction",
    y = "Value"
  )



# Percentile ranking of true vs false compounds
ggplot(eval_results_split$ranks, aes(x = is_true, y = rank_percentile)) +
  geom_boxplot() +
  scale_y_reverse() +
  labs(
    title = "Rank percentile of true vs false compounds (Corrected)",
    x = "True source-pathway compound",
    y = "Rank percentile (top = better)"
  )





# Score distribution of true vs false compounds
ggplot(subset(eval_results_split$ranks, frac < 1), aes(x = is_true, y = compound_score)) +
  geom_boxplot() +
  labs(
    title = "Compound scores for true vs false compounds (Corrected)",
    x = "True source-pathway compound",
    y = "Compound score"
  )



# 

ggplot(subset(eval_results_split$summary, frac < 1), aes(x = frac, y = average_precision)) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE) +
  labs(
    title = "Average precision across sampling fractions (Corrected)",
    x = "Sampling fraction",
    y = "Average precision"
  )



plot_df <- eval_results_split$summary %>%
  filter(frac < 1) %>%
  group_by(frac) %>%
  summarise(
    mean_precision = mean(average_precision, na.rm = TRUE),
    sd_precision = sd(average_precision, na.rm = TRUE),
    n = n(),
    se_precision = sd_precision / sqrt(n),
    .groups = "drop"
  )

ggplot(plot_df, aes(x = frac, y = mean_precision)) +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(
    aes(
      ymin = mean_precision - se_precision,
      ymax = mean_precision + se_precision
    ),
    width = 0.005
  ) +
  labs(
    title = "Average precision across sampling fractions (Degree = 3) (Corrected)",
    x = "Sampling fraction",
    y = "Average precision"
  ) +
  theme_minimal()
  
  