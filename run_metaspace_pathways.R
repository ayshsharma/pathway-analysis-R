# Function to do metaboanalyst queries on all the datasets we have

run_metaspace_metaboanalyst_pathway <- function(
    metaspace_examples_script,
    compound_col = "kegg_compounds",
    organism_col = "kegg_organism",
    dataset_id_col = "dataset_id",
    dataset_name_col = "dataset_name",
    kegg_lib = "current",
    method = "hyperg",
    background = "rbc",
    out_rds = "metaspace_examples_script_with_pathway_results_onsample.rds",
    out_combined_csv = "metaspace_all_pathway_results_onsample.csv",
    progress_every = 10,
    verbose = TRUE
) {
  stopifnot(requireNamespace("MetaboAnalystR", quietly = TRUE))
  
  # ---- checks ----
  required_cols <- c(compound_col, organism_col, dataset_id_col)
  missing_cols <- setdiff(required_cols, colnames(metaspace_examples_script))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in metaspace_examples_script: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  if (!is.list(metaspace_examples_script[[compound_col]])) {
    stop(
      "The compound column must be a list-column, where each row contains ",
      "a character vector of KEGG compound IDs."
    )
  }
  
  total_runs <- nrow(metaspace_examples_script)
  
  if (verbose) {
    cat("Total MetaboAnalystR runs to perform:", total_runs, "\n")
  }
  
  # ---- progress helper ----
  fmt_eta <- function(done, total, t0) {
    if (done == 0L) return(NA_character_)
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    rate <- done / elapsed
    remaining <- (total - done) / rate
    sprintf("%0.1fs elapsed, ETA %0.1fs", elapsed, remaining)
  }
  
  # ---- single dataset runner ----
  run_one_dataset <- function(i) {
    
    dataset_id <- as.character(metaspace_examples_script[[dataset_id_col]][i])
    
    dataset_name <- if (dataset_name_col %in% colnames(metaspace_examples_script)) {
      as.character(metaspace_examples_script[[dataset_name_col]][i])
    } else {
      NA_character_
    }
    
    organism <- as.character(metaspace_examples_script[[organism_col]][i])
    
    compound_ids <- metaspace_examples_script[[compound_col]][[i]]
    compound_ids <- unique(compound_ids)
    compound_ids <- compound_ids[!is.na(compound_ids) & compound_ids != ""]
    compound_ids <- sort(compound_ids)
    
    n_input <- length(compound_ids)
    
    pathway_table <- data.frame()
    ok_run <- FALSE
    err <- NA_character_
    
    if (n_input == 0) {
      return(list(
        pathway_table = data.frame(),
        n_input = 0L,
        n_pathways = 0L,
        ok_run = FALSE,
        error = "Empty compound list"
      ))
    }
    
    if (is.na(organism) || organism == "") {
      return(list(
        pathway_table = data.frame(),
        n_input = n_input,
        n_pathways = 0L,
        ok_run = FALSE,
        error = "Missing KEGG organism code"
      ))
    }
    
    tryCatch({
      
      nSet <- MetaboAnalystR::InitDataObjects(
        "conc",
        "pathora",
        FALSE,
        default.dpi = 300
      )
      
      nSet <- MetaboAnalystR::Setup.MapData(nSet, compound_ids)
      nSet <- MetaboAnalystR::CrossReferencing(nSet, "kegg")
      nSet <- MetaboAnalystR::CreateMappingResultTable(nSet)
      
      nSet <- MetaboAnalystR::SetKEGG.PathLib(nSet, organism, kegg_lib)
      nSet <- MetaboAnalystR::SetMetabolomeFilter(nSet, FALSE)
      nSet <- MetaboAnalystR::CalculateOraScore(nSet, background, method)
      
      ora <- nSet$analSet$ora.mat
      
      if (!is.null(ora) && nrow(ora) > 0) {
        pathway_table <- as.data.frame(ora, stringsAsFactors = FALSE)
        
        # Preserve pathway names if they are row names
        pathway_table$pathway_name <- rownames(ora)
        rownames(pathway_table) <- NULL
        
        # Add dataset metadata to every pathway result row
        pathway_table$dataset_id <- dataset_id
        pathway_table$dataset_name <- dataset_name
        pathway_table$kegg_organism <- organism
        pathway_table$n_input_compounds <- n_input
        
        # Move metadata columns to the front
        front_cols <- c(
          "dataset_id",
          "dataset_name",
          "kegg_organism",
          "n_input_compounds",
          "pathway_name"
        )
        
        pathway_table <- pathway_table[
          c(front_cols, setdiff(colnames(pathway_table), front_cols))
        ]
      }
      
      ok_run <- TRUE
      
    }, error = function(e) {
      err <<- conditionMessage(e)
    })
    
    list(
      pathway_table = pathway_table,
      n_input = n_input,
      n_pathways = nrow(pathway_table),
      ok_run = ok_run,
      error = err
    )
  }
  
  # ---- main loop ----
  pathway_results <- vector("list", total_runs)
  n_input_compounds <- integer(total_runs)
  n_pathways_found <- integer(total_runs)
  ok_run <- logical(total_runs)
  error_msg <- character(total_runs)
  
  t0 <- Sys.time()
  
  for (i in seq_len(total_runs)) {
    
    result_i <- run_one_dataset(i)
    
    pathway_results[[i]] <- result_i$pathway_table
    n_input_compounds[i] <- result_i$n_input
    n_pathways_found[i] <- result_i$n_pathways
    ok_run[i] <- result_i$ok_run
    error_msg[i] <- ifelse(is.na(result_i$error), NA_character_, result_i$error)
    
    if (verbose && (i %% progress_every == 0L || i == total_runs)) {
      cat(sprintf("[%d/%d] %s\n", i, total_runs, fmt_eta(i, total_runs, t0)))
    }
  }
  
  # ---- save results into original dataframe ----
  metaspace_examples_script$pathway_results <- I(pathway_results)
  metaspace_examples_script$n_input_compounds <- n_input_compounds
  metaspace_examples_script$n_pathways_found <- n_pathways_found
  metaspace_examples_script$pathway_ok_run <- ok_run
  metaspace_examples_script$pathway_error <- error_msg
  
  # ---- create combined long pathway table ----
  non_empty <- lengths(pathway_results) > 0
  
  combined_pathway_results <- if (any(non_empty)) {
    do.call(rbind, pathway_results[non_empty])
  } else {
    data.frame()
  }
  
  # ---- write outputs ----
  saveRDS(metaspace_examples_script, out_rds)
  
  utils::write.csv(
    combined_pathway_results,
    out_combined_csv,
    row.names = FALSE
  )
  
  if (verbose) {
    cat("Done.\n")
    cat("Wrote dataframe with nested pathway results to:", out_rds, "\n")
    cat("Wrote combined pathway table to:", out_combined_csv, "\n")
  }
  
  return(metaspace_examples_script)
}






metaspace_examples_script <- run_metaspace_metaboanalyst_pathway(
  metaspace_examples_script,
  compound_col = "kegg_compounds",
  organism_col = "kegg_organism",
  dataset_id_col = "dataset_id",
  dataset_name_col = "dataset_name",
  method = "hyperg",
  background = "rbc",
  progress_every = 10,
  verbose = TRUE
)



# For compounds on sample
metaspace_onsample_tmp <- run_metaspace_metaboanalyst_pathway(
  metaspace_examples_script,
  compound_col = "kegg_compounds_on_sample",
  organism_col = "kegg_organism",
  dataset_id_col = "dataset_id",
  dataset_name_col = "dataset_name",
  method = "hyperg",
  background = "rbc",
  out_rds = "metaspace_examples_script_with_pathway_results_onsample_tmp.rds",
  out_combined_csv = "metaspace_all_pathway_results_onsample.csv",
  progress_every = 10,
  verbose = TRUE
)





metaspace_examples_script$pathway_results_onsample <- 
  metaspace_onsample_tmp$pathway_results

metaspace_examples_script$n_input_compounds_onsample <- 
  metaspace_onsample_tmp$n_input_compounds

metaspace_examples_script$n_pathways_found_onsample <- 
  metaspace_onsample_tmp$n_pathways_found

metaspace_examples_script$pathway_ok_run_onsample <- 
  metaspace_onsample_tmp$pathway_ok_run

metaspace_examples_script$pathway_error_onsample <- 
  metaspace_onsample_tmp$pathway_error



saveRDS(
  metaspace_examples_script,
  file = "metaspace_examples_script_with_comp")

saveRDS(sample_10,
        file = "sample 10")


# metaspace_examples_script <- readRDS(
#   "metaspace_examples_script_with_compound_scores.rds"
# )
# 


# Scoring system

normalize_pathway_name <- function(x) {
  x <- as.character(x)
  x <- sub(" - .*$", "", x)        # remove " - Homo sapiens (human)"
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- trimws(x)
  x
}











# Creating a nested format





find_first_existing_col <- function(df, candidates) {
  hits <- intersect(candidates, colnames(df))
  if (length(hits) == 0) {
    return(NA_character_)
  }
  hits[1]
}




















make_nested_map_with_keggget <- function(
    organism_code,
    pathway_ids = NULL,
    sleep_sec = 0.1,
    verbose = TRUE
) {
  
  if (is.null(pathway_ids)) {
    pathways <- KEGGREST::keggList("pathway", organism_code)
    pathway_ids <- names(pathways)
  }
  
  nested <- setNames(
    vector("list", length(pathway_ids)),
    pathway_ids
  )
  
  for (i in seq_along(pathway_ids)) {
    
    pid <- pathway_ids[i]
    
    if (verbose && (i %% 10 == 0 || i == length(pathway_ids))) {
      message("Fetching ", organism_code, " pathway ", i, " / ", length(pathway_ids), ": ", pid)
    }
    
    kg <- tryCatch(
      KEGGREST::keggGet(pid),
      error = function(e) NULL
    )
    
    if (is.null(kg) || length(kg) == 0) {
      nested[[pid]] <- list(
        pathway_id = pid,
        pathway_name = NA_character_,
        compound_ids = character(0),
        compounds = data.frame(
          compound_id = character(0),
          compound_name = character(0),
          stringsAsFactors = FALSE
        )
      )
      next
    }
    
    entry <- kg[[1]]
    
    pathway_name <- if (!is.null(entry$NAME)) {
      entry$NAME[1]
    } else {
      NA_character_
    }
    
    if (!is.null(entry$COMPOUND)) {
      compound_ids <- names(entry$COMPOUND)
      compound_names <- unname(entry$COMPOUND)
      
      compounds_df <- data.frame(
        compound_id = compound_ids,
        compound_name = compound_names,
        stringsAsFactors = FALSE
      )
      
      compounds_df <- compounds_df[
        grepl("^C\\d{5}$", compounds_df$compound_id),
        ,
        drop = FALSE
      ]
      
      compound_ids <- sort(unique(compounds_df$compound_id))
      
    } else {
      compounds_df <- data.frame(
        compound_id = character(0),
        compound_name = character(0),
        stringsAsFactors = FALSE
      )
      
      compound_ids <- character(0)
    }
    
    nested[[pid]] <- list(
      pathway_id = pid,
      pathway_name = pathway_name,
      compound_ids = compound_ids,
      compounds = compounds_df
    )
    
    Sys.sleep(sleep_sec)
  }
  
  nested
}


organisms_needed <- unique(na.omit(metaspace_examples_script$kegg_organism))
organisms_needed <- organisms_needed[organisms_needed != ""]

kegg_pathway_nested_maps <- setNames(
  vector("list", length(organisms_needed)),
  organisms_needed
)

for (org in organisms_needed) {
  kegg_pathway_nested_maps[[org]] <- make_nested_map_with_keggget(
    organism_code = org,
    sleep_sec = 0.1,
    verbose = TRUE
  )
}



saveRDS(
  kegg_pathway_nested_maps,
  file = "kegg_pathway_nested_maps.rds"
)
# kegg_pathway_nested_maps <- readRDS(file = "kegg_pathway_nested_maps.rds")



score_compounds_using_nested_kegg_maps <- function(
    compound_ids,
    pathway_results,
    organism_code,
    kegg_pathway_nested_maps,
    hits_col = "Hits",
    total_col = "Total",
    pathway_col = "pathway_name",
    verbose = FALSE
) {
  
  compound_ids <- sort(unique(na.omit(compound_ids)))
  compound_ids <- compound_ids[compound_ids != ""]
  
  score_df <- data.frame(
    compound_id = compound_ids,
    compound_score = rep(0, length(compound_ids)),
    stringsAsFactors = FALSE
  )
  
  if (length(compound_ids) == 0) {
    return(score_df)
  }
  
  if (
    is.null(pathway_results) ||
    !is.data.frame(pathway_results) ||
    nrow(pathway_results) == 0
  ) {
    return(score_df)
  }
  
  if (!organism_code %in% names(kegg_pathway_nested_maps)) {
    warning("No nested KEGG pathway map found for organism: ", organism_code)
    return(score_df)
  }
  
  organism_map <- kegg_pathway_nested_maps[[organism_code]]
  
  if (!hits_col %in% colnames(pathway_results)) {
    stop("Could not find Hits column: ", hits_col)
  }
  
  if (!total_col %in% colnames(pathway_results)) {
    stop("Could not find Total column: ", total_col)
  }
  
  if (!pathway_col %in% colnames(pathway_results)) {
    stop("Could not find pathway column: ", pathway_col)
  }
  
  n_pathway_matches <- 0L
  n_compound_pathway_matches <- 0L
  
  for (j in seq_len(nrow(pathway_results))) {
    
    pathway_id_j <- as.character(pathway_results[[pathway_col]][j])
    pathway_id_j <- trimws(pathway_id_j)
    
    hits_j <- suppressWarnings(as.numeric(pathway_results[[hits_col]][j]))
    total_j <- suppressWarnings(as.numeric(pathway_results[[total_col]][j]))
    
    if (is.na(hits_j) || is.na(total_j) || total_j == 0) {
      next
    }
    
    if (!pathway_id_j %in% names(organism_map)) {
      next
    }
    
    n_pathway_matches <- n_pathway_matches + 1L
    
    pathway_compounds <- organism_map[[pathway_id_j]]$compound_ids
    pathway_compounds <- pathway_compounds[
      !is.na(pathway_compounds) & pathway_compounds != ""
    ]
    
    compounds_in_dataset_and_pathway <- intersect(
      compound_ids,
      pathway_compounds
    )
    
    if (length(compounds_in_dataset_and_pathway) == 0) {
      next
    }
    
    contribution_j <- hits_j / total_j
    
    score_df$compound_score[
      score_df$compound_id %in% compounds_in_dataset_and_pathway
    ] <- score_df$compound_score[
      score_df$compound_id %in% compounds_in_dataset_and_pathway
    ] + contribution_j
    
    n_compound_pathway_matches <- n_compound_pathway_matches +
      length(compounds_in_dataset_and_pathway)
  }
  
  if (verbose) {
    message("Pathway rows in result table: ", nrow(pathway_results))
    message("Pathway rows matched to nested KEGG map: ", n_pathway_matches)
    message("Compound-pathway matches found: ", n_compound_pathway_matches)
    message("Compounds with non-zero score: ", sum(score_df$compound_score > 0))
  }
  
  score_df <- score_df[
    order(-score_df$compound_score, score_df$compound_id),
  ]
  
  rownames(score_df) <- NULL
  
  score_df
}





metaspace_examples_script$compound_scores <- I(
  lapply(
    seq_len(nrow(metaspace_examples_script)),
    function(i) {
      score_compounds_using_nested_kegg_maps(
        compound_ids = metaspace_examples_script$kegg_compounds[[i]],
        pathway_results = metaspace_examples_script$pathway_results[[i]],
        organism_code = metaspace_examples_script$kegg_organism[i],
        kegg_pathway_nested_maps = kegg_pathway_nested_maps,
        verbose = FALSE
      )
    }
  )
)


all_ratios <- unlist(
  lapply(metaspace_examples_script$pathway_results, function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$Hits / df$Total
  }),
  use.names = FALSE
)
all_ratios <- all_ratios[is.finite(all_ratios)]
hist(all_ratios, breaks = 25)




# Contaminate molecules again!
# Do pathway analysis, and validate whether or not the scoring system actually works



# Creating a column with only molecules that exist in a pathway
clean_ids <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- toupper(x)
  x <- x[nzchar(x)]
  unique(x)
}

metaspace_examples_script$comp_in_pathways <- lapply(seq_len(nrow(metaspace_examples_script)), function(i) {
  org <- metaspace_examples_script$kegg_organism[i]
  pr  <- metaspace_examples_script$pathway_results[[i]]
  
  # choose one:
  # dataset_ids <- clean_ids(metaspace_examples_script$kegg_compounds[[i]])
  dataset_ids <- clean_ids(metaspace_examples_script$kegg_compounds_on_sample[[i]])
  
  if (is.null(pr) || nrow(pr) == 0 || length(dataset_ids) == 0) return(character(0))
  org_map <- kegg_pathway_nested_maps[[org]]
  if (is.null(org_map)) return(character(0))
  
  pids <- unique(as.character(pr$pathway_name))
  
  pathway_ids <- unlist(lapply(pids, function(pid) {
    node <- org_map[[pid]]
    if (is.null(node) || is.null(node$compound_ids)) return(character(0))
    node$compound_ids
  }), use.names = FALSE)
  
  pathway_ids <- clean_ids(pathway_ids)
  
  intersect(dataset_ids, pathway_ids)
})

summary(lengths(metaspace_examples_script$comp_in_pathways))

metaspace_examples_script$fraction <- vapply(
  metaspace_examples_script$pathway_results,
  function(df) {
    if (is.null(df) || nrow(df) == 0) return(NA_real_)
    sum(df$Hits, na.rm = TRUE) / sum(df$Total, na.rm = TRUE)
  },
  numeric(1)
)

metaspace_examples_script$fraction_list <- I(lapply(
  metaspace_examples_script$pathway_results,
  function(df) {
    if (is.null(df) || nrow(df) == 0) return(numeric(0))
    as.numeric(df$Hits) / as.numeric(df$Total)
  }
))

metaspace_examples_script$fraction_median <- lapply(
  metaspace_examples_script$fraction_list,
  function(df) {
    as.numeric(median(df))
  }
)


hist(unlist(metaspace_examples_script$fraction_median), breaks = 25)
hist(metaspace_examples_script$fraction, breaks = 25)
hist(all_ratios, breaks = 100, main = "Fraction of molecules present in pathways", xlab = "Fraction (F)")

# There is going to be overlap between the molecules that are detected in each pathways so how exactly are we supposed to calculate a fraction for these guys????



sample_10 <- metaspace_examples_script[sample(nrow(metaspace_examples_script), 10), ]





# Contaminating

sample_10$contaminated_list <- lapply(sample_10$comp_in_pathways, function(comp_vec) {
  unlist(
    lapply(comp_vec, function(kegg_id) {
      mass <- getMolActualWeight(kegg_id)
      Sys.sleep(0.1)
      additional <- getAdditionalMol(kegg_id, mass)
      
      # collect original + new ids
      unique(c(kegg_id, additional$kegg_id))
    }),
    use.names = FALSE
  ) |> unique()
})




# Metaboanalyst

library(MetaboAnalystR)

run_pathway_analysis <- function(compound_vec, organism = "hsa") {
  
  # Initialize analysis object
  mSet <- InitDataObjects("conc", "pathora", FALSE, default.dpi = 300)
  
  # Provide compound list
  mSet <- Setup.MapData(mSet, compound_vec)
  
  # Map compounds to KEGG
  mSet <- CrossReferencing(mSet, "kegg")
  mSet <- CreateMappingResultTable(mSet)
  
  # Set pathway library
  mSet <- SetMetabolomeFilter(mSet, FALSE)
  mSet<-SetKEGG.PathLib(mSet, organism, "current")
  
  # Perform enrichment + topology
  mSet <- CalculateOraScore(mSet, "rbc", "hyperg")
  mSet <- CalculateTopoScore(mSet, "degree")
  
  # Return results
  mSet$analSet$path.mat
}


sample_10$contaminated_pathways <- mapply(
  FUN = run_pathway_analysis,
  compound_vec = sample_10$contaminated_list,
  organism = sample_10$kegg_organism,
  SIMPLIFY = FALSE
)


