.kegg_cache <- new.env(parent = emptyenv())
cache_key <- function(mass, tol) paste0(format(mass, digits = 12), "_", tol)

library(KEGGREST)

FindMolEqualWeight <- function(query_weight, tolerance = 1000){
  key <- cache_key(query_weight, tolerance)
  if (exists(key, envir = .kegg_cache, inherits = FALSE)) return(get(key, envir = .kegg_cache))
  
  tolerance_da = (query_weight/1e6)*tolerance # Calculates tolerance in daltons
  print(tolerance_da)
  lower_limit = query_weight - tolerance_da
  upper_limit = query_weight + tolerance_da
  similar <- keggFind("compound", floor(lower_limit):ceiling(upper_limit), "exact_mass")
  
  # Defines dataframe for result
  df <- data.frame( 
    kegg_id = names(similar),
    mol_mass = as.numeric(as.vector(similar)),
    stringsAsFactors = FALSE
  )
  
  df$delta_da <- abs(df$mol_mass - query_weight)
  
  df <- df[df$delta_da <= tolerance_da, , drop = FALSE]

  df <- df[order(df$delta_da), , drop = FALSE]
  # if (nrow(df) > max_return) df <- df[seq_len(max_return), , drop = FALSE]

  df$query_weight <- rep(query_weight, nrow(df))
  df$tolerance_ppm <- rep(tolerance, nrow(df))
  df$tolerance_da <- rep(tolerance_da, nrow(df))
  rownames(df) <- NULL
  
  assign(key, df, envir = .kegg_cache)
  return(df)
}





# Removes the compound that is already in the pathway from the list
getAdditionalMol <- function(query_kegg, query_weight, tolerance = 1000, max_return = 5){
  df <- FindMolEqualWeight(query_weight, tolerance = tolerance)
  df <- df[df$kegg_id != query_kegg, , drop = FALSE]
  
  return(df)
}






# Creating a dataframe with all combinations of tolerances and percentage of pathway compounds
fracs <- list(0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1)
tolerances <- list(1, 3, 5, 10, 30, 60, 100)









sampleCompounds <- function(pathway_id, frac, pathway_db = all_pathway_lists) {
  df <- pathway_db[[pathway_id]][, c("kegg_id", "mono_mass")]
  
  # remove missing/blank ids
  df <- df[!is.na(df$kegg_id) & df$kegg_id != "", , drop = FALSE]
  
  # ensure unique compounds (keep first row per kegg_id)
  df <- df[!duplicated(df$kegg_id), , drop = FALSE]
  
  n_total <- nrow(df)
  if (n_total == 0) return(df)
  
  size <- max(1L, floor(n_total * frac))
  df <- df[sample.int(n_total, size), , drop = FALSE]
  
  data.frame(
    kegg_id = df$kegg_id,
    mol_mass = as.numeric(df$mono_mass),
    stringsAsFactors = FALSE
  )
}

sampleCompounds("hsa00260", 0.5)




testing_main_db <- list()

fracs <- c(1, 0.75, 0.5, 0.2, 0.1, 0.05, 0.01)
tolerances <- c(1, 3, 5, 10, 30, 60, 100)
for (pathway in names(all_pathway_lists)) {
  testing_main_db[[pathway]] <- list()
  
  for (frac in fracs) {
    frac_key <- as.character(frac)
    testing_main_db[[pathway]][[frac_key]] <- list()
    
    tempSampleCompounds <- sampleCompounds(pathway, frac = frac)
    testing_main_db[[pathway]][[frac_key]][["sample"]] <- tempSampleCompounds
    
    for (tol in tolerances) {
      tol_key <- as.character(tol)
      
      out <- tempSampleCompounds
      
      for (i in seq_len(nrow(tempSampleCompounds))) {
        pathway_compound <- tempSampleCompounds$kegg_id[i]
        comp_mass <- tempSampleCompounds$mol_mass[i]
        if (!is.finite(comp_mass)) next
        
        add_mol_max <- getAdditionalMol(pathway_compound, comp_mass, tolerance = 100)
        tol_da <- (comp_mass / 1e6) * tol
        
        add_mol <- add_mol_max[add_mol_max$delta_da <= tol_da, , drop = FALSE]
        add_mol$delta_da <- NULL
        add_mol$query_weight <- NULL
        add_mol$tolerance_ppm <- NULL
        add_mol$tolerance_da <- NULL
        
        out <- rbind(out, add_mol)
      }
      
      out <- out[!duplicated(out$kegg_id), , drop = FALSE]
      testing_main_db[[pathway]][[frac_key]][[tol_key]] <- out
    }
  }
}

saveRDS(testing_main_db, file="Main_Testing.Rda")


# Get a graph of number of molecules against T for a set of 50 random query molecules


# Adding an in_pathway column
for (pathway in names(testing_main_db)) {
  for (frac in names(testing_main_db[[pathway]])){
    for (tol in names(testing_main_db[[pathway]][[as.character(frac)]])){
      
      df <- testing_main_db[[pathway]][[frac]][[tol]]
      
      df$in_pathway <- df$kegg_id %in% testing_main_db[[pathway]][[frac]]$sample$kegg_id
      
      testing_main_db[[pathway]][[frac]][[tol]] <- df
    }
  }
}

library(MetaboAnalystR)

sampled_pathways = c("hsa00250", "hsa00480", "hsa00330", "hsa00430", "hsa00470", "hsa00650", "hsa00130", "hsa00410", "hsa00010", "hsa00860", "hsa00600", "hsa00120") # "hsa00785" "hsa00350" "hsa00730" "hsa00240" "hsa00220" "hsa00260", "hsa00670", "hsa00400", "hsa00290", "hsa00360" "hsa00270" "hsa00770" "hsa00630" "hsa00280", "hsa00020", "hsa00640", "hsa00620"
sampled_tol = as.character(c(1, 3, 5, 10, 30, 60, 100))
sampled_fracs = as.character(c(1, 0.75, 0.5, 0.2, 0.1, 0.05, 0.01)) # 1, 0.75, 0.5, 0.2, 0.1, 0.05, 0.01

# all_pathway_results <- vector("list", length(sampled_pathways))
# names(all_pathway_results) <- sampled_pathways

for (pathway in sampled_pathways) {
  for (frac in sampled_fracs) {
    for (tol in sampled_tol) {
      
      ids <- testing_main_db[[pathway]][[frac]][[tol]]$kegg_id
      ids <- unique(ids[!is.na(ids)])
      
      # Ensure nested list exists
      if (is.null(all_pathway_results[[pathway]]) || !is.list(all_pathway_results[[pathway]]))
        all_pathway_results[[pathway]] <- list()
      if (is.null(all_pathway_results[[pathway]][[frac]]) || !is.list(all_pathway_results[[pathway]][[frac]]))
        all_pathway_results[[pathway]][[frac]] <- list()
      if (is.null(all_pathway_results[[pathway]][[frac]][[tol]]) || !is.list(all_pathway_results[[pathway]][[frac]][[tol]]))
        all_pathway_results[[pathway]][[frac]][[tol]] <- list()
      # if (is.list(all_pathway_results[[pathway]]$`5`))
      #   all_pathway_results[[pathway]]$`5` <- NULL
      # I had added the other line to fix an earlier error in my code
      
      
      if (length(ids) < 3) {
        all_pathway_results[[pathway]][[frac]][[tol]]$pathway_anal <- NULL
        all_pathway_results[[pathway]][[frac]][[tol]]$status <- "skipped"
        all_pathway_results[[pathway]][[frac]][[tol]]$reason <- paste0("Only ", length(ids), " unique non-NA compounds (< 3)")
        next
      }
      
      res <- tryCatch({
        oSet <- InitDataObjects("conc", "pathora", FALSE, default.dpi = 300)
        oSet <- Setup.MapData(oSet, as.vector(ids))
        oSet <- CrossReferencing(oSet, "kegg")
        oSet <- CreateMappingResultTable(oSet)
        oSet <- SetKEGG.PathLib(oSet, "hsa", "current")
        oSet <- SetMetabolomeFilter(oSet, FALSE)
        oSet <- CalculateOraScore(oSet, "rbc", "hyperg")
        oSet$analSet$ora.mat
      }, error = function(e) e)
      
      if (inherits(res, "error")) {
        all_pathway_results[[pathway]][[frac]][[tol]]$pathway_anal <- NULL
        all_pathway_results[[pathway]][[frac]][[tol]]$status <- "error"
        all_pathway_results[[pathway]][[frac]][[tol]]$reason <- conditionMessage(res)
      } else {
        all_pathway_results[[pathway]][[frac]][[tol]]$pathway_anal <- res
        all_pathway_results[[pathway]][[frac]][[tol]]$status <- "ok"
        all_pathway_results[[pathway]][[frac]][[tol]]$n_compounds <- length(ids)
      }
    }
  }
}

save(all_pathway_results, file="Pathway analysis with parameters.Rda")

# Select 50 random compounds
# Run FindMolEqualWeight() for each with very high value of tolerance
# Filter for smaller values of T (saving a bunch of KEGGfind requests)

all_compounds <- data.frame(kegg_id=character(), mol_mass=numeric(), stringsAsFactors = FALSE)

for (pathway in names(all_pathway_lists)) {
  temp <- data.frame(
    kegg_id  = all_pathway_lists[[pathway]][["kegg_id"]],
    mol_mass = all_pathway_lists[[pathway]][["mono_mass"]],
    stringsAsFactors = FALSE
  )
  all_compounds <- rbind(all_compounds, temp)
}


df_ok <- all_compounds[!is.na(all_compounds$mol_mass), c("kegg_id", "mol_mass")] # Creating a dataframe with no compounds that have molecular mass NA
unique_df_ok <- df_ok[!duplicated(df_ok$kegg_id), ] # Making sure that all compounds are unique

# Sampling 50 compounds from the list
set.seed(1)
n <- min(50, nrow(unique_df_ok))
sampled_for_z_t <- unique_df_ok[sample(seq_len(nrow(unique_df_ok)), n), ]


Z_T_relationship <- list()

for (i in seq_len(nrow(sampled_for_z_t))) {
  
  pathway_compound <- sampled_for_z_t$kegg_id[i]
  comp_mass <- sampled_for_z_t$mol_mass[i]
  
  Z_T_relationship[[pathway_compound]] <- list()
  
  at_max_tol <- FindMolEqualWeight(comp_mass, tolerance = 2000)
  t <- seq(2000, 0, by = -1)
  
  Z <- numeric(length(t))
  
  for (j in seq_along(t)) {
    tolerance <- t[j]
    tolerance_da <- (comp_mass / 1e6) * tolerance
    Z[j] <- nrow(at_max_tol[at_max_tol$delta_da <= tolerance_da, , drop = FALSE])
  }
  
  Z_T_relationship[[pathway_compound]] <- list(T = t, Z = Z)
}






kegg_ids <- names(Z_T_relationship)

all_df_zt <- data.frame(
  kegg_id = kegg_ids,
  stringsAsFactors = FALSE
)

all_df_zt$tbl <- lapply(kegg_ids, function(k) {
  data.frame(T = Z_T_relationship[[k]]$T,
             Z = Z_T_relationship[[k]]$Z)
})

# all_df_zt
# # access one:
# all_df_zt$tbl[[ all_df_zt$kegg_id == "C00380" ]]





#lm(formula = Z ~ T) # for linear regression (do for 1<T<400)





# Extracting all pathways for which FDR<=0.05



all_pathway_results <- lapply(all_pathway_results, function(x) {
  lapply(x, function(y) {
    lapply(y, function(z) {
      
      mat <- z$pathway_anal
      
      if (!is.null(mat)) {
        sig <- mat[mat[, "FDR"] <= 0.05, , drop = FALSE]
        
        if (nrow(sig) > 0) {
          sig <- sig[order(sig[, "FDR"]), , drop = FALSE]
        }
        
        z$sig_pathways <- sig
      } else {
        z$sig_pathways <- NULL
      }
      
      z$n_sig <- if (!is.null(mat)) sum(mat[, "FDR"] <= 0.05) else 0
      z
    })
  })
})



# Get ordered names (already sorted earlier)
frac_names <- names(all_pathway_results[[1]])
tol_names  <- names(all_pathway_results[[1]][[1]])

# Initialize matrices
mat_sig_count   <- matrix(0, nrow = length(tol_names), ncol = length(frac_names),
                          dimnames = list(tol_names, frac_names))

mat_top1        <- mat_sig_count
mat_top5        <- mat_sig_count
mat_total_sig   <- mat_sig_count
mat_top1_exfdr  <- mat_sig_count

# Loop over everything
for (pathway in names(all_pathway_results)) {
  for (f in frac_names) {
    for (t in tol_names) {
      
      entry   <- all_pathway_results[[pathway]][[f]][[t]]
      sig     <- entry$sig_pathways
      pathways <- entry$pathway_anal
      
      sig_names <- rownames(sig)
      
      ## (extra) Top 1 even when FDR isn't significant
      if (!is.null(pathways) && NROW(pathways) > 0) {
        pathway_names <- rownames(pathways)
        if (!is.null(pathway_names) && length(pathway_names) >= 1 && pathway_names[1] == pathway) {
          mat_top1_exfdr[t, f] <- mat_top1_exfdr[t, f] + 1
        }
      }
      
      
      if (!is.null(sig) && NROW(sig) > 0 && length(sig_names) > 0) {
        
        
        
        ## 1. Count if THIS pathway appears in its own sig list
        if (pathway %in% sig_names) {
          mat_sig_count[t, f] <- mat_sig_count[t, f] + 1
        }
        
        ## 2. Top 1
        if (length(sig_names) >= 1 && sig_names[1] == pathway) {
          mat_top1[t, f] <- mat_top1[t, f] + 1
        }
        
        ## 3. Top 5
        top_n <- min(5, length(sig_names))
        if (pathway %in% sig_names[1:top_n]) {
          mat_top5[t, f] <- mat_top5[t, f] + 1
        }
        
        ## 4. Total number of significant entries
        mat_total_sig[t, f] <- mat_total_sig[t, f] + nrow(sig)
      }
    }
  }
}


## Finding the pathway that seems to be missing at frac = 1 and tol = 100

all_pathways <- names(all_pathway_results)

present <- sapply(all_pathways, function(pw) {
  sig <- all_pathway_results[[pw]][["frac_0.5"]][["tol_1"]]$sig_pathways
  
  if (is.null(sig) || nrow(sig) == 0) return(FALSE)
  
  pw %in% rownames(sig)
})

missing_pathways <- all_pathways[!present]

missing_pathways







# Function to recursively extract all data frames with columns kegg_id and mol_mass
extract_kegg_dfs <- function(x) {
  out <- list()
  
  if (is.data.frame(x) && all(c("kegg_id", "mol_mass") %in% names(x))) {
    out <- list(x[, c("kegg_id", "mol_mass")])
  } else if (is.list(x)) {
    out <- unlist(lapply(x, extract_kegg_dfs), recursive = FALSE)
  }
  
  out
}

# Extract all matching data frames
all_dfs <- extract_kegg_dfs(testing_main_db)

# Combine into one table
all_kegg <- bind_rows(all_dfs)

# Remove exact duplicate rows (same kegg_id + same mol_mass)
unique_kegg <- distinct(all_kegg, kegg_id, mol_mass)

# View result
unique_kegg



library(dplyr)

one_row_per_kegg <- all_kegg %>%
  group_by(kegg_id) %>%
  summarise(
    mol_mass = {
      vals <- unique(na.omit(mol_mass))
      if (length(vals) == 0) NA_real_ else vals[1]
    },
    n_unique_masses = n_distinct(na.omit(mol_mass)),
    .groups = "drop"
  )

one_row_per_kegg

unique_kegg <- one_row_per_kegg





# Plot a histogram


hist(
  na.omit(one_row_per_kegg$mol_mass),
  breaks = 600,
  main = "Histogram of Molecular Masses",
  # xlim = range(258:261),
  xlab = "Molecular Mass",
  ylab = "Number of Compounds",
  col = "skyblue",
  border = "white"
)





library(ggplot2)

ggplot(one_row_per_kegg, aes(x = mol_mass)) +
  geom_histogram(bins = 500, fill = "skyblue", color = "black", na.rm = TRUE) +
  labs(
    title = "Histogram of Molecular Masses",
    x = "Molecular Mass",
    y = "Number of Compounds"
  ) +
  theme_minimal()

unique_kegg_200_300 <- unique_kegg[unique_kegg$mol_mass>520 & unique_kegg$mol_mass<550, ]

choose(30, 2)




num_mol_chance <- function(num_candidates, pathway_size, db_size = 19572) {
  return((num_candidates * pathway_size) / db_size)
}

frac_names <- names(all_pathway_results[[1]])
tol_names  <- names(all_pathway_results[[1]][[1]])

for (pathway in names(all_pathway_results)) {
  for (f in frac_names) {
    for (t in tol_names) {
      
      pathway_size <- NROW(testing_main_db[[pathway]][["frac_1"]][["tol_sample"]])
      
      all_pathway_results[[pathway]][[f]][[t]]$expected_values <- data.frame(
        num_candidates = 1:4000,
        expected_mols = num_mol_chance(1:4000, pathway_size),
        frac_mols = (1:4000)/19572
      )
      
    }
  }
}

# To access:
# all_pathway_results[[pathway]][[f]][[t]]$expected_values[450, "expected_mols"]


# Plotting
# Binding all tables together:

all_ev <- do.call(rbind, lapply(names(all_pathway_results), function(pathway) {
  do.call(rbind, lapply(names(all_pathway_results[[pathway]]), function(f) {
    do.call(rbind, lapply(names(all_pathway_results[[pathway]][[f]]), function(t) {
      
      df <- all_pathway_results[[pathway]][[f]][[t]]$expected_values
      
      df$pathway <- pathway
      df$f <- f
      df$t <- t
      
      return(df)
      
    }))
  }))
}))


library(dplyr)

summary_df <- all_ev %>%
  group_by(num_candidates) %>%
  summarise(
    median = median(expected_mols, na.rm = TRUE),
    mean   = mean(expected_mols, na.rm = TRUE),
    sd     = sd(expected_mols, na.rm = TRUE),
    q10    = quantile(expected_mols, 0.10, na.rm = TRUE),
    q90    = quantile(expected_mols, 0.90, na.rm = TRUE)
  )



library(ggplot2)

ggplot(summary_df, aes(x = num_candidates)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
              fill = "grey", alpha = 0.4) +
  geom_line(aes(y = median), color = "blue", size = 1) +
  labs(
    title = "Median Expected Molecules with Variability",
    x = "Number of Candidates",
    y = "Expected Molecules"
  ) +
  theme_minimal()




# i = 30
# db = 20000
# 
# 30/20000 * ((30/20000)*(29/19999) + (19970/20000)*(30/19999))




# Look at metspace, take 20 datasets
num_compounds <- c(629, 138, 2593, 645, 3639, 141, 1378, 1278, 11, 609)

metaspace_examples <- data.frame(
  num_compounds = num_compounds,
  expected_mols = as.vector(all_pathway_results[["hsa00770"]][["frac_1"]][["tol_100"]]$expected_values[num_compounds, "expected_mols"]),
  expected_fracs = as.vector(all_pathway_results[["hsa00770"]][["frac_1"]][["tol_100"]]$expected_values[num_compounds, "frac_mols"]),
  dataset_id = c("2026-05-29_21h54m50s", 
                 "2026-05-13_14h12m22s", 
                 "2026-04-28_09h20m56s", 
                 "2026-04-24_09h36m39s", 
                 "2025-12-15_21h10m33s", 
                 "2026-04-10_01h17m48s", 
                 "2026-04-08_09h12m38s", 
                 "2026-04-08_09h10m11s", 
                 "2026-04-07_00h56m38s",
                 "2026-02-25_14h16m03s")
  # link = c("")
)
metaspace_examples








