library(MetaboAnalystR)
tmp.vec <- c("Acetoacetic acid", "Beta-Alanine", "Creatine", "Dimethylglycine", "Fumaric acid", "Glycine", "Homocysteine", "L-Cysteine", "L-Isoleucine", "L-Phenylalanine", "L-Serine", "L-Threonine", "L-Tyrosine", "L-Valine", "Phenylpyruvic acid", "Propionic acid", "Pyruvic acid", "Sarcosine")
mSet<-InitDataObjects("conc", "pathora", FALSE, default.dpi = 300)
mSet<-Setup.MapData(mSet, tmp.vec);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet);

unmatched_indices <- which(mSet$name.map$match.state == 0)
unmatched_queries <- mSet$name.map$query.vec[unmatched_indices]


# Loop through each unmatched compound
for (i in seq_along(unmatched_indices)) {
  query <- unmatched_queries[i]
  
  # Get candidates for THIS specific compound
  mSet <- PerformDetailMatch(mSet, query)
  mSet <- GetCandidateList(mSet)
  
  # Now extract candidates from the matrix (column 1)
  candidates <- mSet$name.map$hits.candidate.list[, 1]
  candidates <- candidates[candidates != ""]  # Remove empty rows
  
  print(paste("Candidates for", query, ":"))
  print(candidates)
  choice <- as.integer(readline(prompt = "Enter the candidate number you want to use: "))
  
  if (!is.na(choice) && choice > 0 && choice <= length(candidates)) {
    chosen_candidate <- candidates[choice]
    mSet <- SetCandidate(mSet, query, chosen_candidate)
    print(paste("Selected:", chosen_candidate))
  } else {
    cat("Invalid choice, skipping...\n")
  }
}

# # Print detailed mapping info
# for (i in seq_along(unmatched_indices)) {
#   idx <- unmatched_indices[i]
#   query <- unmatched_queries[i]
#   candidates <- mSet$name.map$hits.candidate.list[[i]]
#   
#   print(paste("Index:", idx, "| Query:", query, "| Candidates:", paste(candidates, collapse=", ")))
# }



mSet$name.map$hits.candidate.list <- NULL

mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateOraScore(mSet, "rbc", "hyperg")

mSet$analSet$ora.mat

library(KEGGREST)

# Specify your pathway ID
kegg_id <- "hsa00260"  # Example: Glycine, serine and threonine metabolism

# Query KEGG for pathway details
pw_info <- keggGet(kegg_id)[[1]]

# Extract the COMPOUND list (named vector: KEGG ID = compound name)
kegg_comp_ids <- names(pw_info$COMPOUND)
kegg_comp_names <- unname(pw_info$COMPOUND)

# Print all KEGG compound IDs in the pathway
print(kegg_comp_ids)

df <- data.frame(
  KEGG_ID = names(pw_info$COMPOUND),
  Name = as.vector(pw_info$COMPOUND),
  row.names = NULL
)
print(df)

kegg_ids <- rownames(mSet$analSet$ora.mat)  # List of pathway IDs

print(kegg_ids)

# Convert to data.frame if needed (character matrix may need this)
map_table_df <- as.data.frame(mSet$dataSet$map.table, stringsAsFactors = FALSE)

# Your input compound KEGG IDs (removing blanks)
input_keggs <- unique(map_table_df$KEGG[map_table_df$KEGG != ""])




chunk10 <- function(x) split(x, ceiling(seq_along(x)/10))

get_kegg_formulas <- function(compound_ids, sleep = 0.0) {
  compound_ids <- unique(compound_ids)
  compound_ids <- compound_ids[!is.na(compound_ids) & compound_ids != ""]
  if (length(compound_ids) == 0) return(setNames(character(0), character(0)))
  
  batches <- chunk10(compound_ids)
  
  recs <- unlist(lapply(batches, function(b) {
    if (sleep > 0) Sys.sleep(sleep)
    # keggGet returns a list of records (one per ID)
    out <- KEGGREST::keggGet(b)
    out
  }), recursive = FALSE)
  
  # Build a mapping: id -> FORMULA
  formulas <- sapply(recs, function(rec) {
    if (!is.null(rec$FORMULA)) rec$FORMULA else NA_character_
  })
  
  ids <- sapply(recs, function(rec) rec$ENTRY)
  # ENTRY looks like "C00037                      Compound"
  ids <- sub("\\s+.*$", "", ids)
  
  setNames(as.character(formulas), ids)
}







all_pathway_lists <- lapply(kegg_ids, function(pid) {
  pw <- keggGet(pid)[[1]]
  if (is.null(pw$COMPOUND)) return(data.frame())
  
  kegg_ids_this <- names(pw$COMPOUND)
  comp_names <- as.vector(pw$COMPOUND)
  in_input <- kegg_ids_this %in% input_keggs
  
  # Fetch formulas for all compounds in this pathway (batched by 10)
  formula_map <- get_kegg_formulas(kegg_ids_this, sleep = 0.1)
  
  mol_formula <- unname(formula_map[kegg_ids_this])  # aligned to same order
  # If any are missing, they become NA (that's fine)
  
  data.frame(
    kegg_id = kegg_ids_this,
    compound_name = comp_names,
    in_input = in_input,
    mol_formula = mol_formula,
    stringsAsFactors = FALSE
  )
})

names(all_pathway_lists) <- kegg_ids




# # 1) collect all compounds across pathways
# all_comp_ids <- unique(unlist(lapply(kegg_ids, function(pid) {
#   pw <- keggGet(pid)[[1]]
#   if (is.null(pw$COMPOUND)) return(character(0))
#   names(pw$COMPOUND)
# })))
# 
# # 2) fetch all formulas once (batched)
# formula_map_all <- get_kegg_formulas(all_comp_ids, sleep = 0.1)
# 
# # 3) build pathway tables using the cached map
# all_pathway_lists <- lapply(kegg_ids, function(pid) {
#   pw <- keggGet(pid)[[1]]
#   if (is.null(pw$COMPOUND)) return(data.frame())
# 
#   kegg_ids_this <- names(pw$COMPOUND)
#   data.frame(
#     kegg_id = kegg_ids_this,
#     compound_name = as.vector(pw$COMPOUND),
#     in_input = kegg_ids_this %in% input_keggs,
#     mol_formula = unname(formula_map_all[kegg_ids_this]),
#     stringsAsFactors = FALSE
#   )
# })
# names(all_pathway_lists) <- kegg_ids

save(all_pathway_lists, file="data.Rda")

# library(OrgMassSpecR)
# OrgMassSpecR::monoisotopic("C2H5NO2")
# 
# all_pathway_lists <- lapply(all_pathway_lists, function(df) {
#   if (nrow(df) == 0) return(df)
#   df$mono_mass <- vapply(df$mol_formula, function(f) {
#     if (is.na(f) || f == "") return(NA_real_)
#     as.numeric(OrgMassSpecR::monoisotopic(f))
#   }, numeric(1))
#   df
# })


library(Rdisop)


calc_exact_mass <- function(formula) {
  if (is.na(formula) || formula == "") return(NA_real_)
  formula <- gsub("\\s+", "", formula)  # remove spaces if any
  # getMolecule() returns a list with $exactmass
  out <- tryCatch(Rdisop::getMolecule(formula)$exactmass, error = function(e) NA_real_)
  as.numeric(out)
}

all_pathway_lists <- lapply(all_pathway_lists, function(df) {
  if (nrow(df) == 0) return(df)
  df$mono_mass <- vapply(df$mol_formula, calc_exact_mass, numeric(1))
  df
})
