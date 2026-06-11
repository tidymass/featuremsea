# featureMSEA <img src="man/figures/logo.png" align="right" alt="" width="120" />

**featureMSEA** implements Feature-Based Metabolite Set Enrichment Analysis (fMSEA) — a method that runs metabolic pathway enrichment directly from mass spectrometry features, bypassing the need for complete metabolite identification. It is part of the [TidyMass](https://github.com/tidymass) ecosystem.

## Installation

```r
# Install TidyMass dependencies first
remotes::install_gitlab("tidymass/masstools")
remotes::install_gitlab("tidymass/massdataset")
remotes::install_gitlab("tidymass/metid")

# Install featuremsea
remotes::install_github("tidymass/featuremsea")
```


## Step 1 — Data Preparation

Three datasets are required:

- **Feature table**: [Demo data](https://drive.google.com/file/d/1UwjizHDok-k9yYI407IZ_IF-YHkZjVjI/view?usp=drive_link)
- **MS1 database** and **Metabolite set databases**: [Download here](https://www.tidymass.org/databases/)

---

## Step 2 — Annotate Features

`annotate_feature_table()` combines MS2 annotations with isotope/adduct grouping (MFC) to generate a scored annotation table.

```r
annotation_table <- annotate_feature_table(
  feature_table        = feature_table,
  annotation_table_ms2 = annotation_table_ms2,
  column               = "rp",          # "rp" (reverse-phase) or "hilic"
  database_type        = "KEGG",        # "KEGG" or "HMDB"
  metabolite_database  = kegg_database, # metid database object
  ms1_match_ppm        = 15,
  mfc_rt_tol           = 5,
  isotope_number       = 3
)
```

---

## Step 3 — Process Annotation Table

`process_annotation_table()` reshapes the scored annotation into a feature × metabolite matrix and extracts the ranking weights used by the enrichment algorithm.

```r
processed <- process_annotation_table(
  annotation_table_final2 = annotation_table,
  database_type            = "KEGG"   # must match Step 2
)

score_matrix   <- processed$original_score_annotation  # features × metabolites
ranking_table  <- processed$ranking_table              # features ranked by condition
```

---

## Step 4 — Run fMSEA Analysis

`perform_fmsea_analysis()` runs the iterative enrichment algorithm against a pathway database.

```r
# Supported databases: KEGG, SMPDB, IMETPD, Reactome, Wikipathway
fmsea_result <- perform_fmsea_analysis(
  pathway_database    = kegg_pathway_db,  # S4 pathway database object
  annotation_table    = score_matrix,
  ranking_table       = ranking_table,
  threads             = 4,
  min.compounds.num   = 15,
  max.compounds.num   = 300,
  id.col              = "KEGG_ID",
  perm.num            = 1000,
  seed                = 123,
  fdr.thr             = 0.05,
  max.iter.num        = 3,
  verbose             = TRUE
)

# Significant pathways are in:
fmsea_result@significant_modules
```

---

## Step 5 — LLM-Assisted Evaluation (Optional)

### 5a. Matrix Confidence Scoring

Assesses how reliably metabolites in your sample matrix (urine, plasma, etc.) indicate each pathway's activity.

```r
# Using OpenAI
fmsea_result <- analyze_matrix_relevance(
  results       = fmsea_result,
  sample_source = "urine",          # "urine", "plasma", "serum", "blood", "feces"
  api_key       = "sk-openai-xxx",
  provider      = "openai"
)

# Using SiliconFlow / Qwen (recommended for users in China)
fmsea_result <- analyze_matrix_relevance(
  results       = fmsea_result,
  sample_source = "plasma",
  api_key       = "sk-siliconflow-xxx",
  provider      = "siliconflow"
)
```

Adds columns to `significant_modules`:
- `matrix_confidence_score` — integer: 0 / 25 / 50 / 75 / 100
- `matrix_confidence_reason` — brief explanation
- `matrix_source` — the sample source used

### 5b. Topic Relevance via PubMed + LLM

Links significant pathways to a research topic using PubMed literature search and embedding-based re-ranking.

```r
fmsea_result <- analyze_topic_relevance(
  results          = fmsea_result,
  research_topic   = "type 2 diabetes",
  api_key          = "sk-siliconflow-xxx",
  provider         = "siliconflow",
  pubmed_api_key   = NULL,          # optional, increases PubMed rate limit
  similarity_cutoff = 0.6           # cosine similarity threshold for fuzzy matches
)
```

Adds columns to `significant_modules`:
- `literature_pmids_exact` — PMIDs from exact PubMed search
- `literature_pmids_fuzzy` — PMIDs from fuzzy search filtered by embedding similarity
- `topic_confidence_score` — LLM score (0/25/50/75/100) when no literature found
- `research_topic` — the topic used

> **API provider choice:** Both functions require explicit `provider` selection — either `"openai"` or `"siliconflow"`. Default models: `gpt-5.5` / `Qwen/Qwen3-32B` (chat); `text-embedding-3-small` / `Qwen/Qwen3-Embedding-8B` (embedding). Custom models can be specified via the `model` argument.

---

## Step 6 — Visualization

```r
# Plot enrichment results for a specific pathway
plot_fmsea_plot(
  fmsea_obj  = fmsea_result,
  pathway_id = "hsa00010",          # pathway ID from significant_modules
  title      = "Glycolysis / Gluconeogenesis"
)
```

The plot shows the enrichment score curve, ranked feature hits, annotation weights, and pair weights across all features.

---

## Citation

If you use featureMSEA in your research, please cite:

> Shen X, Liu Y, et al. featureMSEA: Feature-Based Metabolite Set Enrichment Analysis. *TidyMass Project* (2026). <https://github.com/tidymass/featuremsea>

## Issues

Please report bugs at <https://github.com/tidymass/featuremsea/issues>.
