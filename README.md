# FeatureMSEA <img src="man/figures/logo.png" align="right" alt="" width="120" />

**FeatureMSEA** implements Feature-Based Metabolite Set Enrichment Analysis — a method that performs functional interpretation directly from ranked mass spectrometry features, without requiring complete metabolite identification or an arbitrary significance threshold to preselect features.

## Installation

```r
# Install TidyMass dependencies first
remotes::install_gitlab("tidymass/masstools")
remotes::install_gitlab("tidymass/massdataset")
remotes::install_gitlab("tidymass/metid")

# Install featuremsea
remotes::install_github("tidymass/featuremsea")

library(masstools)
library(massdataset)
library(metid)
library(featuremsea)
```


## Step 1 — Data Preparation

The primary input for featureMSEA is a feature table with variable_id (unique identifier), m/z, rt (retention time in seconds), condition (phenotype-associated ranking statistic such as absolute signal-to-noise ratio (SNR), absolute log2 (fold change) (log2FC), and absolute correlation coefficient), polarity mode ("positive" or "negative"), and mean_intensity (average feature intensity across samples). The analysis additionally requires an MS1 metabolite database and a metabolite set database, which can be downloaded from TidyMass website.

- **Feature table**: [Demo data](https://drive.google.com/file/d/1UwjizHDok-k9yYI407IZ_IF-YHkZjVjI/view?usp=drive_link)
- **MS1 database** and **Metabolite set databases**: [Download here](https://www.tidymass.org/databases/)

---

## Step 2 — MS1-based Annotation


```r
annotation_table_final <-
  annotate_feature_table(
    feature_table = feature_table,
    column = "rp",
    metabolite_database = kegg_compound_ms1, # use downloaded ms1 database
    database_type = "KEGG", # select database type :KEGG or HMDB
    ms1_match_ppm = 15, 
    mfc_rt_tol = 10,
    isotope_number = 3
  )
```

---

## Step 3 — Remove redundant annotations

```r
annotation_table_final2 =
  featuremsea::remove_redundancy(
    annotation_table = annotation_table_final)

results <- process_annotation_table(
  annotation_table_final2 = annotation_table_final2,
  database_type = "KEGG" # KEGG or HMDB here
)

ranking_table <- results$ranking_table
annotation_table <- results$original_score_annotation

```

---

## Step 4 — FeatureMSEA Analysis


```r
# Supported databases: KEGG, SMPDB, IMETPD, Reactome, Wikipathway
results <- perform_fmsea_analysis(
  pathway_database, # use downloaded metabolite set database
  annotation_table,         
  ranking_table,
  threads = 6,
  min.compounds.num = 15,
  max.compounds.num = 300,
  id.col = "KEGG_ID",
  perm.num = 10000,
  seed = 123,
  fdr.thr = 0.05,
  max.iter.num = 3,
  verbose = TRUE
)


# Significant pathways are in:
results@significant_modules
```

---

## Step 5 — LLM-assisted interpretation (Optional)

### 5a. Matrix confidence assessment

Assesses how reliably metabolites in your sample matrix (urine, plasma, etc.) indicate each pathway's activity.

```r
# Using OpenAI
fmsea_result <- analyze_matrix_relevance(
  results,
  sample_source = "urine",          # matrix source of features
  api_key       = "sk-openai-xxx",  # OpenAI API key
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

Added columns to `significant_modules`:
- `matrix_confidence_score` — integer: 0 / 25 / 50 / 75 / 100
- `matrix_confidence_reason` — brief explanation
- `matrix_source` — the sample source used

### 5b. Topic relevance assessment

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

Added columns to `significant_modules`:
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

> FeatureMSEA: Metabolic Feature-based Metabolite Set Enrichment Analysis

