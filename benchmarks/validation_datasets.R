# ==============================================================================
#  GeneSelectR 2.0 — PRE-REGISTERED validation dataset registry
# ==============================================================================
#
#  Five GEO cohorts, chosen and fully specified BEFORE any of them was run.
#  Everything that could be tuned after seeing a result is fixed here:
#
#    * which samples enter the analysis (filters, exclusions)
#    * which contrast is tested, and which class is "positive"
#    * the confounders that are residualised out
#    * the grouping variable for paired / repeated-measures designs
#    * `disease_term` and `target_go_terms`, chosen from the DISEASE NAME
#
#  The last one is the point of this file. `disease_term` seeds the Open
#  Targets/STRING biology pillar and `target_go_terms` seeds the semantic
#  pillar. Picking them after seeing which genes a method returned is invisible
#  tuning of the exact pillar the paper claims is external evidence — the
#  biology score stops being independent of the result it is used to justify.
#  They are written down here, once, and the benchmark reads them.
#
#  COMMITMENT: run all five, report all five, including the losses. A dataset
#  that is dropped after it is run is a dataset that was dropped for its result.
#
#  Consumed by:
#    benchmarks/validation_prepare.R    -- builds the analysis-ready matrices
#    benchmarks/validation_benchmark.R  -- runs the benchmark
#
#  This file defines data and pure helper functions ONLY. It is sourced by both
#  scripts and must have no side effects.
# ==============================================================================


# Null-coalescing operator, also defined in the callers; harmless to repeat and
# it lets this file be sourced on its own for inspection.
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x


# ------------------------------------------------------------------------------
#  Characteristic extraction by KEY, not by column position
# ------------------------------------------------------------------------------
#
#  GEO stores sample characteristics as "key: value" strings spread across
#  Sample_characteristics_ch1_* columns, and SERIES DO NOT GUARANTEE THE SAME
#  KEY ORDER FOR EVERY SAMPLE. fetch_geo_data.py assigns those columns
#  positionally, which is correct for the cohorts where the order is stable and
#  silently wrong where it is not.
#
#  On GSE20194 the positional parse put response values ("RD", "pCR") into the
#  ER-status column. Reading the outcome from there would have benchmarked the
#  wrong variable while looking entirely normal.
#
#  This reads the raw characteristic columns and pulls a key by name, per
#  sample. Returns NA where the key is absent.
extract_characteristic <- function(metadata, key) {
  chr_cols <- grep("^Sample_characteristics", colnames(metadata), value = TRUE)
  if (length(chr_cols) == 0) {
    stop(sprintf("No Sample_characteristics columns to read '%s' from.", key))
  }
  pattern <- paste0("^\\s*", key, "\\s*:")
  vapply(seq_len(nrow(metadata)), function(i) {
    vals <- as.character(unlist(metadata[i, chr_cols], use.names = FALSE))
    hit  <- grep(pattern, vals, value = TRUE, ignore.case = TRUE)
    if (length(hit) == 0) NA_character_ else trimws(sub("^[^:]*:", "", hit[1]))
  }, character(1))
}


# ------------------------------------------------------------------------------
#  Registry constructor
# ------------------------------------------------------------------------------
#
#  Every field is named so that reading one entry tells you the whole analysis
#  plan for that cohort without opening another file.
#
#  @param accession      GEO series accession, and the data/<accession>/ folder.
#  @param label          Human-readable cohort name for figures and logs.
#  @param scale          What the shipped matrix IS: "counts", "log_normalized",
#                        "linear_normalized", or "linear_background_corrected".
#                        VERIFIED against the matrix by
#                        validation_prepare.R, which errors on a mismatch rather
#                        than normalising the wrong thing quietly.
#  @param derive         Optional function(metadata) -> metadata, adding columns
#                        that must be parsed out of free text (patient IDs in a
#                        sample title, for instance). Runs FIRST, before filters.
#  @param sample_filters List of list(column=, keep=): rows whose value is not in
#                        `keep` are dropped. Applied in order, before the outcome
#                        mapping, so the log shows where each sample went.
#  @param outcome        list(column=, negative=, positive=, negative_label=,
#                        positive_label=). Raw values not listed in either arm
#                        are dropped. The factor is built with negative FIRST,
#                        because every metric in the benchmark treats level 2 as
#                        the positive class.
#  @param group_column   Column holding the unit that must not be split across a
#                        CV fold (patient, subject). NULL for independent
#                        samples. When set, the benchmark uses grouped folds and
#                        grouped stability subsamples.
#  @param confounders    Categorical metadata columns residualised out inside
#                        each training split. Empty = no residualisation.
#  @param disease_term   PRE-REGISTERED. Open Targets disease query string.
#  @param target_go_terms PRE-REGISTERED. GO BP terms the semantic pillar scores
#                        against. Chosen from the disease name and the tissue,
#                        never from a result.
#  @param expected_n     Named integer c(negative=, positive=). Checked after
#                        filtering; a mismatch is reported loudly. This is the
#                        pre-registration's arithmetic, and if the data does not
#                        reproduce it, something changed and you need to know.
#  @param excluded_samples_file Optional path to a one-accession-per-line file of
#                        samples to drop (an authors' erratum). If named and
#                        absent, preparation STOPS -- running the uncorrected
#                        cohort and calling it corrected is the worse failure.
#  @param analysis_role  "nested_cv" for a method-comparison cohort or
#                        "external_holdout" for a cohort that receives panels
#                        fixed in another dataset.
#  @param source_accession Dataset in which panels are selected before an
#                        external holdout is read. NULL for nested-CV cohorts.
#  @param locked_panel_file Panel manifest that must exist before an external
#                        holdout can be prepared.
#  @param notes          Everything a reader needs to interpret the entry.
validation_dataset <- function(accession, label, scale,
                               outcome,
                               disease_term, target_go_terms,
                               expected_n,
                               derive = NULL,
                               sample_filters = list(),
                               group_column = NULL,
                               confounders = character(0),
                               excluded_samples_file = NULL,
                               analysis_role = "nested_cv",
                               source_accession = NULL,
                               locked_panel_file = NULL,
                               notes = "") {
  if (!analysis_role %in% c("nested_cv", "external_holdout")) {
    stop("analysis_role must be 'nested_cv' or 'external_holdout'.")
  }
  if (analysis_role == "external_holdout" &&
      (is.null(source_accession) || is.null(locked_panel_file))) {
    stop("external_holdout datasets require source_accession and locked_panel_file.")
  }
  list(accession = accession, label = label, scale = scale,
       derive = derive, sample_filters = sample_filters, outcome = outcome,
       group_column = group_column, confounders = confounders,
       disease_term = disease_term, target_go_terms = target_go_terms,
       expected_n = expected_n,
       excluded_samples_file = excluded_samples_file,
       analysis_role = analysis_role, source_accession = source_accession,
       locked_panel_file = locked_panel_file,
       notes = notes)
}


# ------------------------------------------------------------------------------
#  Shared GO vocabulary
# ------------------------------------------------------------------------------
#
#  Terms reused across cohorts, spelled out once so a typo cannot silently make
#  two datasets score against different sets. Every ID is a GO Biological
#  Process term.

GO_IMMUNE_SYSTEM_PROCESS  <- "GO:0002376"
GO_IMMUNE_RESPONSE        <- "GO:0006955"
GO_DEFENSE_RESPONSE       <- "GO:0006952"
GO_INFLAMMATORY_RESPONSE  <- "GO:0006954"
GO_INNATE_IMMUNE_RESPONSE <- "GO:0045087"
GO_ADAPTIVE_IMMUNE_RESP   <- "GO:0002250"
GO_CYTOKINE_PRODUCTION    <- "GO:0001816"
GO_RESPONSE_TO_BACTERIUM  <- "GO:0009617"
GO_RESPONSE_TO_LPS        <- "GO:0032496"
GO_TYPE_I_IFN_SIGNALING   <- "GO:0060337"
GO_RESPONSE_TO_TYPE_II_IFN<- "GO:0034341"
GO_DEFENSE_GRAM_POSITIVE  <- "GO:0050830"
GO_TYPE_2_IMMUNE_RESPONSE <- "GO:0042092"
GO_TH17_IMMUNE_RESPONSE   <- "GO:0072538"
GO_NEUTROPHIL_CHEMOTAXIS  <- "GO:0030593"
GO_LYMPHOCYTE_CHEMOTAXIS  <- "GO:0048247"
GO_EPIDERMIS_DEVELOPMENT  <- "GO:0008544"
GO_KERATINOCYTE_DIFF      <- "GO:0030216"
GO_CELL_CYCLE             <- "GO:0007049"
GO_MITOTIC_CELL_CYCLE     <- "GO:0000278"
GO_CELL_PROLIFERATION     <- "GO:0008283"
GO_DNA_REPAIR             <- "GO:0006281"
GO_DNA_DAMAGE_RESPONSE    <- "GO:0006974"
GO_APOPTOTIC_PROCESS      <- "GO:0006915"
GO_RESPONSE_TO_DRUG       <- "GO:0042493"
GO_GI_EPITHELIUM_MAINT    <- "GO:0030277"
GO_CELLULAR_RESPONSE_TNF  <- "GO:0071356"


# ------------------------------------------------------------------------------
#  The five cohorts
# ------------------------------------------------------------------------------

validation_datasets <- list(

  # ============================================================================
  #  GSE65682 — Sepsis, MARS cohort, whole blood (Affymetrix U219)
  # ============================================================================
  #
  #  Contrast: 28-day mortality among sepsis patients (non-survivor = positive),
  #  NOT sepsis vs healthy. Sepsis vs healthy is trivially separable in whole
  #  blood and would measure nothing; mortality is the prediction problem the
  #  cohort exists for, and it is genuinely hard.
  #
  #  Healthy controls carry mortality_event_28days = NA and are therefore
  #  dropped by the outcome mapping, not by a filter -- the same mechanism that
  #  drops patients with no follow-up.
  GSE65682 = validation_dataset(
    accession = "GSE65682",
    label     = "Sepsis, MARS (whole blood) — 28-day mortality",
    scale     = "log_normalized",
    outcome   = list(
      column         = "ch_mortality_event_28days",
      negative       = c("0"),
      positive       = c("1"),
      negative_label = "survivor",
      positive_label = "non_survivor"
    ),
    group_column = NULL,
    # Sex and age are recorded and both are associated with sepsis mortality.
    # Only sex is residualised: the residualisation design matrix is built from
    # CATEGORICAL columns, and binning age would impose an arbitrary cut.
    confounders  = c("ch_gender"),
    disease_term = "sepsis",
    target_go_terms = c(
      GO_IMMUNE_SYSTEM_PROCESS,
      GO_IMMUNE_RESPONSE,
      GO_INFLAMMATORY_RESPONSE,
      GO_INNATE_IMMUNE_RESPONSE,
      GO_RESPONSE_TO_LPS,
      GO_RESPONSE_TO_BACTERIUM,
      GO_CYTOKINE_PRODUCTION
    ),
    expected_n = c(negative = 365, positive = 114),
    notes = paste(
      "802 samples in the series; only the ~479 sepsis patients with 28-day",
      "follow-up enter the analysis. ~114 events, so the positive class is",
      "~24% -- report balanced accuracy and MCC alongside AUC.",
      "The series matrix currently in data/GSE65682/raw is TRUNCATED (gzip",
      "reports 'unexpected end of file'); re-download before preparing."
    )
  ),

  # ============================================================================
  #  GSE69683 — Asthma, U-BIOPRED, peripheral blood (Affymetrix HT HG-U133+ PM)
  # ============================================================================
  #
  #  Contrast: severe vs moderate asthma, NON-SMOKERS ONLY.
  #
  #  The smoking arm is dropped rather than adjusted for. There are 88 severe
  #  smokers and ZERO moderate smokers, so smoking status is perfectly
  #  confounded with severity within the smoking stratum: any model could get
  #  those 88 right by detecting tobacco. Residualisation cannot fix a
  #  confounder with no overlap. Dropping them costs sample size and buys a
  #  contrast that means what it says.
  #
  #  Healthy controls are excluded by the outcome mapping for the same reason as
  #  in sepsis: healthy-vs-disease is the easy question.
  GSE69683 = validation_dataset(
    accession = "GSE69683",
    label     = "Asthma, U-BIOPRED (blood) — severe vs moderate",
    scale     = "log_normalized",
    sample_filters = list(
      # Non-smokers only. Note this keeps healthy non-smokers too; the outcome
      # mapping below drops them.
      list(column = "ch_cohort",
           keep   = c("Severe asthma, non-smoking",
                      "Moderate asthma, non-smoking",
                      "Healthy, non-smoking"))
    ),
    outcome = list(
      column         = "ch_cohort",
      negative       = c("Moderate asthma, non-smoking"),
      positive       = c("Severe asthma, non-smoking"),
      negative_label = "moderate",
      positive_label = "severe"
    ),
    group_column = NULL,
    confounders  = c("ch_gender"),
    disease_term = "asthma",
    target_go_terms = c(
      GO_IMMUNE_RESPONSE,
      GO_INFLAMMATORY_RESPONSE,
      GO_ADAPTIVE_IMMUNE_RESP,
      GO_TYPE_2_IMMUNE_RESPONSE,
      GO_CYTOKINE_PRODUCTION,
      GO_NEUTROPHIL_CHEMOTAXIS,
      GO_LYMPHOCYTE_CHEMOTAXIS
    ),
    expected_n = c(negative = 77, positive = 246),
    notes = paste(
      "Published AUC for this contrast is around 0.70, so this is a cohort",
      "where a real signal exists and methods can be separated -- unlike",
      "SOS-ALL. Class ratio is 3.2:1 toward severe."
    )
  ),

  # ============================================================================
  #  GSE13355 — Psoriasis, lesional vs uninvolved skin (Affymetrix U133 Plus 2.0)
  # ============================================================================
  #
  #  POSITIVE CONTROL, NOT A DISCRIMINATOR. Lesional versus uninvolved psoriatic
  #  skin is one of the largest effect sizes in human transcriptomics; every
  #  method should reach a near-ceiling AUC. It is included to show the pipeline
  #  detects a signal that is unambiguously there. A method that FAILS here is
  #  broken; a method that wins here has won nothing, and the paper must say so.
  #
  #  PAIRED: each of the 58 patients contributes one PP (involved) and one PN
  #  (uninvolved) biopsy. Splitting those two samples across a CV boundary lets
  #  the model memorise the patient rather than the lesion, which is worth many
  #  AUC points on its own. `group_column` forces both the outer folds and the
  #  stability subsamples to move whole patients.
  #
  #  The 64 NN normal-control samples are dropped: they come from different
  #  individuals and mixing them in would break the pairing.
  GSE13355 = validation_dataset(
    accession = "GSE13355",
    label     = "Psoriasis (skin) — lesional vs uninvolved [PAIRED]",
    scale     = "log_normalized",
    # Sample_title is "Individual_<patient>_<NN|PN|PP>_sample". The patient ID
    # and the biopsy type exist nowhere else in this series' metadata, so they
    # are parsed here rather than guessed downstream. The regex is anchored and
    # the parse is checked, because a silent NA in the group column would
    # degrade grouped CV back to ordinary CV without any error.
    derive = function(metadata) {
      parsed <- regmatches(
        metadata$Sample_title,
        regexec("^Individual_(.+)_(NN|PN|PP)_sample$", metadata$Sample_title)
      )
      metadata$patient_id <- vapply(parsed, function(p)
        if (length(p) == 3) p[2] else NA_character_, character(1))
      metadata$skin_type  <- vapply(parsed, function(p)
        if (length(p) == 3) p[3] else NA_character_, character(1))

      n_unparsed <- sum(is.na(metadata$patient_id))
      if (n_unparsed > 0) {
        stop(sprintf(paste0("GSE13355: %d Sample_title values did not match ",
                            "'Individual_<id>_<NN|PN|PP>_sample'. The patient ",
                            "grouping would be silently wrong; fix the parse ",
                            "before running."), n_unparsed))
      }
      metadata
    },
    sample_filters = list(
      list(column = "skin_type", keep = c("PN", "PP"))
    ),
    outcome = list(
      column         = "skin_type",
      negative       = c("PN"),
      positive       = c("PP"),
      negative_label = "uninvolved",
      positive_label = "lesional"
    ),
    group_column = "patient_id",
    # No residualisation. Patient is the obvious blocking factor and it IS
    # orthogonal to the outcome here (every patient contributes one of each
    # class), but encoding it costs 57 design columns against 116 samples. The
    # pairing is handled by grouped CV instead, which spends no degrees of
    # freedom. This series also has a known batch effect that no recorded
    # column identifies, which is a further reason to treat it as a control.
    confounders  = character(0),
    disease_term = "psoriasis",
    target_go_terms = c(
      GO_IMMUNE_RESPONSE,
      GO_INFLAMMATORY_RESPONSE,
      GO_INNATE_IMMUNE_RESPONSE,
      GO_ADAPTIVE_IMMUNE_RESP,
      GO_TH17_IMMUNE_RESPONSE,
      GO_EPIDERMIS_DEVELOPMENT,
      GO_KERATINOCYTE_DIFF
    ),
    expected_n = c(negative = 58, positive = 58),
    notes = paste(
      "58 complete PP/PN pairs verified in the metadata. Perfectly balanced.",
      "Expect near-ceiling AUC from everything including Random at large k --",
      "read AUC-minus-Random at matched k, never raw AUC."
    )
  ),

  # ============================================================================
  #  GSE107994 — Tuberculosis, Leicester cohort, whole blood RNA-seq
  # ============================================================================
  #
  #  Contrast: active TB vs latent infection (LTBI). Both arms are infected, so
  #  this asks the clinically real question -- who has progressed -- rather than
  #  infected vs uninfected. The 50 uninfected controls are dropped.
  #
  #  TWO pre-registered restrictions, both needed and both decided here:
  #
  #    * BASELINE VISIT ONLY. The series is longitudinal: some patients are
  #      sampled repeatedly and would otherwise appear in both the training and
  #      the test half of a split. At baseline every patient appears exactly
  #      once (verified: 160 baseline samples, 160 distinct patient IDs), so
  #      this removes the repeated-measures problem outright instead of
  #      managing it. `group_column` stays set as a belt-and-braces check.
  #
  #    * AUTHOR-FLAGGED OUTLIERS DROPPED (ch_outlier == "Yes"). This flag is the
  #      submitters' own QC call, made before and independently of anything in
  #      this benchmark, so honouring it is not a choice we are making about our
  #      own results.
  #
  #  Note this yields 42 vs 52 = 94 samples, NOT the 125 quoted when the cohort
  #  was picked. That 125 counted every timepoint of every TB/LTBI patient,
  #  i.e. repeated measures of the same people. 94 independent patients is the
  #  smaller and correct number.
  GSE107994 = validation_dataset(
    accession = "GSE107994",
    label     = "Tuberculosis, Leicester (blood) — active vs latent",
    scale     = "counts",
    sample_filters = list(
      list(column = "ch_timepoint_months", keep = c("Baseline")),
      list(column = "ch_outlier",          keep = c("No"))
    ),
    outcome = list(
      column         = "ch_group",
      # LTBI_Progressor is latent AT BASELINE -- progression happens later --
      # so it belongs in the latent arm. Its 8 baseline samples are far too few
      # to be a third class, and treating a not-yet-progressed patient as
      # active TB would label the outcome with information from the future.
      negative       = c("LTBI", "LTBI_Progressor"),
      positive       = c("Active_TB"),
      negative_label = "latent",
      positive_label = "active"
    ),
    group_column = "ch_patient_id",
    confounders  = c("ch_gender", "ch_ethnicity"),
    disease_term = "tuberculosis",
    target_go_terms = c(
      GO_IMMUNE_SYSTEM_PROCESS,
      GO_IMMUNE_RESPONSE,
      GO_INNATE_IMMUNE_RESPONSE,
      GO_INFLAMMATORY_RESPONSE,
      GO_TYPE_I_IFN_SIGNALING,
      GO_RESPONSE_TO_TYPE_II_IFN,
      GO_RESPONSE_TO_BACTERIUM
    ),
    expected_n = c(negative = 52, positive = 42),
    notes = paste(
      paste(
        "RNA-seq counts; the benchmark applies the fixed abundance rule,",
        "TMM and log-CPM within each training fold."
      ),
      "Expression ships in a supplementary .xlsx, not the series matrix:",
      "fetch GSE107994_Raw_counts_Leicester_with_progressor_longitudinal.xlsx",
      "and convert it to data/GSE107994/expression.tsv before preparing.",
      "The type-I interferon signature is the established active-TB marker, so",
      "GO:0060337 is the term to watch in the biology pillar."
    )
  ),

  # ============================================================================
  #  GSE20194 — Breast cancer, MAQC-II, neoadjuvant chemotherapy
  # ============================================================================
  #
  #  Contrast: pathological complete response (pCR) vs residual disease (RD)
  #  after neoadjuvant chemotherapy.
  #
  #  Chosen because this series exists to BENCHMARK PREDICTION METHODS -- it is
  #  the MAQC-II reference breast-cancer set. That makes it the least
  #  cherry-pickable venue available for the claim, and a null here would be as
  #  informative as a win.
  #
  #  No confounders are residualised. ER status is by far the strongest clinical
  #  covariate of pCR, but it is BIOLOGY, not nuisance: ER-negative tumours
  #  genuinely respond more often. Residualising it would remove real signal the
  #  method is supposed to find. Nothing in the metadata identifies a batch or
  #  site, so there is no nuisance term to remove.
  GSE20194 = validation_dataset(
    accession = "GSE20194",
    label     = "Breast cancer, MAQC-II — pCR vs residual disease",
    scale     = "log_normalized",
    derive = function(metadata) {
      metadata$response <- extract_characteristic(metadata, "pcr_vs_rd")
      n_missing <- sum(is.na(metadata$response))
      if (n_missing > 0) {
        stop(sprintf("GSE20194: %d samples have no 'pcr_vs_rd' characteristic.",
                     n_missing))
      }
      metadata
    },
    outcome = list(
      column         = "response",
      negative       = c("RD"),
      positive       = c("pCR"),
      negative_label = "residual_disease",
      positive_label = "pCR"
    ),
    group_column = NULL,
    confounders  = character(0),
    disease_term = "breast carcinoma",
    target_go_terms = c(
      GO_CELL_CYCLE, GO_MITOTIC_CELL_CYCLE, GO_CELL_PROLIFERATION,
      GO_DNA_REPAIR, GO_DNA_DAMAGE_RESPONSE, GO_APOPTOTIC_PROCESS,
      GO_RESPONSE_TO_DRUG
    ),
    expected_n = c(negative = 222, positive = 56),
    notes = paste(
      "Characteristic fields are NOT in a consistent order across samples --",
      "the positional parse put response values into the ER-status column, so",
      "the outcome is re-read by key. Class ratio 4:1 toward residual disease."
    )
  ),

  # ============================================================================
  #  GSE25066 — Breast cancer, neoadjuvant taxane-anthracycline
  # ============================================================================
  #
  #  The same pCR vs RD contrast as GSE20194, in a larger independent cohort.
  #  Run together they give a replication rather than a single result.
  #
  #  `source` IS residualised here: this cohort pools four sites (MDACC, ISPY,
  #  LBJ/IN/GEI, USO), which is a genuine batch effect rather than biology. Four
  #  levels, so the design stays small -- unlike the 14-level ethnicity variable
  #  that broke GSE107994.
  GSE25066 = validation_dataset(
    accession = "GSE25066",
    label     = "Breast cancer, neoadjuvant chemo (n=488) — pCR vs residual disease",
    scale     = "log_normalized",
    derive = function(metadata) {
      metadata$response <- extract_characteristic(metadata,
                              "pathologic_response_pcr_rd")
      metadata$site     <- extract_characteristic(metadata, "source")
      metadata$site[is.na(metadata$site) | !nzchar(metadata$site)] <- "unknown"
      metadata
    },
    outcome = list(
      column         = "response",
      # "NA" is a literal string in this series, not a missing value; it is
      # listed nowhere below, so those samples are dropped by the mapping.
      negative       = c("RD"),
      positive       = c("pCR"),
      negative_label = "residual_disease",
      positive_label = "pCR"
    ),
    group_column = NULL,
    confounders  = c("site"),
    disease_term = "breast carcinoma",
    target_go_terms = c(
      GO_CELL_CYCLE, GO_MITOTIC_CELL_CYCLE, GO_CELL_PROLIFERATION,
      GO_DNA_REPAIR, GO_DNA_DAMAGE_RESPONSE, GO_APOPTOTIC_PROCESS,
      GO_RESPONSE_TO_DRUG
    ),
    expected_n = c(negative = 389, positive = 99),
    notes = paste(
      "508 samples in the series; 20 carry a literal 'NA' response and are",
      "dropped, leaving 488. Same characteristic-order problem as GSE20194, so",
      "the outcome and site are re-read by key."
    )
  ),

  # ============================================================================
  #  GSE16879 — IBD, response to infliximab (mucosal biopsy)
  # ============================================================================
  #
  #  Contrast: responder vs non-responder to infliximab, predicted from
  #  PRE-TREATMENT mucosa.
  #
  #  Chosen as the closest structural match to IMvigor210 -- "will this patient
  #  respond to this biologic" -- and because anti-TNF response has resisted
  #  prediction for years, so the signal should be weak. That is the regime the
  #  method has looked strongest in.
  #
  #  BE HONEST ABOUT WHY THIS COHORT WAS PICKED: it was selected AFTER seeing
  #  that GeneSelectR wins on low-signal problems and loses on high-signal ones.
  #  It is a hypothesis test of that pattern, NOT pre-registered validation, and
  #  it must not be reported as the latter.
  #
  #  THE PRE-TREATMENT FILTER IS NOT OPTIONAL. The series contains paired
  #  before/after biopsies from the same patients. Post-treatment mucosa reflects
  #  whether the drug worked, so including it would let a classifier read the
  #  answer off the outcome itself. Only the 61 "Before first infliximab
  #  treatment" samples are used, which also means each patient appears once and
  #  no grouping is needed.
  GSE16879 = validation_dataset(
    accession = "GSE16879",
    label     = "IBD, infliximab response (pre-treatment mucosa)",
    # MAS5 output on a linear scale (max ~68,000), NOT logged -- prepare applies
    # log2(x + 1). The declared value is checked against the matrix.
    scale     = "linear_normalized",
    sample_filters = list(
      list(column = "ch_before_or_after_first_infliximab_treatment",
           keep   = c("Before first infliximab treatment"))
    ),
    outcome = list(
      column         = "ch_response_to_infliximab",
      negative       = c("No"),
      positive       = c("Yes"),
      negative_label = "non_responder",
      positive_label = "responder"
    ),
    group_column = NULL,
    # Tissue only. Ileal and colonic mucosa differ enormously for reasons that
    # have nothing to do with drug response, so it is nuisance. Disease (CD vs
    # UC) is NOT residualised: response rates genuinely differ between them
    # (54% vs 33% here), so it is case-mix biology rather than noise, and
    # removing it would strip signal. One design column keeps this safe at
    # n = 61.
    confounders  = c("ch_tissue"),
    disease_term = "inflammatory bowel disease",
    target_go_terms = c(
      GO_INFLAMMATORY_RESPONSE, GO_IMMUNE_SYSTEM_PROCESS, GO_IMMUNE_RESPONSE,
      GO_INNATE_IMMUNE_RESPONSE, GO_CYTOKINE_PRODUCTION,
      GO_CELLULAR_RESPONSE_TNF, GO_GI_EPITHELIUM_MAINT
    ),
    expected_n = c(negative = 33, positive = 28),
    notes = paste(
      "61 pre-treatment samples, 28 responders / 33 non-responders -- the most",
      "balanced cohort in this registry. Small: 5-fold CV leaves ~12 samples",
      "per test fold, so per-fold AUC will be noisy and the error bars wide.",
      "Mixed CD (37) and UC (24)."
    )
  ),

  # ============================================================================
  #  GSE57945 — Crohn's disease, RISK cohort, ileal biopsies (RNA-seq)
  # ============================================================================
  #
  #  Contrast: treatment-naive paediatric CD vs non-IBD control ileum. UC is
  #  dropped -- it is a third disease, not a severity grade of the same one.
  #
  #  THE AUTHORS ISSUED A CORRECTION excluding a set of mislabelled samples.
  #  `excluded_samples_file` names the file that must list them. If the file is
  #  absent, preparation STOPS. It would be trivial to run the uncorrected 218
  #  vs 42 and describe it as the corrected cohort; that is exactly the failure
  #  this stop exists to prevent. Populate the file from the published erratum
  #  (one GSM accession per line) or drop the dataset from the report.
  #
  #  `scale` is declared "counts" but the series' own documentation is
  #  ambiguous and RPKM is also distributed. validation_prepare.R verifies the
  #  matrix against this declaration and errors on a mismatch -- TMM-normalising
  #  RPKM is wrong rather than merely suboptimal, and nothing downstream would
  #  complain about it.
  GSE57945 = validation_dataset(
    accession = "GSE57945",
    label     = "Crohn's disease, RISK (ileum) — CD vs non-IBD",
    scale     = "counts",
    outcome = list(
      column         = "ch_diagnosis",
      # "Not IBD" and "not IBD" both occur in the metadata; both are listed
      # rather than case-folded, so a third spelling errors instead of being
      # quietly absorbed.
      negative       = c("Not IBD", "not IBD"),
      positive       = c("CD"),
      negative_label = "not_IBD",
      positive_label = "CD"
    ),
    group_column = NULL,
    confounders  = c("ch_sex"),
    excluded_samples_file = "data/GSE57945/authors_correction_excluded.txt",
    disease_term = "Crohn's disease",
    target_go_terms = c(
      GO_IMMUNE_SYSTEM_PROCESS,
      GO_IMMUNE_RESPONSE,
      GO_INFLAMMATORY_RESPONSE,
      GO_INNATE_IMMUNE_RESPONSE,
      GO_ADAPTIVE_IMMUNE_RESP,
      GO_RESPONSE_TO_BACTERIUM,
      GO_GI_EPITHELIUM_MAINT
    ),
    # Uncorrected the series holds 218 CD and 42 non-IBD. The erratum removes
    # 26 CD and 2 controls, giving the numbers below. If the prepared cohort
    # does not match, the exclusion list is not the one the erratum specifies.
    expected_n = c(negative = 40, positive = 192),
    notes = paste(
      "RNA-seq counts, ileal biopsy, treatment-naive at sampling so no drug",
      "effect confounds the contrast. Heavily imbalanced (~5:1 toward CD).",
      "Expression ships in a supplementary file (GSE57945_RAW.tar or the RPKM",
      "table), not the series matrix. Use the COUNTS, not the RPKM: the RPKM",
      "table will fail the scale check, which is the intended behaviour."
    )
  ),

  # ============================================================================
  #  GSE92415 — Ulcerative colitis, golimumab induction response
  # ============================================================================
  #
  #  Contrast: week-6 clinical response to golimumab predicted from the
  #  pretreatment colonic mucosal biopsy. GEO contains healthy samples, placebo
  #  samples, and week-6 biopsies. The two filters below retain the 59 baseline
  #  golimumab-treated participants with a recorded response.
  GSE92415 = validation_dataset(
    accession = "GSE92415",
    label     = "Ulcerative colitis, golimumab — week-6 response",
    scale     = "log_normalized",
    derive = function(metadata) {
      metadata$subject       <- extract_characteristic(metadata, "subject")
      metadata$treatment     <- extract_characteristic(metadata, "treatment")
      metadata$visit         <- extract_characteristic(metadata, "visit")
      metadata$wk6_response  <- extract_characteristic(metadata, "wk6response")
      metadata
    },
    sample_filters = list(
      list(column = "treatment", keep = c("golimumab")),
      list(column = "visit", keep = c("Week 0"))
    ),
    outcome = list(
      column         = "wk6_response",
      negative       = c("No"),
      positive       = c("Yes"),
      negative_label = "non_responder",
      positive_label = "responder"
    ),
    group_column = "subject",
    confounders  = character(0),
    disease_term = "ulcerative colitis",
    target_go_terms = c(
      GO_INFLAMMATORY_RESPONSE,
      GO_IMMUNE_SYSTEM_PROCESS,
      GO_IMMUNE_RESPONSE,
      GO_ADAPTIVE_IMMUNE_RESP,
      GO_CYTOKINE_PRODUCTION,
      GO_TH17_IMMUNE_RESPONSE,
      GO_CELLULAR_RESPONSE_TNF,
      GO_GI_EPITHELIUM_MAINT
    ),
    expected_n = c(negative = 27, positive = 32),
    notes = paste(
      "GEO contains 87 baseline UC biopsies. Restricting to the golimumab arm",
      "gives 32 week-6 responders and 27 non-responders. Every retained",
      "participant contributes one baseline biopsy. This replaces the mixed",
      "disease and tissue analysis in GSE16879."
    )
  ),

  # ============================================================================
  #  GSE206285 — Ulcerative colitis, ustekinumab induction response
  # ============================================================================
  #
  #  Contrast: clinical remission at week 8 among participants randomized to
  #  ustekinumab, predicted from the week-0 sigmoid-colon biopsy. Healthy and
  #  placebo samples are excluded. All weight-based 6 mg/kg dose labels belong
  #  to the same randomized regimen and are retained with the 130 mg regimen.
  GSE206285 = validation_dataset(
    accession = "GSE206285",
    label     = "Ulcerative colitis, ustekinumab — week-8 remission",
    scale     = "log_normalized",
    derive = function(metadata) {
      metadata$donor_id <- extract_characteristic(metadata, "donor id")
      metadata$visit    <- extract_characteristic(metadata, "visit")
      metadata$diagnosis <- extract_characteristic(metadata, "diagnosis")
      metadata$treatment <- extract_characteristic(metadata, "treatment")
      metadata$clinical_remission <- extract_characteristic(
        metadata, "clinical remission at week 8"
      )
      metadata$treatment_arm <- ifelse(
        grepl("^Ustekinumab 6 mg/kg", metadata$treatment),
        "ustekinumab_6mgkg",
        ifelse(metadata$treatment == "Ustekinumab 130 mg IV",
               "ustekinumab_130mg", NA_character_)
      )
      metadata
    },
    sample_filters = list(
      list(column = "diagnosis", keep = c("ulcerative colitis")),
      list(column = "visit", keep = c("WEEK I-0")),
      list(column = "treatment_arm",
           keep = c("ustekinumab_130mg", "ustekinumab_6mgkg"))
    ),
    outcome = list(
      column         = "clinical_remission",
      negative       = c("N"),
      positive       = c("Y"),
      negative_label = "no_remission",
      positive_label = "remission"
    ),
    group_column = "donor_id",
    confounders  = character(0),
    disease_term = "ulcerative colitis",
    target_go_terms = c(
      GO_INFLAMMATORY_RESPONSE,
      GO_IMMUNE_SYSTEM_PROCESS,
      GO_IMMUNE_RESPONSE,
      GO_ADAPTIVE_IMMUNE_RESP,
      GO_CYTOKINE_PRODUCTION,
      GO_TH17_IMMUNE_RESPONSE,
      GO_GI_EPITHELIUM_MAINT
    ),
    expected_n = c(negative = 315, positive = 49),
    notes = paste(
      "The analysis includes all 364 ustekinumab-treated baseline biopsies:",
      "49 participants achieved clinical remission at week 8 and 315 did not.",
      "Treatment arm is retained in the derived metadata for auditing.",
      "Baseline expression precedes assigned treatment. Treatment-stratified",
      "results can be reported as a sensitivity analysis."
    )
  ),

  # ============================================================================
  #  GSE101794 — Treatment-naive paediatric Crohn disease
  # ============================================================================
  #
  #  Positive control replacing GSE57945. GEO supplies Kallisto TPM, and the
  #  source builder writes log2(TPM + 1). The array/count normalization branch
  #  therefore receives a verified log-normalized matrix and does not use TMM.
  GSE101794 = validation_dataset(
    accession = "GSE101794",
    label     = "Paediatric Crohn disease (ileum) — CD vs non-IBD",
    scale     = "log_normalized",
    outcome = list(
      column         = "diagnosis",
      negative       = c("Non-IBD"),
      positive       = c("CD"),
      negative_label = "non_IBD",
      positive_label = "CD"
    ),
    group_column = NULL,
    confounders  = character(0),
    disease_term = "Crohn's disease",
    target_go_terms = c(
      GO_IMMUNE_SYSTEM_PROCESS,
      GO_IMMUNE_RESPONSE,
      GO_INFLAMMATORY_RESPONSE,
      GO_INNATE_IMMUNE_RESPONSE,
      GO_ADAPTIVE_IMMUNE_RESP,
      GO_RESPONSE_TO_BACTERIUM,
      GO_GI_EPITHELIUM_MAINT
    ),
    expected_n = c(negative = 50, positive = 254),
    notes = paste(
      "Treatment-naive ileal biopsies quantified as TPM. Run",
      "benchmarks/gse101794_build_expression.R before validation_prepare.R.",
      "This is a high-signal positive control and is excluded from the primary",
      "difficult-endpoint average."
    )
  ),

  # ============================================================================
  #  GSE19442 — South African TB external holdout
  # ============================================================================
  #
  #  Panels are selected and fixed using GSE107994. GSE19442 is then read once
  #  for external evaluation. Running nested feature-selection CV on this
  #  cohort would consume the holdout and is blocked in validation_benchmark.R.
  GSE19442 = validation_dataset(
    accession = "GSE19442",
    label     = "Tuberculosis, South Africa — external active vs latent holdout",
    scale     = "linear_background_corrected",
    derive = function(metadata) {
      sample_title <- as.character(metadata$Sample_title)
      matched <- regexec("^(PTB|LTB)_SA_val[0-9]+$", sample_title)
      parsed <- regmatches(sample_title, matched)
      metadata$tb_group <- vapply(
        parsed,
        function(x) if (length(x) == 2) x[2] else NA_character_,
        character(1)
      )
      if (anyNA(metadata$tb_group)) {
        stop("GSE19442: one or more sample titles do not encode PTB/LTB.")
      }
      metadata
    },
    outcome = list(
      column         = "tb_group",
      negative       = c("LTB"),
      positive       = c("PTB"),
      negative_label = "latent",
      positive_label = "active"
    ),
    group_column = NULL,
    confounders  = character(0),
    disease_term = "tuberculosis",
    target_go_terms = c(
      GO_IMMUNE_SYSTEM_PROCESS,
      GO_IMMUNE_RESPONSE,
      GO_INNATE_IMMUNE_RESPONSE,
      GO_INFLAMMATORY_RESPONSE,
      GO_TYPE_I_IFN_SIGNALING,
      GO_RESPONSE_TO_TYPE_II_IFN,
      GO_RESPONSE_TO_BACTERIUM
    ),
    expected_n = c(negative = 31, positive = 20),
    analysis_role = "external_holdout",
    source_accession = "GSE107994",
    locked_panel_file = "locked_panels/GSE107994_to_GSE19442.csv",
    notes = paste(
      "All 51 participants were sampled before antimycobacterial treatment.",
      "Panels must be locked in GSE107994 before this cohort is prepared or",
      "scored. This cohort is excluded from nested-CV method ranking."
    )
  ),

  # ============================================================================
  #  GSE91061 — Melanoma, anti-PD-1 (nivolumab), tumour RNA-seq (Riaz 2017)
  # ============================================================================
  #
  #  Added 2026-08-13 to test the difficulty hypothesis in prereg_difficulty.md,
  #  which predicts a positive GeneSelectR margin wherever the best competitor
  #  scores below 0.80 AUC. Checkpoint-blockade response was predicted hard
  #  before the run, on the same grounds as IMvigor210.
  #
  #  PRE-TREATMENT SAMPLES ONLY. The series pairs each patient's pre-treatment
  #  and on-treatment biopsy (51 Pre, 58 On). On-treatment expression reflects
  #  the drug response that is being predicted, and keeping both timepoints
  #  would put the same patient in the training and test half of a split. At
  #  the Pre visit each patient appears once. `group_column` is set as a check.
  #
  #  RESPONSE CODING FOLLOWS IMvigor210. That benchmark scores CR/PR as
  #  responder and SD/PD as non-responder, so PRCR is positive here and SD and
  #  PD are both negative. UNK is listed nowhere below and is therefore dropped.
  #  Coding SD as a responder, or discarding it, would each make this cohort
  #  incomparable to IMvigor210, which is the cohort it is paired with.
  #  Pre-treatment counts: PRCR 10, SD 16, PD 23, UNK 2 -> 10 vs 39, n=49.
  #  Ten events across 5 folds means two per test fold, so per-fold AUC will be
  #  very noisy; this is the smallest cohort in the registry.
  GSE91061 = validation_dataset(
    accession = "GSE91061",
    label     = "Melanoma, anti-PD-1 (n=49) — responder vs non-responder",
    scale     = "counts",
    derive = function(metadata) {
      metadata$visit    <- extract_characteristic(metadata,
                              "visit \\(pre or on treatment\\)")
      metadata$response <- extract_characteristic(metadata, "response")
      # Patient ID is the leading Pt<N> field of the sample title; it is the
      # only patient key the series exposes.
      metadata$patient  <- sub("_.*$", "", metadata$Sample_title)
      metadata
    },
    sample_filters = list(
      list(column = "visit", keep = c("Pre"))
    ),
    outcome = list(
      column         = "response",
      negative       = c("SD", "PD"),
      positive       = c("PRCR"),
      negative_label = "non_responder",
      positive_label = "responder"
    ),
    group_column = "patient",
    confounders  = NULL,
    disease_term = "melanoma",
    target_go_terms = c(
      GO_IMMUNE_SYSTEM_PROCESS,
      GO_IMMUNE_RESPONSE,
      GO_ADAPTIVE_IMMUNE_RESP,
      GO_INFLAMMATORY_RESPONSE,
      GO_CYTOKINE_PRODUCTION,
      GO_RESPONSE_TO_TYPE_II_IFN,
      GO_LYMPHOCYTE_CHEMOTAXIS
    ),
    expected_n = c(negative = 39, positive = 10),
    notes = paste(
      paste(
        "RNA-seq counts; the benchmark applies the fixed abundance rule,",
        "TMM and log-CPM within each training fold. Expression ships"
      ),
      "in a supplementary file (GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz),",
      "not the series matrix; its columns are sample TITLES and its rows are",
      "Entrez IDs, both translated when data/GSE91061/expression.tsv is built.",
      "22068 of 22187 Entrez IDs mapped to a symbol, no duplicate symbols."
    )
  )
)


# ------------------------------------------------------------------------------
#  Lookup helper
# ------------------------------------------------------------------------------
#
#  Fails with the list of valid accessions rather than returning NULL, so a
#  typo on the command line stops immediately instead of one screen later.
get_validation_dataset <- function(accession) {
  if (!accession %in% names(validation_datasets)) {
    stop(sprintf("Unknown dataset '%s'. Registered: %s",
                 accession, paste(names(validation_datasets), collapse = ", ")))
  }
  validation_datasets[[accession]]
}


# Prints the plan for one dataset. Called by both scripts at startup so every
# log begins with the pre-registered specification the run is executing.
print_validation_dataset <- function(dataset) {
  cat(sprintf("Dataset:        %s (%s)\n", dataset$label, dataset$accession))
  cat(sprintf("  scale:        %s\n", dataset$scale))
  cat(sprintf("  contrast:     %s (positive) vs %s (negative), column '%s'\n",
              dataset$outcome$positive_label, dataset$outcome$negative_label,
              dataset$outcome$column))
  cat(sprintf("  expected n:   %d %s / %d %s\n",
              dataset$expected_n[["positive"]], dataset$outcome$positive_label,
              dataset$expected_n[["negative"]], dataset$outcome$negative_label))
  cat(sprintf("  grouping:     %s\n", dataset$group_column %||% "none (independent samples)"))
  cat(sprintf("  confounders:  %s\n",
              if (length(dataset$confounders) == 0) "none"
              else paste(dataset$confounders, collapse = ", ")))
  cat(sprintf("  role:         %s\n", dataset$analysis_role))
  if (dataset$analysis_role == "external_holdout") {
    cat(sprintf("  source:       %s (panels fixed before holdout access)\n",
                dataset$source_accession))
    cat(sprintf("  panel file:   %s\n", dataset$locked_panel_file))
  }
  cat(sprintf("  disease term: %s   [PRE-REGISTERED]\n", dataset$disease_term))
  cat(sprintf("  GO terms:     %s   [PRE-REGISTERED]\n",
              paste(dataset$target_go_terms, collapse = ", ")))
  if (!is.null(dataset$excluded_samples_file)) {
    cat(sprintf("  exclusions:   %s\n", dataset$excluded_samples_file))
  }
  if (nzchar(dataset$notes)) {
    cat(sprintf("  notes:        %s\n", dataset$notes))
  }
  invisible(NULL)
}
