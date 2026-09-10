config_path <- file.path("analysis", "config.R")

analysis_config <- list(
    package_ref = "621be0c1",
    results_root = file.path("redesign", "results_corrected"),
    development_dataset = "sosall",
    validation_datasets = c(
        "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683",
        "imvigor210"
    ),
    exploratory_datasets = c("GSE16879", "GSE91061", "GSE92415", "GSE206285"),
    repeats = 3L,
    folds = 5L,
    gene_set_sizes = c(10L, 20L, 50L, 100L, 200L, 500L),
    candidate_genes = 2000L,
    fit_resamples = 50L,
    fit_folds = 5L,
    outcome_permutations = 20L,
    null_fits = 20L,
    seed = 42L,
    workers = as.integer(Sys.getenv("GENESELECTR_N_CORES", "2"))
)

biology_older7_config <- data.frame(
    dataset = c(
        "GSE101794", "GSE107994", "GSE13355", "GSE65682", "GSE69683",
        "imvigor210", "sosall"
    ),
    disease_label = c(
        "Crohn's disease", "Tuberculosis", "Psoriasis", "Sepsis", "Asthma",
        "Urinary bladder cancer treated with anti-PD-L1", "Atopic dermatitis"
    ),
    ontology_ids = c(
        "MONDO_0005265", "MONDO_0018076", "MONDO_0005083", "HP_0100806",
        "MONDO_0004979", "MONDO_0004986", "MONDO_0004980"
    ),
    association_files = file.path(
        "data", "r_user_cache", "R", "GeneSelectR",
        paste0("ot_seeds_", c(
            "MONDO0005265", "MONDO0018076", "MONDO0005083", "HP0100806",
            "MONDO0004979", "MONDO0004986", "MONDO0004980"
        ), "_n100_s0.1.rds")
    ),
    association_rule = "single disease score",
    target_go_terms = c(
        "GO:0002376;GO:0006955;GO:0006954;GO:0045087;GO:0002250;GO:0009617;GO:0030277",
        "GO:0002376;GO:0006955;GO:0045087;GO:0006954;GO:0060337;GO:0034341;GO:0009617",
        "GO:0006955;GO:0006954;GO:0045087;GO:0002250;GO:0072538;GO:0008544;GO:0030216",
        "GO:0002376;GO:0006955;GO:0006954;GO:0045087;GO:0032496;GO:0009617;GO:0001816",
        "GO:0006955;GO:0006954;GO:0002250;GO:0042092;GO:0001816;GO:0030593;GO:0048247",
        "GO:0006955;GO:0002376;GO:0045087;GO:0002250;GO:0042110;GO:0050863;GO:0007179;GO:0071559;GO:0002682;GO:0050776",
        "GO:0006955;GO:0002376;GO:0006952;GO:0002250;GO:0045087;GO:0002682;GO:0050776"
    ),
    stringsAsFactors = FALSE
)

biology_external_config <- data.frame(
    dataset = c("GSE16879", "GSE91061", "GSE92415", "GSE206285"),
    disease_label = c(
        "Inflammatory bowel disease (Crohn disease and ulcerative colitis)",
        "Melanoma", "Ulcerative colitis", "Ulcerative colitis"
    ),
    ontology_ids = c(
        "MONDO_0005265;MONDO_0005101", "MONDO_0005105", "MONDO_0005101",
        "MONDO_0005101"
    ),
    association_files = c(
        paste(
            file.path("data", "r_user_cache", "R", "GeneSelectR",
                      "ot_seeds_MONDO0005265_n100_s0.1.rds"),
            file.path("data", "r_user_cache", "R", "GeneSelectR",
                      "ot_seeds_MONDO0005101_n100_s0.1.rds"),
            sep = ";"
        ),
        file.path("data", "r_user_cache", "R", "GeneSelectR",
                  "ot_seeds_MONDO0005105_n100_s0.1.rds"),
        file.path("data", "r_user_cache", "R", "GeneSelectR",
                  "ot_seeds_MONDO0005101_n100_s0.1.rds"),
        file.path("data", "r_user_cache", "R", "GeneSelectR",
                  "ot_seeds_MONDO0005101_n100_s0.1.rds")
    ),
    association_rule = c(
        "maximum score across the two diseases", rep("single disease score", 3)
    ),
    stringsAsFactors = FALSE
)

biology_disease_config <- data.frame(
    dataset = c(
        "GSE65682", "GSE69683", "GSE13355", "GSE107994", "GSE101794",
        "imvigor210", "sosall"
    ),
    disease_label = c(
        "Sepsis", "Asthma", "Psoriasis", "Tuberculosis", "Crohn's disease",
        "Urinary bladder cancer", "Atopic dermatitis"
    ),
    ontology_id = c(
        "HP_0100806", "MONDO_0004979", "MONDO_0005083", "MONDO_0018076",
        "MONDO_0005265", "MONDO_0004986", "MONDO_0004980"
    ),
    association_file = file.path(
        "cache", "opentargets",
        paste0(c(
            "HP_0100806", "MONDO_0004979", "MONDO_0005083",
            "MONDO_0018076", "MONDO_0005265", "MONDO_0004986",
            "MONDO_0004980"
        ), ".tsv")
    ),
    stringsAsFactors = FALSE
)

get_biology_config <- function(name = c("older7", "external", "disease")) {
    name <- match.arg(name)
    switch(
        name,
        older7 = biology_older7_config,
        external = biology_external_config,
        disease = biology_disease_config
    )
}
