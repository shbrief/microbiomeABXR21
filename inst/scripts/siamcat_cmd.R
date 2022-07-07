## Required inputs: studyName, abx
## studyName = "VincentC_2016"
## abx = "cephalosporins"



# ## Metadata
# meta <- curatedMetagenomicData::sampleMetadata
# submeta <- meta %>% filter(study_name == studyName)
#
# ## Download data
# dat <- returnSamples(submeta, "relative_abundance")
# assay(dat) <- assay(dat)/100 # percentage to decimals
#
# ## Clean ABX info
# colData(dat)$abx <- colData(dat)$antibiotics_family
# colData(dat)$abx[which(is.na(colData(dat)$abx))] <- "control"
# abx_ls <- c("carbapenems", "cephalosporins", "fluoroquinolones",
#             "glycopeptides", "macrolides", "nitroimidazoles",
#             "penicillin", "sulfonamides")
# for (abx_name in abx_ls) {
#     ind <- which(colData(dat)$antibiotics_family %>% str_detect(., abx_name))
#     colData(dat)$abx[ind] <- abx_name
# }

## Subset data
dat <- subset(allData, select = colData(allData)$study_name == studyName)

## Convert to phyloseq object
phylo <- makePhyloseqFromTreeSummarizedExperiment(dat, abund_values = "relative_abundance")

## Extract sample and relative abundance tables <<<<<<<<<<<<<<<<<<<<<<< REVISE!!
sample.phylo <- sample_data(phylo)
feat.phylo <- otu_table(phylo)

## Set variable of interest
var_lab <- create.label(meta = sample.phylo,
                        label = "abx",
                        case = abx)

## Create SIAMCAT object
sc.obj <- siamcat(feat = feat.phylo,
                  label = var_lab,
                  meta = sample.phylo)

## Filter by abundance
sc.filt <- filter.features(sc.obj,
                           filter.method = "abundance",
                           cutoff = 0.001)

# ## Check association
# sc.assc <- check.associations(sc.filt, test = "lm", feature.type = "filtered")

# ## Check confounders and save the result in PDF
# fname <- paste0(studyName, "_confounder_", abx)
# if (length(metadata(dat)) != 0) {
#     fname <- paste0(fname, "_", metadata(dat))}
# sc.conf <- check.confounders(sc.filt,
#                              fn.plot = paste0(fname, ".pdf"),
#                              meta.in = NULL,
#                              feature.type = "filtered",
#                              verbose = 1)

## Normalization
sc.norm <- normalize.features(sc.filt, norm.method = "rank.unit")

## Split tranining and test datasets
sc.split <- create.data.split(sc.norm,
                              num.folds = 6,
                              num.resample = 2,
                              inseparable = "subject_id")

## Random Forest
system.time(sc.mod <- train.model(sc.split,
                                  method = "randomForest")) # this takes 5~10 min
sc.pred <- make.predictions(sc.mod)
pred_matrix <- pred_matrix(sc.pred)

## Model evaluation plot as PDF
fname <- paste0(studyName, "_evaluation_", abx)
if (length(metadata(dat)) != 0) {
    fname <- paste0(fname, "_", metadata(dat))}

sc.eval <- evaluate.predictions(sc.pred)
saveRDS(sc.eval, file = paste0(fname, ".rds")) # Save the model evaluation
model.evaluation.plot(sc.eval,
                      fn.plot = paste0(fname, ".pdf"))

## Export final model interpretation plot as PDF
fname <- paste0(studyName, "_interpretation_", abx)
if (length(metadata(dat)) != 0) {
    fname <- paste0(fname, "_", metadata(dat))}
model.interpretation.plot(sc.eval,
                          fn.plot = paste0(fname, ".pdf"),
                          consens.thres = 0.01,
                          limits = c(-3, 3),
                          heatmap.type = "zscore")
