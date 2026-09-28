#!/usr/bin/env Rscript
repo <- normalizePath(Sys.getenv('GENESELECTR_BENCHMARK_REPO', getwd()),mustWork=TRUE)
commit <- 'a840a6bf6ab36c87ac69ed65a5d506b6f29de21d'
out <- Sys.getenv('GENESELECTR_STABL_OUT',file.path(repo,'redesign/stabl_comparison'))
scratch <- tempdir()
datasets <- c('GSE101794','GSE107994','GSE13355','GSE65682','GSE69683','imvigor210','sosall')
methods <- c('GS_full_ungrouped','Stabl')
sizes <- c(10L,20L,50L,100L,200L,500L)
extract <- function(path,dest) {
  status <- system2('git',c('-C',shQuote(repo),'show',
    shQuote(paste0(commit,':',path))),stdout=dest)
  if(!identical(status,0L))stop('Missing ',path)
}
new <- read.csv(file.path(out,'evaluation.csv'))
stopifnot(nrow(new)==length(datasets)*15L*length(sizes),
          !anyDuplicated(new[c('dataset','repeat_idx','fold_idx','method','k')]),
          all(is.finite(new$AUC)),all(new$n_panel==new$k))
base <- list()
for(ds in datasets) {
  dir <- if(ds %in% c('imvigor210','sosall'))
    paste0('redesign/results_corrected/full_recipe/',ds)
  else paste0('redesign/results_corrected/validation_benchmark/',ds)
  temp <- file.path(scratch,'baseline_current.csv')
  extract(paste0(dir,'/eval_deterministic.csv'),temp)
  d <- read.csv(temp)
  d <- d[d$arm %in% c('GS_full_ungrouped','Boruta','DGE','ElasticNet',
                        'LASSO','RF_importance','mRMR'),]
  d$dataset <- ds
  d$method <- d$arm
  base[[ds]] <- d[c('dataset','repeat_idx','fold_idx','method','k','AUC')]
}
all <- rbind(do.call(rbind,base),new[c('dataset','repeat_idx','fold_idx','method','k','AUC')])
all <- all[all$k %in% sizes,]
stopifnot(!anyDuplicated(all[c('dataset','repeat_idx','fold_idx','method','k')]))
by_dataset_size <- aggregate(AUC~dataset+method+k,all,mean)
by_dataset_primary <- aggregate(AUC~dataset+method,all[all$k %in% c(10,20,50),],mean)
by_validation_size <- aggregate(AUC~method+k,
  by_dataset_size[by_dataset_size$dataset!='sosall',],mean)
by_validation_primary <- aggregate(AUC~method,
  by_dataset_primary[by_dataset_primary$dataset!='sosall',],mean)
tables <- list(by_dataset_size=by_dataset_size,
               by_dataset_primary=by_dataset_primary,
               by_validation_size=by_validation_size,
               by_validation_primary=by_validation_primary)
for(name in names(tables))
  write.csv(tables[[name]],file.path(out,paste0(name,'.csv')),row.names=FALSE)
# Fixed-k Nogueira stability on the candidate-pool intersection per dataset.
rows <- list()
for(ds in datasets) {
  split_dir <- if(ds %in% c('imvigor210','sosall'))
    paste0('redesign/results_corrected/grouped_benchmark/',ds)
  else paste0('redesign/results_corrected/validation_benchmark/',ds)
  rank_dir <- if(ds %in% c('imvigor210','sosall'))
    paste0('redesign/results_corrected/full_recipe/',ds)
  else split_dir
  pools <- vector('list',15)
  ranking <- lapply(methods,function(m) vector('list',15))
  names(ranking) <- methods
  for(i in 1:15) {
    r <- (i-1L)%/%5L+1L; f <- (i-1L)%%5L+1L
    temp <- file.path(scratch,'stability_current.rds')
    extract(sprintf('%s/split_r%d_f%d.rds',split_dir,r,f),temp)
    sp <- readRDS(temp)
    pools[[i]] <- sp$pools$var2000
    temp <- file.path(scratch,'stability_gs_current.csv')
    gs_name <- if(ds %in% c('imvigor210','sosall')) 'full_ungrouped'
               else 'GS_full_ungrouped'
    extract(sprintf('%s/ranking_r%d_f%d_%s.csv',rank_dir,r,f,gs_name),temp)
    ranking[['GS_full_ungrouped']][[i]] <- read.csv(temp)$gene
    for(m in 'Stabl') {
      path <- file.path(out,sprintf('%s_r%d_f%d_%s.csv',ds,r,f,m))
      ranking[[m]][[i]] <- read.csv(path)$gene
    }
  }
  common <- Reduce(intersect,pools)
  p <- length(common); stopifnot(p>500)
  for(m in methods) for(k in sizes) {
    Z <- do.call(rbind,lapply(ranking[[m]],function(rank)
      as.integer(common %in% head(rank,k))))
    pj <- colMeans(Z); kbar <- sum(pj)
    variation <- sum((15/14)*pj*(1-pj))
    value <- 1-variation/(kbar*(1-kbar/p))
    rows[[length(rows)+1L]] <- data.frame(dataset=ds,method=m,k=k,
      Nogueira=value,common_genes=p,mean_projected_size=kbar)
  }
  cat(ds,'common candidates',p,'\n')
}
stability <- do.call(rbind,rows)
write.csv(stability,file.path(out,'stability_by_dataset_size.csv'),row.names=FALSE)
write.csv(aggregate(Nogueira~dataset+method,
  stability[stability$k %in% c(10,20,50),],mean),
  file.path(out,'stability_by_dataset_primary.csv'),row.names=FALSE)
cat('All comparison summaries complete\n')
