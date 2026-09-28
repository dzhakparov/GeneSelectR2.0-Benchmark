#!/usr/bin/env Rscript
extra_lib <- Sys.getenv('GENESELECTR_R_LIB', '')
if(nzchar(extra_lib)) .libPaths(c(extra_lib, .libPaths()))
Sys.setenv(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
           MKL_NUM_THREADS='1', MPLCONFIGDIR='/private/tmp/mplconfig')
repo <- normalizePath(Sys.getenv('GENESELECTR_BENCHMARK_REPO', getwd()),mustWork=TRUE)
commit <- 'a840a6bf6ab36c87ac69ed65a5d506b6f29de21d'
out <- Sys.getenv('GENESELECTR_STABL_OUT',file.path(repo,'redesign/stabl_comparison'))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
scratch <- tempdir()
py <- Sys.getenv('GENESELECTR_STABL_PYTHON','python3')
source(file.path(repo, 'redesign/R/evaluator.R'))
args <- commandArgs(trailingOnly=TRUE)
datasets <- if(length(args)) args[[1]] else c('GSE101794','GSE107994','GSE13355','GSE65682','GSE69683','imvigor210','sosall')
max_splits <- if(length(args)>1) as.integer(args[[2]]) else 15L
extract <- function(path, dest) {
  status <- system2('git', c('-C', shQuote(repo), 'show',
                             shQuote(paste0(commit, ':', path))), stdout=dest)
  if(!identical(status, 0L)) stop('Failed to extract: ', path)
}
for(ds in datasets) {
  cohort <- ds %in% c('sosall','imvigor210')
  dir <- if(cohort) paste0('redesign/results_corrected/grouped_benchmark/',ds)
         else paste0('redesign/results_corrected/validation_benchmark/',ds)
  base_path <- if(cohort) paste0(dir,'/',ds,'_base.rds') else paste0(dir,'/base_data.rds')
  base_file <- file.path(scratch,paste0(ds,'_base.rds'))
  if(!file.exists(base_file)) extract(base_path,base_file)
  base <- readRDS(base_file)
  stopifnot(is.factor(base$outcome))
  for(r in 1:3) for(f in 1:5) {
    if((r-1L)*5L+f > max_splits) next
    key <- sprintf('%s_r%d_f%d',ds,r,f)
    split_file <- file.path(scratch,'current_split.rds')
    extract(sprintf('%s/split_r%d_f%d.rds',dir,r,f),split_file)
    sp <- readRDS(split_file)
    pool <- sp$pools$var2000
    stopifnot(length(pool)==2000L)
    x <- sp$train_raw[,pool,drop=FALSE]
    z <- sp$test_raw[,pool,drop=FALSE]
    mu <- colMeans(x); sdv <- apply(x,2,sd)
    sdv[!is.finite(sdv)|sdv==0] <- 1
    x <- sweep(sweep(x,2,mu,'-'),2,sdv,'/')
    z <- sweep(sweep(z,2,mu,'-'),2,sdv,'/')
    y <- droplevels(base$outcome[sp$train_idx])
    yt <- droplevels(base$outcome[sp$test_idx])
    stopifnot(identical(colnames(x),pool), nrow(x)==length(y),
              nrow(z)==length(yt), all(is.finite(x)), all(is.finite(z)))
    xcsv <- file.path(scratch,'current_X.csv')
    ycsv <- file.path(scratch,'current_y.csv')
    gcsv <- file.path(scratch,'current_groups.csv')
    write.csv(x,xcsv,row.names=FALSE)
    write.csv(data.frame(y=as.integer(y)-1L),ycsv,row.names=FALSE)
    has_groups <- ds=='GSE13355'
    if(has_groups) write.csv(data.frame(group=base$groups[sp$train_idx]),gcsv,row.names=FALSE)
    for(method in 'Stabl') {
      rank_file <- file.path(out,paste0(key,'_',method,'.csv'))
      if(!file.exists(rank_file)) {
        seed <- 520000L + r*1000L + f*100L + 1L
        status <- system2(py,c(shQuote(file.path(out,'rank_stabl.py')),
          '--x',shQuote(xcsv),'--y',shQuote(ycsv),'--groups',
          shQuote(if(has_groups) gcsv else 'none'),'--output',
          shQuote(rank_file),'--seed',seed))
        if(!identical(status,0L)) stop('Ranking failed: ',key,' ',method)
      }
      genes <- read.csv(rank_file,check.names=FALSE)$gene
      stopifnot(length(genes)==2000L, !anyDuplicated(genes), setequal(genes,pool))
      for(k in c(10L,20L,50L,100L,200L,500L)) {
        eval_file <- file.path(out,'evaluation.csv')
        if(file.exists(eval_file)) {
          prior <- read.csv(eval_file)
          if(any(prior$dataset==ds & prior$repeat_idx==r & prior$fold_idx==f &
                 prior$method==method & prior$k==k)) next
        }
        panel <- head(genes,k)
        eval_seed <- 420000L + r*1000L + f*100L +
          match(k,c(10L,20L,50L,100L,200L,500L))
        probs <- predict_with_ensemble(x[,panel,drop=FALSE],y,
                                       z[,panel,drop=FALSE],random_seed=eval_seed)
        row <- data.frame(dataset=ds,repeat_idx=r,fold_idx=f,method=method,
                          k=k,AUC=bench_auc(yt,probs),n_panel=length(panel),
                          evaluation_seed=eval_seed,
                          evaluator_components=paste(attr(probs,'components_used'),collapse='+'))
        write.table(row,eval_file,append=file.exists(eval_file),sep=',',
                    row.names=FALSE,col.names=!file.exists(eval_file))
      }
      cat(sprintf('[%s] %s done\n',key,method))
      flush.console()
    }
    unlink(c(split_file,xcsv,ycsv,gcsv))
  }
}
