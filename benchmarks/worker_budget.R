# ==============================================================================
#  worker_budget.R -- a worker cap shared across concurrently running R processes
# ==============================================================================
#
#  WHY THIS FILE EXISTS
#
#  Each benchmark script caps itself at MAX_PARALLEL_WORKERS (7). That cap is
#  per-process, so it does nothing when more than one script runs at a time:
#  three concurrent runs took 7 workers each, 21 in total, filled 24 GB of RAM,
#  drove swap to 15 GB and crashed the machine. Nothing in any script could
#  detect this, because no script could see the others.
#
#  The registry below makes the cap global. Every run records how many workers
#  it holds in a small file named after its own PID. A starting run sums the
#  live claims and takes only what is left. Entries whose PID is no longer a
#  running R process are deleted, so a killed run or a machine crash cannot
#  leave the budget permanently consumed.
#
#  Usage, immediately before a cluster is created:
#
#      source("benchmarks/worker_budget.R")
#      n_parallel_jobs <- gs_claim_workers(n_parallel_jobs, label = "GSE20194")
#
#  gs_claim_workers() returns the number granted, which is <= what was asked
#  for. It stops with an explanatory message when nothing is available, rather
#  than starting anyway, because starting anyway is what crashed the machine.
# ==============================================================================

# Total workers permitted across ALL concurrent runs on this machine. Raise it
# only after measuring per-worker RSS with `ps`; the previous crash happened at
# 21 workers on 24 GB.
GS_GLOBAL_WORKER_CAP <- as.integer(Sys.getenv("GS_GLOBAL_WORKER_CAP", "7"))

# Kept outside the project tree so that a `git clean` or a fresh clone cannot
# wipe the registry while runs are using it.
GS_WORKER_REGISTRY <- Sys.getenv(
  "GS_WORKER_REGISTRY",
  file.path(path.expand("~"), ".geneselectr_workers")
)

# --- helpers ------------------------------------------------------------------

# TRUE only when this PID is currently a running R process. The command check
# matters because PIDs are recycled: without it, an unrelated process inheriting
# a dead run's PID would keep that run's workers reserved forever.
.gs_pid_is_live_r <- function(pid) {
  out <- suppressWarnings(system2("ps", c("-p", pid, "-o", "command="),
                                  stdout = TRUE, stderr = FALSE))
  length(out) > 0 && any(grepl("exec/R|Rscript", out))
}

# dir.create() is atomic on POSIX, so it doubles as a mutex. Without it two
# scripts starting together could both read "0 in use" and both claim 7.
.gs_with_lock <- function(expr, timeout_s = 30) {
  lock <- file.path(GS_WORKER_REGISTRY, ".lock")
  deadline <- Sys.time() + timeout_s
  repeat {
    if (suppressWarnings(dir.create(lock, showWarnings = FALSE))) break
    # A lock older than the timeout belongs to a process that died holding it.
    age <- difftime(Sys.time(), file.info(lock)$mtime, units = "secs")
    if (!is.na(age) && age > timeout_s) { unlink(lock, recursive = TRUE); next }
    if (Sys.time() > deadline) {
      warning("worker registry lock is stuck; proceeding without it")
      break
    }
    Sys.sleep(0.2)
  }
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  force(expr)
}

# Drop registry entries for runs that are no longer alive, and return the
# surviving ones as a data frame.
.gs_reap <- function() {
  files <- list.files(GS_WORKER_REGISTRY, pattern = "^claim_[0-9]+$",
                      full.names = TRUE)
  if (length(files) == 0)
    return(data.frame(pid = integer(0), workers = integer(0),
                      label = character(0), started = character(0),
                      stringsAsFactors = FALSE))
  rows <- lapply(files, function(f) {
    rec <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(rec) || !.gs_pid_is_live_r(rec$pid)) { unlink(f); return(NULL) }
    data.frame(pid = rec$pid, workers = rec$workers, label = rec$label,
               started = rec$started, stringsAsFactors = FALSE)
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0)
    return(data.frame(pid = integer(0), workers = integer(0),
                      label = character(0), started = character(0),
                      stringsAsFactors = FALSE))
  do.call(rbind, rows)
}

# --- public interface ---------------------------------------------------------

#' Report what is currently holding workers.
gs_worker_status <- function() {
  dir.create(GS_WORKER_REGISTRY, showWarnings = FALSE, recursive = TRUE)
  held <- .gs_with_lock(.gs_reap())
  used <- if (nrow(held) == 0) 0L else sum(held$workers)
  cat(sprintf("Worker budget: %d of %d in use, %d free\n",
              used, GS_GLOBAL_WORKER_CAP, GS_GLOBAL_WORKER_CAP - used))
  if (nrow(held) > 0) {
    for (i in seq_len(nrow(held)))
      cat(sprintf("  pid %-7d %2d worker(s)  %-14s since %s\n",
                  held$pid[i], held$workers[i], held$label[i], held$started[i]))
  }
  invisible(held)
}

#' Claim up to `requested` workers from the machine-wide budget.
#'
#' Returns the number granted. Stops when the budget is exhausted, listing the
#' runs that hold it, so the caller can decide what to stop rather than
#' oversubscribing RAM.
gs_claim_workers <- function(requested, label = "run", min_workers = 1L) {
  requested <- as.integer(requested)
  dir.create(GS_WORKER_REGISTRY, showWarnings = FALSE, recursive = TRUE)

  granted <- .gs_with_lock({
    held <- .gs_reap()
    used <- if (nrow(held) == 0) 0L else sum(held$workers)
    free <- GS_GLOBAL_WORKER_CAP - used

    if (free < min_workers) {
      msg <- sprintf(paste0(
        "\nCannot start: the machine-wide worker budget is full.\n",
        "  cap %d, in use %d, free %d\n%s\n",
        "Stop one of the runs above, or wait for it to finish.\n",
        "Raise the cap only after measuring RAM: GS_GLOBAL_WORKER_CAP=N\n"),
        GS_GLOBAL_WORKER_CAP, used, free,
        paste(sprintf("  pid %-7d %2d worker(s)  %s (since %s)",
                      held$pid, held$workers, held$label, held$started),
              collapse = "\n"))
      stop(msg, call. = FALSE)
    }

    g <- min(requested, free)
    saveRDS(list(pid = Sys.getpid(), workers = g, label = label,
                 started = format(Sys.time(), "%H:%M")),
            file.path(GS_WORKER_REGISTRY, sprintf("claim_%d", Sys.getpid())))
    g
  })

  if (granted < requested)
    cat(sprintf(paste0(
      "  NOTE: asked for %d workers, granted %d -- other runs hold the rest.\n",
      "        (`Rscript -e 'source(\"benchmarks/worker_budget.R\"); ",
      "gs_worker_status()'` to see them)\n"), requested, granted))

  # Release on normal exit. A hard kill or a crash leaves the file behind, which
  # is why .gs_reap() checks PID liveness rather than trusting cleanup.
  reg.finalizer(globalenv(), function(e) gs_release_workers(), onexit = TRUE)
  granted
}

#' Release this process's claim.
gs_release_workers <- function() {
  f <- file.path(GS_WORKER_REGISTRY, sprintf("claim_%d", Sys.getpid()))
  if (file.exists(f)) unlink(f)
  invisible(NULL)
}
