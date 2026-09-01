suppressPackageStartupMessages(library(data.table))
msg <- function(...) cat(..., '\n', file = stderr())
set.seed(42)
src <- readLines('scripts/55_gse154378_subgroup_quasideconv.R', encoding = 'UTF-8')
end <- grep("^## ---------- 4", src)
eval(parse(text = paste(src[1:(end-1)], collapse = '\n')))
msg('functions defined')
GROUPS <- c("Normal", "PreX", "GDM", "cHTN")
for (g in GROUPS) {
  msg('run_group', g, 'start')
  r <- run_group(dcast_obs, Rlm, g)
  msg('run_group', g, 'OK')
  gc()
}
msg('run_group x4 done')
for (g in c('PreX','GDM','cHTN')) {
  msg('contrast', g, 'start')
  r <- run_contrast(dcast_obs, Rlm, g)
  msg('contrast', g, 'OK nrow', nrow(r))
  gc()
}
msg('ALL DONE')
