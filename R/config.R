# automatically create a bib database for R packages
knitr::write_bib(c(
	.packages(), 'bookdown', 'knitr', 'rmarkdown'
), 'packages.bib')


## cache is only valid with a specific version of R and session info
## cache will be kept for at most a month (re-compute the next month)
knitr::opts_chunk$set(cache.extra = list(
	R.version, sessionInfo(), format(Sys.Date(), '%Y-%m'), knitr::rand_seed
), 
collapse = TRUE,
comment = '#>'
)
