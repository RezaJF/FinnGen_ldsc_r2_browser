require("rmarkdown")

# build
setwd("../site")
rmarkdown::render_site()

render('correlation_plots.Rmd', output_dir='../docs')
render('correlation_plots_agglo.Rmd', output_dir='../docs')

render('index.Rmd', output_dir='../docs')

# yaml for per-pheno sites
load("../Rdata_outputs/geno_correlation_sig.Rdata")
pheno <- unique(geno_corr_df$p1)

for(i in pheno) {
	render('../site/rg_phenotype_template.Rmd', params=list(pheno=i, dat='../Rdata_outputs/geno_correlation_sig.Rdata'),
		output_file=paste0('rg_summary_', i, '.html'), output_dir='../docs')
}

## Removed sex-stratified renders
