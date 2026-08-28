.PHONY: all simulations application validate assets pdf smoke

all:
	Rscript run_all.R

simulations:
	Rscript R/run_simulations.R

application:
	Rscript R/run_application.R

validate:
	Rscript R/validate_outputs.R

assets:
	Rscript R/make_manuscript_assets.R

pdf:
	CS_SKIP_SIMULATIONS=1 Rscript run_all.R

smoke:
	Rscript R/smoke_test.R
