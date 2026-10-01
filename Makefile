.PHONY: all dependencies simulations comparator calibrated-power variance confounding curves application validate assets pdf smoke audit compact

all:
	Rscript run_all.R

dependencies:
	Rscript R/check_dependencies.R

compact:
	Rscript R/compress_results.R

simulations:
	Rscript R/run_simulations.R

comparator:
	Rscript R/run_comparator.R
	Rscript R/test_comparator_outputs.R

calibrated-power:
	Rscript R/test_calibrated_power.R
	Rscript R/run_calibrated_power.R
	Rscript R/summarize_calibrated_power.R
	Rscript R/test_calibrated_power_outputs.R
	Rscript R/make_calibrated_power_assets.R

variance:
	Rscript R/run_variance_sensitivity.R
	Rscript R/test_variance_sensitivity.R

application:
	Rscript R/run_application.R
	Rscript R/run_confounding_application.R

confounding:
	Rscript R/run_confounding_simulation.R
	Rscript R/summarize_confounding.R
	Rscript R/run_confounding_application.R
	Rscript R/test_confounding_sensitivity.R
	Rscript R/test_confounding_outputs.R
	Rscript R/make_confounding_assets.R

validate:
	Rscript R/test_results_io.R
	Rscript R/check_dependencies.R
	Rscript R/test_effect_map.R
	Rscript R/test_application.R
	Rscript R/test_calibrated_power.R
	Rscript R/test_calibrated_power_outputs.R
	Rscript R/validate_outputs.R
	Rscript R/test_comparator_outputs.R
	Rscript R/test_variance_sensitivity.R
	Rscript R/test_confounding_sensitivity.R
	Rscript R/test_confounding_outputs.R
	Rscript R/test_confounding_curve.R
	Rscript R/test_confounding_curve_outputs.R

curves:
	Rscript R/run_confounding_curve.R
	Rscript R/summarize_confounding_curve.R
	Rscript R/test_confounding_curve.R
	Rscript R/test_confounding_curve_outputs.R
	Rscript R/make_confounding_curve_assets.R

assets:
	Rscript R/summarize_calibrated_power.R
	Rscript R/make_calibrated_power_assets.R
	Rscript R/paired_fit_audit.R
	JKSS_REUSE_RESULTS=1 Rscript R/make_manuscript_assets.R
	JKSS_VARIANCE_ASSETS_ONLY=1 Rscript R/run_variance_sensitivity.R
	Rscript R/summarize_confounding.R
	Rscript R/make_confounding_assets.R
	Rscript R/summarize_confounding_curve.R
	Rscript R/make_confounding_curve_assets.R

pdf:
	JKSS_REUSE_RESULTS=1 Rscript run_all.R

smoke:
	Rscript R/test_results_io.R
	Rscript R/test_effect_map.R
	Rscript R/smoke_test.R

audit:
	Rscript R/test_calibrated_power_outputs.R
	Rscript R/validate_outputs.R
	Rscript R/test_comparator_outputs.R
	Rscript R/paired_fit_audit.R
