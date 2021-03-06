.PHONY: clean data lint dependencies

#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
PROJECT_NAME = tls1.q
R_INTERPRETER = /usr/bin/env Rscript

#################################################################################
# COMMANDS                                                                      #
#################################################################################

## Installs R Dependencies
dependencies: test_environment
	$(R_INTERPRETER) -e 'renv::restore()'
	$(R_INTERPRETER) -e 'if (tinytex::tinytex_root() == "") tinytex::install_prebuilt()'

## Prepares dataset from individual experiments
0 data: dependencies
	@echo ">>> Preparing dataset..."
	$(R_INTERPRETER) R/osemn/obtain.R
	$(R_INTERPRETER) R/osemn/scrub.R

## Makes Exploratory Data Analysis and outputs correlogram, pairs plot and EDA notebooks
1 eda: data
	@echo ">>> Exploring data for models..."
	$(R_INTERPRETER) R/osemn/explore.R
	#@echo ">>> Exploring data by individual dataset..."
	#$(R_INTERPRETER) R/eda.R
	@echo ">>> Done! See results at $(PROJECT_DIR)/data/output and $(PROJECT_DIR)/vignettes"

## Reproduces research modelling survival and regression
2 models: data
	@echo ">>> Modeling data and summarising results..."
	$(R_INTERPRETER) R/osemn/model.R
	@echo ">>> Done! See results at $(PROJECT_DIR)/data/output"

## Deletes all compiled R ouputs, Rmarkdown and (La)TeX files
clean:
	find ./data/output -type f -name "*.txt" -delete
	find ./data/output -type f -name "*.png" -delete
	find ./data/output -type f -name "*.pdf" -delete
	find ./vignettes -type d -name "figure" -exec rm -rv {} +
	find ./vignettes -type d -name "img" -exec rm -rv {} +
	find ./vignettes -type f -name "*.pdf" -delete
	find ./vignettes -type f -name "*.Rnw" -delete
	find ./vignettes -type f -name "*.log" -delete

## Lints the project using lintr
lint:
	$(R_INTERPRETER) -e 'lintr::lint_package()'

## Sets up R interpreter environment
setup:
	@echo ">>> Installing renv if not already installed."
	$(R_INTERPRETER) setup.R
	@echo ">>> New renv created! It automatically activates running R in $(PROJECT_DIR)"

## Test if R environment is setup correctly
test_environment:
	$(R_INTERPRETER) test_environment.R

#################################################################################
# PROJECT RULES                                                                 #
#################################################################################



#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

# Inspired by <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>
.PHONY: help
help:
	@echo "$$(tput bold)Available commands:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')