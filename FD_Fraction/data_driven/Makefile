
###################################
# Makefile for FD fraction evaluation
###################################

# Configuration variables
NAME ?= OO
BASE_CONFIG_DIR = ./configs

TRIAL_DIR = $(NAME)

# Output directories
OUTPUT_CONFIG = $(TRIAL_DIR)/configs
OUTPUT_PROJ = $(TRIAL_DIR)/Projections
OUTPUT_PROJ_PROMPT_ENHANCED = $(OUTPUT_PROJ)/prompt_enhanced
OUTPUT_FIT = $(TRIAL_DIR)/RawYields
OUTPUT_EFF = $(TRIAL_DIR)/Efficiencies
OUTPUT_FRAC = $(TRIAL_DIR)/Fractions
OUTPUT_FRAC_DS = $(TRIAL_DIR)/Fractions/Ds
OUTPUT_FRAC_DPLUS = $(TRIAL_DIR)/Fractions/Dplus
OUTPUT_CONFIG_CUTSET = $(OUTPUT_CONFIG)/cutsets
OUTPUT_CONFIG_PROJ = $(OUTPUT_CONFIG)/Projections
OUTPUT_CONFIG_FIT = $(OUTPUT_CONFIG)/RawYields
OUTPUT_CONFIG_EFF = $(OUTPUT_CONFIG)/Efficiencies

# Configuration files (adjust names based on data type)
CONFIG_MAKE_CUTS = config_make_cuts.yml
CONFIG_PROJ_PROMPT_ENHANCED = config_projection_prompt_enhanced.yml
CONFIG_PROJ = config_projection.yml
CONFIG_FIT_PROMPT_ENHANCED = config_fit_prompt_enhanced.yml
CONFIG_FIT = config_fit.yml
CONFIG_EFF = config_efficiency.yml
CONFIG_CUTSET_PROMPT_ENHANCED = cutset_prompt_enhanced.yml
CONFIG_CUTSET = cutset.yml
CONFIG_CUTVAR = config_cutvar.json

# Source configuration files
SRC_CONFIG_MAKE_CUTS = $(BASE_CONFIG_DIR)/$(CONFIG_MAKE_CUTS)
SRC_CONFIG_PROJ_PROMPT_ENHANCED = $(BASE_CONFIG_DIR)/$(CONFIG_PROJ_PROMPT_ENHANCED)
SRC_CONFIG_PROJ = $(BASE_CONFIG_DIR)/$(CONFIG_PROJ)
SRC_CONFIG_FIT_PROMPT_ENHANCED = $(BASE_CONFIG_DIR)/$(CONFIG_FIT_PROMPT_ENHANCED)
SRC_CONFIG_FIT = $(BASE_CONFIG_DIR)/$(CONFIG_FIT)
SRC_CONFIG_EFF = $(BASE_CONFIG_DIR)/$(CONFIG_EFF)
SRC_CONFIG_CUTSET_PROMPT_ENHANCED = $(BASE_CONFIG_DIR)/$(CONFIG_CUTSET_PROMPT_ENHANCED)
SRC_CONFIG_CUTSET = $(BASE_CONFIG_DIR)/$(CONFIG_CUTSET)
SRC_CONFIG_CUTVAR = $(BASE_CONFIG_DIR)/$(CONFIG_CUTVAR)

# Copied configuration files - FIXED: Use unique destination paths
DEST_CONFIG_MAKE_CUTS = $(TRIAL_DIR)/$(CONFIG_MAKE_CUTS)
DEST_CONFIG_PROJ = $(TRIAL_DIR)/$(CONFIG_PROJ)
DEST_CONFIG_PROJ_PROMPT_ENHANCED = $(TRIAL_DIR)/$(CONFIG_PROJ_PROMPT_ENHANCED)
DEST_CONFIG_FIT = $(TRIAL_DIR)/$(CONFIG_FIT)
DEST_CONFIG_FIT_PROMPT_ENHANCED = $(TRIAL_DIR)/$(CONFIG_FIT_PROMPT_ENHANCED)
DEST_CONFIG_EFF = $(TRIAL_DIR)/$(CONFIG_EFF)
DEST_CONFIG_CUTSET = $(TRIAL_DIR)/$(CONFIG_CUTSET)
DEST_CONFIG_CUTSET_PROMPT_ENHANCED = $(TRIAL_DIR)/$(CONFIG_CUTSET_PROMPT_ENHANCED)
DEST_CONFIG_CUTVAR = $(TRIAL_DIR)/$(CONFIG_CUTVAR)

# .gitignore file for trial directory
TRIAL_GITIGNORE = $(TRIAL_DIR)/.gitignore

# Scripts (adjust paths as needed)
CUTS_SCRIPT = MakeCutsFilesForSyst.py
CONFIG_GENERATOR_SCRIPT = generate_configs.py
PROJECTIONS_SCRIPT = RunProjections.sh
RAWYIELDS_SCRIPT = RunRawYields.sh
EFF_SCRIPT = RunEfficiency.sh
FRAC_SCRIPT = RunCutVariation.sh

# Output files (adjust based on what your scripts actually produce)
CONFIG_GENERATED = $(OUTPUT_CONFIG)/configs.done
CUTS_OUTPUT = $(OUTPUT_CONFIG_CUTSET)/cuts.done
PROJ_OUTPUT = $(OUTPUT_PROJ)/projections.done
FIT_OUTPUT = $(OUTPUT_FIT)/rawYields.done
EFF_OUTPUT = $(OUTPUT_EFF)/efficiencies.done
FRAC_OUTPUT = $(OUTPUT_FRAC)/fd_fraction.done
SOURCE_FILE = $(TRIAL_DIR)/source.sh

# Define the targets
.PHONY: all clean-output help setup generate-configs

# Default target
all: setup generate-configs project fit efficiency fractions

# Help target
help:
	@echo "Available targets:"
	@echo "  all           - Run complete pipeline (setup -> generate-configs -> weight -> efficiency)"
	@echo "  setup         - Create directories and copy config files"
	@echo "  generate-configs - Generate config files for each cutset"
	@echo "  weight        - Run weight calculations"
	@echo "  efficiency    - Run efficiency evaluation"
	@echo "  clean-output  - Remove output files but keep configs"
	@echo "  help          - Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  NAME=$(NAME) - Trial name (override with NAME=your_trial)"
	@echo ""
	@echo "Examples:"
	@echo "  make NAME=PP    - Run analysis for PP trial"
	@echo "  make NAME=OO    - Run analysis for OO trial"
	@echo ""
	@echo "Notes:"
	@echo "  - Config files are automatically generated for each cutset"
	@echo "  - Weights config is only populated if 'apply: True' is found in efficiency config"
	@echo "  - If weights are not needed, the weights script is skipped"

# All directories that need to be created
DIRECTORIES = $(OUTPUT_CONFIG) $(OUTPUT_CONFIG_CUTSET) $(OUTPUT_CONFIG_PROJ) $(OUTPUT_CONFIG_FIT) $(OUTPUT_CONFIG_EFF) $(OUTPUT_PROJ) $(OUTPUT_PROJ_PROMPT_ENHANCED) $(OUTPUT_FIT) $(OUTPUT_EFF) $(OUTPUT_FRAC_DS) $(OUTPUT_FRAC_DPLUS) $(TRIAL_DIR)

# All config files that need to be copied
CONFIG_FILES = $(DEST_CONFIG_CUTSET) $(DEST_CONFIG_CUTSET_PROMPT_ENHANCED) $(DEST_CONFIG_PROJ) $(DEST_CONFIG_PROJ_PROMPT_ENHANCED) $(DEST_CONFIG_FIT) $(DEST_CONFIG_FIT_PROMPT_ENHANCED) $(DEST_CONFIG_EFF) $(DEST_CONFIG_MAKE_CUTS) $(TRIAL_GITIGNORE) $(DEST_CONFIG_CUTVAR)

# Setup directories and copy configs
setup: $(DIRECTORIES) $(CONFIG_FILES) generate-configs $(SOURCE_FILE)

# Create all directories with a single rule
$(DIRECTORIES):
	@mkdir -p $@

# Generic rule for copying config files
$(TRIAL_DIR)/%: $(BASE_CONFIG_DIR)/%
	@echo "Copying config: $< -> $@"
	@cp $< $@

# Create .gitignore in trial directory
$(TRIAL_GITIGNORE): $(TRIAL_DIR)
	@echo "Creating .gitignore in trial directory..."
	@echo "*" > $@

# Generate config files for each cutset
generate-configs: $(CUTS_OUTPUT) $(CONFIG_GENERATED)

# Create source file with environment variables
$(SOURCE_FILE):
	@echo "#!/bin/bash" > $@
	@echo "export NAME=$(NAME)" >> $@

# Create cutset files
$(CUTS_OUTPUT): $(DEST_CONFIG_CUTSET) $(CUTS_SCRIPT) $(DEST_CONFIG_MAKE_CUTS)
	@sed -i 's|in_dir:.*|in_dir: "$(abspath $(TRIAL_DIR))/"|g' $(word 3,$^)
	@sed -i 's|cut_file_central:.*|cut_file_central: "$(CONFIG_CUTSET)"|g' $(word 3,$^)
	@sed -i 's|out_dir:.*|out_dir: "$(abspath $(OUTPUT_CONFIG_CUTSET))/"|g' $(word 3,$^)
	@python3 $(abspath $(CUTS_SCRIPT)) -c $(abspath $(DEST_CONFIG_MAKE_CUTS))
	@touch $@

$(CONFIG_GENERATED): $(OUTPUT_CONFIG_CUTSET) $(OUTPUT_CONFIG_PROJ) $(OUTPUT_CONFIG_FIT) $(OUTPUT_CONFIG_EFF) $(DEST_CONFIG_CUTSET) $(DEST_CONFIG_CUTSET_PROMPT_ENHANCED) $(DEST_CONFIG_PROJ) $(DEST_CONFIG_PROJ_PROMPT_ENHANCED) $(DEST_CONFIG_FIT) $(DEST_CONFIG_FIT_PROMPT_ENHANCED) $(DEST_CONFIG_EFF)
	@echo "Generating config files for each cutset..."
	@python $(CONFIG_GENERATOR_SCRIPT) --trial-dir $(TRIAL_DIR)
	@touch $@

project: $(PROJ_OUTPUT)

$(PROJ_OUTPUT): $(CONFIG_GENERATED)
	@echo "Running projections..."
	@bash $(PROJECTIONS_SCRIPT) $(NAME)
	@touch $@

fit: $(FIT_OUTPUT)

$(FIT_OUTPUT): $(PROJECTIONS_DONE) $(CONFIG_GENERATED)
	@echo "Running raw yields extraction..."
	@bash $(RAWYIELDS_SCRIPT) $(NAME)
	@touch $@

efficiency: $(EFF_OUTPUT)

$(EFF_OUTPUT): $(CONFIG_GENERATED)
	@echo "Running efficiency evaluation..."
	@bash $(EFF_SCRIPT) $(NAME)
	@touch $@

fractions: $(FRAC_OUTPUT)

$(FRAC_OUTPUT): $(FIT_OUTPUT) $(EFF_OUTPUT) $(CONFIG_GENERATED)
	@echo "Calculating FD fraction..."
	@bash $(FRAC_SCRIPT) $(NAME)
	@touch $@

# Clean that preserves configs
clean-output:
	find $(OUTPUT_PROJ) \( -name "*.done" -o -name "*.root" -o -name "*.pdf" -o -name "*.png" -o -name "*.parquet" \) -delete
	find $(OUTPUT_FIT) \( -name "*.done" -o -name "*.root" -o -name "*.pdf" -o -name "*.png" -o -name "*.parquet" \) -delete
	# Add specific output files to clean as needed
