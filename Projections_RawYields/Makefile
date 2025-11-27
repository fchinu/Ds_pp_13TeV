###################################
# Makefile for projections and fitting
###################################

# Configuration variables
NAME ?= OO
ISDATA ?= true
BASE_CONFIG_DIR = ./configs

# Set base directory based on data type
ifeq ($(ISDATA),true)
    BASE_DIR = ./data
else
    BASE_DIR = ./MC
endif

TRIAL_DIR = $(BASE_DIR)/$(NAME)

# Output directories
OUTPUT_PROJ_BASE = $(TRIAL_DIR)/Projections
OUTPUT_PROJ_TEMPLATE = $(OUTPUT_PROJ_BASE)/wo_bdt
OUTPUT_FIT = $(TRIAL_DIR)/RawYields

# Configuration files (adjust names based on data type)
ifeq ($(ISDATA),true)
    CONFIG_PROJ = config_projection_data.yaml
    CONFIG_FIT = config_fit_data.yml
	OUTPUT_PROJ = $(OUTPUT_PROJ_BASE)
else
    CONFIG_PROJ = config_projection_mc.yaml
    CONFIG_PROJ_TEMPLATE = config_projection_mc_template.yaml
    CONFIG_FIT = config_fit_mc.yml
	OUTPUT_PROJ = $(OUTPUT_PROJ_BASE)/w_bdt
endif
CONFIG_CUTSET = cutset.yml
CONFIG_CUTSET_TEMPLATE = cutset_template.yml

# Source configuration files
SRC_CONFIG_PROJ = $(BASE_CONFIG_DIR)/$(CONFIG_PROJ)
SRC_CONFIG_FIT = $(BASE_CONFIG_DIR)/$(CONFIG_FIT)
SRC_CONFIG_CUTSET = $(BASE_CONFIG_DIR)/$(CONFIG_CUTSET)
SRC_CONFIG_CUTSET_TEMPLATE = $(BASE_CONFIG_DIR)/$(CONFIG_CUTSET_TEMPLATE)

# Copied configuration files
DEST_CONFIG_PROJ = $(OUTPUT_PROJ)/$(CONFIG_PROJ)
DEST_CONFIG_PROJ_TEMPLATE = $(OUTPUT_PROJ_TEMPLATE)/$(CONFIG_PROJ_TEMPLATE)
DEST_CONFIG_FIT = $(OUTPUT_FIT)/$(CONFIG_FIT)
DEST_CONFIG_CUTSET = $(OUTPUT_PROJ)/$(CONFIG_CUTSET)
DEST_CONFIG_CUTSET_TEMPLATE = $(OUTPUT_PROJ_TEMPLATE)/$(CONFIG_CUTSET_TEMPLATE)

# Scripts (adjust paths as needed)
PROJ_SCRIPT = project_data_from_sparse.py
FIT_SCRIPT = get_raw_yields.py
TEMPLATE_SCRIPT = get_templates.py

# Output files (adjust based on what your scripts actually produce)
PROJ_OUTPUT = $(OUTPUT_PROJ)/projections.done
FIT_OUTPUT = $(OUTPUT_FIT)/rawYields.done
TEMPLATE_OUTPUT = $(OUTPUT_PROJ_TEMPLATE)/dplus_corr_bkg_template.root
SOURCE_FILE = $(TRIAL_DIR)/source.sh

# Define the targets
.PHONY: all clean-output help setup

# Default target
ifeq ($(ISDATA),true)
all: all-data
else
all: all-mc
endif

# Data
all-data: setup project fit

# MC
all-mc: setup project template

# Help target
help:
	@echo "Available targets:"
	@echo "  all        	- Run complete pipeline (setup -> project -> fit)"
	@echo "  setup      	- Create directories and copy config files"
	@echo "  project    	- Run projections"
	@echo "  fit        	- Run fitting"
	@echo "  clean-output	- Remove output files but keep configs"
	@echo "  help       	- Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  NAME=$(NAME) - Trial name (override with NAME=your_trial)"
	@echo "  ISDATA=$(ISDATA) - Data type flag (true for data, false for MC)"
	@echo ""
	@echo "Examples:"
	@echo "  make NAME=PP ISDATA=true    - Run data analysis for PP trial"
	@echo "  make NAME=OO ISDATA=false   - Run MC analysis for OO trial"
	@echo ""
	@echo "Note: Config files are automatically modified to use correct output directories"

# Setup directories and copy configs
ifeq ($(ISDATA),true)
setup: $(OUTPUT_PROJ) $(OUTPUT_FIT) $(DEST_CONFIG_PROJ) $(DEST_CONFIG_FIT) $(DEST_CONFIG_CUTSET) $(SOURCE_FILE)
else
setup: $(OUTPUT_PROJ) $(DEST_CONFIG_PROJ) $(DEST_CONFIG_CUTSET) $(OUTPUT_PROJ_TEMPLATE) $(DEST_CONFIG_PROJ_TEMPLATE) $(DEST_CONFIG_CUTSET_TEMPLATE) $(SOURCE_FILE)
endif

# Create directories
$(OUTPUT_PROJ):
	@mkdir -p $@

$(OUTPUT_FIT):
	@mkdir -p $@

$(OUTPUT_PROJ_TEMPLATE):
	@mkdir -p $@

# Create source file with environment variables
$(SOURCE_FILE):
	@echo "#!/bin/bash" > $@
	@echo "export NAME=$(NAME)" >> $@
	@echo "export ISDATA=$(ISDATA)" >> $@

# Copy and modify configuration files
$(DEST_CONFIG_PROJ): $(SRC_CONFIG_PROJ) $(DEST_CONFIG_CUTSET) | $(OUTPUT_PROJ)
	@cp $< $@
	@sed -i 's|directory: .*|directory: "$(abspath $(OUTPUT_PROJ))"|g' $@
	@sed -i 's|cutset: .*|cutset: "$(abspath $(DEST_CONFIG_CUTSET))"|g' $@

$(DEST_CONFIG_PROJ_TEMPLATE): $(DEST_CONFIG_PROJ)
	@cp $< $@
	@sed -i 's|directory: .*|directory: "$(abspath $(OUTPUT_PROJ_TEMPLATE))"|g' $@
	@sed -i 's|cutset: .*|cutset: "$(abspath $(DEST_CONFIG_CUTSET_TEMPLATE))"|g' $@

$(DEST_CONFIG_FIT): $(SRC_CONFIG_FIT) $(DEST_CONFIG_CUTSET) | $(OUTPUT_FIT)
	@cp $< $@
	@sed -i 's|directory: .*|directory: "$(abspath $(OUTPUT_FIT))"|g' $@
	@sed -i 's|cutset: .*|cutset: "$(abspath $(DEST_CONFIG_CUTSET))"|g' $@
	@sed -i 's|data: .*|data: "$(abspath $(OUTPUT_PROJ)/projections.root)"|g' $@

$(DEST_CONFIG_CUTSET): $(SRC_CONFIG_CUTSET) | $(OUTPUT_PROJ)
	@cp $< $@

$(DEST_CONFIG_CUTSET_TEMPLATE): $(SRC_CONFIG_CUTSET_TEMPLATE) | $(OUTPUT_PROJ_TEMPLATE)
	@cp $< $@

# Run projections
ifeq ($(ISDATA),false)
$(PROJ_OUTPUT): $(DEST_CONFIG_PROJ) $(DEST_CONFIG_CUTSET) $(PROJ_SCRIPT) $(DEST_CONFIG_PROJ_TEMPLATE)
	@cd $(OUTPUT_PROJ) && find . -name "*.root" -delete && python3 $(abspath $(PROJ_SCRIPT)) $(CONFIG_PROJ) && hadd -f ./projections.root ./*.root && hadd -f ./dplus.root ./dplus_nonprompt.root ./dplus_prompt.root
	@cd $(OUTPUT_PROJ_TEMPLATE) && find . -name "*.root" -delete && python3 $(abspath $(PROJ_SCRIPT)) $(CONFIG_PROJ_TEMPLATE) && hadd -f ./projections.root ./*.root && hadd -f ./dplus.root ./dplus_nonprompt.root ./dplus_prompt.root
	@touch $@
else
$(PROJ_OUTPUT): $(DEST_CONFIG_PROJ) $(DEST_CONFIG_CUTSET) $(PROJ_SCRIPT)
	@cd $(OUTPUT_PROJ) && find . -name "*.root" -delete && python3 $(abspath $(PROJ_SCRIPT)) $(CONFIG_PROJ) && hadd -f ./projections.root ./*.root
	@touch $@
endif

# Run fitting
$(FIT_OUTPUT): $(PROJ_OUTPUT) $(DEST_CONFIG_FIT) $(FIT_SCRIPT)
	@cd $(OUTPUT_FIT) && python3 $(abspath $(FIT_SCRIPT)) $(CONFIG_FIT)
	@touch $@

$(TEMPLATE_OUTPUT): $(DEST_CONFIG_PROJ_TEMPLATE) $(TEMPLATE_SCRIPT)
	@cd $(OUTPUT_PROJ_TEMPLATE) && python3 $(abspath $(TEMPLATE_SCRIPT)) dplus_bkg.root dplus_corr_bkg_template.root
	@touch $@

# Convenient aliases
source: $(SOURCE_FILE)
project: $(PROJ_OUTPUT)
fit: $(FIT_OUTPUT)
template: $(TEMPLATE_OUTPUT)

# Alternative clean that preserves configs
clean-output:
	find $(OUTPUT_PROJ) \( -name "*.done" -o -name "*.root" -o -name "*.pdf" -o -name "*.png" -o -name "*.parquet" \) -delete
	find $(OUTPUT_FIT) \( -name "*.done" -o -name "*.root" -o -name "*.pdf" -o -name "*.png" -o -name "*.parquet" \) -delete
	# Add specific output files to clean as needed

# Debug target to show variables
debug:
	@echo "NAME: $(NAME)"
	@echo "ISDATA: $(ISDATA)"
	@echo "BASE_DIR: $(BASE_DIR)"
	@echo "TRIAL_DIR: $(TRIAL_DIR)"
	@echo "OUTPUT_PROJ: $(OUTPUT_PROJ)"
	@echo "OUTPUT_FIT: $(OUTPUT_FIT)"
	@echo "CONFIG_PROJ: $(CONFIG_PROJ)"
	@echo "CONFIG_FIT: $(CONFIG_FIT)"
	@echo "SRC_CONFIG_PROJ: $(SRC_CONFIG_PROJ)"
	@echo "DEST_CONFIG_PROJ: $(DEST_CONFIG_PROJ)"