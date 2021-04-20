
# Use bash
SHELL := /bin/bash

# Default directory paths
MAIN = ./main
CORE = ./source
LOGS = ./logs
BIN  = ./bin

# Args to programs
GCC_ARGS ?= -w -std=c++17
PED_ARGS ?=
REC_ARGS ?=
DIF_ARGS ?=

# Default: compile all and run interactively
interactive : all
	@./run_interactively

# Run
run : all
	@(( bin/mkped $(PED_ARGS) | tee >( ./chop_ped >&3 ) | \
	bin/recgen $(REC_ARGS) ) 3>&1 ) | \
	bin/treediff $(DIF_ARGS)

# Debug recipe
debug : debug_flag all
debug_flag :
	$(eval GCC_ARGS = $(GCC_ARGS) -g -pg)

# Compile all
all : $(BIN)/mkped $(BIN)/recgen $(BIN)/treediff $(BIN)/treeinfo
	@mkdir -p $(LOGS)

# Recipe for compiling main files into bin
$(BIN)/% : $(MAIN)/%_main.cpp $(CORE)/*
	@echo "Compiling $(@F) into $(BIN)"
	@mkdir -p $(BIN)
	@g++ $(GCC_ARGS) $(CORE)/*.cpp $(MAIN)/$(@F)_main.cpp \
	-Ofast -o $@
