
# Use bash
SHELL := /bin/bash

# Default director paths
MAIN = ./main
CORE = ./source
BIN  = ./bin

# Args to programs
PED_ARGS ?= ''
REC_ARGS ?= ''
DIF_ARGS ?= ''

# Default: compile all and run interactively
interactive : all
	@./run_interactively

# Run
run : all
	@(( bin/mkped $(PED_ARGS) | tee >( ./chop_ped >&3 ) | \
	bin/recgen $(REC_ARGS) ) 3>&1 ) | \
	bin/treediff $(DIF_ARGS)

# Compile all
all : $(BIN)/mkped $(BIN)/recgen $(BIN)/treediff $(BIN)/treeinfo

# Recipe for compiling main files into bin
$(BIN)/% : $(MAIN)/%_main.cpp $(CORE)/*
	@echo "Compiling $(@F) into $(BIN)"
	@mkdir -p $(BIN)
	@g++ -w -std=c++17 $(CORE)/*.cpp $(MAIN)/$(@F)_main.cpp \
	-Ofast -o $@
