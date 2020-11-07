
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
all : $(BIN)/mkped $(BIN)/recgen $(BIN)/treediff

# Recipe to compile pedigree_gen_main
$(BIN)/mkped : $(MAIN)/pedigree_gen_main.cpp $(CORE)/*
	@echo "Compiling mkped into $(BIN)"
	@mkdir -p $(BIN)
	@g++ -w -std=c++17 $(CORE)/*.cpp $(MAIN)/pedigree_gen_main.cpp \
	-o $(BIN)/mkped

# Recipe to compile rec_gen_main
$(BIN)/recgen : $(MAIN)/rec_gen_main.cpp $(CORE)/*
	@echo "Compiling recgen into $(BIN)"
	@mkdir -p $(BIN)
	@g++ -std=c++17 $(CORE)/*.cpp $(MAIN)/rec_gen_main.cpp \
	-o $(BIN)/recgen

# Recipe to compile tree_diff_main
$(BIN)/treediff : $(MAIN)/tree_diff_main.cpp $(CORE)/*
	@echo "Compiling treediff into $(BIN)"
	@mkdir -p $(BIN)
	@g++ -std=c++17 $(CORE)/*.cpp $(MAIN)/tree_diff_main.cpp \
	-o $(BIN)/treediff
