
# Use bash
SHELL := /bin/bash

# Default director paths
MAIN = ./main
CORE = ./source
BIN  = ./bin

# Run
run : all
	@(( bin/mkped | tee >( ./chop_ped >&3 ) | bin/recgen ) 3>&1 )\
	| bin/treediff

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

# Recipe to compile tree_diff__main
$(BIN)/treediff : $(MAIN)/tree_diff_main.cpp $(CORE)/*
	@echo "Compiling treediff into $(BIN)"
	@mkdir -p $(BIN)
	@g++ -std=c++17 $(CORE)/*.cpp $(MAIN)/tree_diff_main.cpp \
	-o $(BIN)/treediff
