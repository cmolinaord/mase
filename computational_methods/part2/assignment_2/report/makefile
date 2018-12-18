OUTPUT_DIR = pdf
SRC_DIR = tex
TEX = latexmk -cd- -pdf -jobname=$(OUTPUT_DIR)/$(NAME) --shell-escape
NAME = p02_wing

.PHONY: all $(NAME) $(NAME)_prev clean clean_all count_words debug_label_chapters debug_warnings
.SILENT: all $(NAME)_prev count_words debug_label_chapters debug_warnings debug_figures debug_figures_unused

all:
	echo "make $(NAME):                      Compile the code from latex to pdf into pdf/ dir"
	echo "make $(NAME)_prev:                 Show the pdf result"
	echo "make clean:                        Remove the auxiliary files and logs in the output dir"
	echo "make clean_all:                    Remove all files in output dir, incluiding .pdf"
	echo "count_words:                       Count number of words in the whole document"
	echo "debug_warnings:                    Review '%WARNINGS' messages during the text put by the writer"
	echo "debug_figures:                     Show all the figures calls or references"
	echo "debug_figures_unused:              Show all the images in the repo that are not actually used in the document"
	echo "debug_label_chapters:              Show all the chapters that do not have a \label"

# Main compilation
##############################
$(NAME): $(OUTPUT_DIR)/$(NAME).pdf

$(OUTPUT_DIR)/$(NAME).pdf: $(SRC_DIR)/main.tex
	mkdir -p $(OUTPUT_DIR)
	$(TEX) $(SRC_DIR)/main.tex

$(NAME)_prev: $(NAME)
	$(TEX) $(SRC_DIR)/main.tex -pv

# Miscelaneous
##############################

clean:
	rm $(OUTPUT_DIR)/*{aux,log,toc,out,fls,lof,fdb*} 2> /dev/null

clean_all:
	rm $(OUTPUT_DIR)/*

count_words:
	pdftotext "$(OUTPUT_DIR)/$(FILENAME).pdf" - |grep -v "^[0-9]" | wc -w

debug_label_chapters:
	echo "Mostrando capitulos que no tienen ningun /label"
	grep -n -E '(\\chapter|\\(sub)*section)' $(SRC_DIR)/*.tex | grep -v -E '\\label' | sed 's/:/:\t/g'

debug_warnings:
	echo "Mostrando lineas con WARNINGS:"
	grep -n -A1 "^%WARNING" $(SRC_DIR)/*.tex | grep -v "WARNING" | sed 's/-%/: /'
	echo " "

debug_figures:
	echo "Mostrando todas las llamadas a etiquetas o referencias a figuras"
	grep -n "fig:" $(SRC_DIR)/*.tex

debug_figures_unused:
	./scripts/search_unused.sh
