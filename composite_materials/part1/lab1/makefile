OUTPUT_DIR = pdf
SRC_DIR = tex
TEX = latexmk -cd- -pdf -jobname=$(OUTPUT_DIR)/$(SHORT_NAME) --shell-escape
FILENAME = Lab_1_report
SHORT_NAME = lab1

.PHONY: all $(SHORT_NAME) $(SHORT_NAME)_prev clean clean_all count_words debug_label_chapters debug_warnings
.SILENT: all $(SHORT_NAME)_prev count_words debug_label_chapters debug_warnings debug_figures debug_figures_unused
.SILENT: list_bibliography_available list_bibliography_used update_bibtex

all:
	echo "make $(SHORT_NAME):                Compile the code from latex to pdf into pdf/ dir"
	echo "make $(SHORT_NAME)_prev:           Show the pdf result"
	echo "make clean:                        Remove the auxiliary files and logs in the output dir"
	echo "make clean_all:                    Remove all files in output dir, incluiding .pdf"
	echo "count_words:                       Count number of words in the whole document"
	echo "update_bibtex:                     Update the file of bibliography from mendeley database"
	echo "debug_warnings:                    Review '%WARNINGS' messages during the text put by the writer"
	echo "debug_figures:                     Show all the figures calls or references"
	echo "debug_figures_unused:              Show all the images in the repo that are not actually used in the document"
	echo "debug_label_chapters:              Show all the chapters that do not have a \label"
	echo "list_bibliography_available:       List all the entries in bibliography available in database"
	echo "list_bibliography_used:            List the actual bibliography referenced within the document"

# Main compilation
##############################
$(SHORT_NAME): $(OUTPUT_DIR)/$(SHORT_NAME).pdf

$(OUTPUT_DIR)/$(SHORT_NAME).pdf: $(SRC_DIR)/main.tex
	mkdir -p $(OUTPUT_DIR)
	$(TEX) $(SRC_DIR)/main.tex

$(SHORT_NAME)_prev: $(SHORT_NAME)
	$(TEX) $(SRC_DIR)/main.tex -pv

# Miscelaneous
##############################

clean:
	rm $(OUTPUT_DIR)/*{aux,log,toc,out,fls,lof,fdb*}

clean_all:
	rm $(OUTPUT_DIR)/*

count_words:
	pdftotext "$(OUTPUT_DIR)/$(FILENAME).pdf" - |grep -v "^[0-9]" | wc -w

update_bibtex:
	# Sed command is used for unscape special characters in URL lines
	#   For example: {\_} or {\&} will be replaced by _ and &
	# Also lines starting with "file" are deleted
	echo "Update BibTeX from Mendeley database"
	sed '/url/s/{\\\(.\)}/\1/g;/^file/d' ~/.mendeley/$(SHORT_NAME).bib > doc/$(BIBNAME).bib

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

list_bibliography_available:
	cat doc/$(BIBNAME).bib |grep ^@ |cut -d{ -f2 |cut -d, -f1

list_bibliography_used:
	grep -n "\cite" tex/*.tex
