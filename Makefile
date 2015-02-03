.PHONY: default all

default: all

all: legacy-850um-paper.pdf

legacy-850um-paper.pdf: legacy-850um-paper.tex legacy-850um-paper.bib
	pdflatex legacy-850um-paper
	bibtex legacy-850um-paper
	pdflatex legacy-850um-paper
	pdflatex legacy-850um-paper
