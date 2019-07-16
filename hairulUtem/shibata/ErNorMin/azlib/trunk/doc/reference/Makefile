################################################################################
# Makefile for LaTeX
#
# $Id: Makefile 1031 2016-01-04 05:37:08Z hayashi $
#
################################################################################

## PROJECT NAME
PRONAME = azlib-manual

## SOURCE
REPSRC  = main.tex # BIBTEX UNUSED
#THESRC  = main.tex # BIBTEX USED

## for Linux
TEX     = platex
BIBTEX  = pbibtex
DVIPDF  = dvipdfmx
RM      = rm -f
OPEN    = evince

## for backup
TGZNAME = backup-$(PRONAME)_`date +%Y-%m-%d`

################################################################################

SRC     = $(REPSRC) $(THESRC)
ALL     = $(SRC:.tex=)
DUM     = $(SRC:.tex=.dum)
AUX     = $(SRC:.tex=.aux)
LOG     = $(SRC:.tex=.log)
BBL     = $(THESRC:.tex=.bbl)
DVI     = $(SRC:.tex=.dvi)
PDF     = $(SRC:.tex=.pdf)

$(ALL) : $(SRC) $(INC) $(AUX) $(BBL) $(LOG) $(DVI) $(DUM)
	@ echo =======================================================
	@ echo  END : $(SRC)
	@ echo =======================================================


## SUFFIX RULES

.SUFFIXES : .dum .aux .log .bbl .dvi .pdf

.tex.aux :
	@ echo ==================== Make AUX [$@] ====================
	$(TEX) $< ; $(RM) $*.log

.tex.bbl :
	@ echo ==================== Make BBL [$@] ====================
	$(BIBTEX) $*

.tex.log :
	@ echo ==================== Make LOG [$@] ====================
	$(TEX) $<

.tex.dvi : 
	@ echo ==================== Make DVI [$@] ====================
	$(TEX) $<

.tex.dum : 
	@ echo ==================== Make DUM [$@] ====================
	$(TEX) $<

.dvi.pdf :
	@ echo ==================== Make PDF [$@] ====================
	$(DVIPDF) $<


## COMMANDS

all : $(ALL)

edit :
	@ $(EDITOR) $(SRC)

pdf : $(ALL) $(PDF)

show :
	@ $(OPEN) $(PDF) &

backup :
	@ tar zcvf $(TGZNAME).tgz *
	@ echo "$(TGZNAME).tgz created."

clean :
	@ $(RM) $(DVI) $(AUX) $(DUM) $(LOG) $(BBL) \
	*.blg *.lof *.lot *.toc *.out *.nav *.snm

## EOF
