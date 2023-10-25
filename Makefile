## Makefile by
## Christian Panse <cp@fgcz.ethz.ch> 
##
##
## 2023-10-23 
## 
## usage:
##   # for timsTOF 
##   make stageinput && make unzip && make diann && stageoutput
##
##   # for Orbitrap
##   make stageinput && make convert && make diann && stageoutput

SHELL=/bin/bash

.SUFFIXES: .mzML .raw .csv .fp-manifest .fasta .fas .zip
.PHONY : dep clean .clean.raw .clean.zip .clean.fasta all fasta zip __dia __dda test check stagelog zipfasta .check-diann .check-R .check-FASTA diann qc qc_result/proteinAbundances.html diannopts

RAW               = $(shell find . -type f -name "*.raw")
MZML              = $(RAW:.raw=.mzML)

MSCONVERTOPTS     = chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert --mzML --64 --zlib --filter "peakPicking true 1-"


## DIANN specific
DIANN       = /usr/diann/1.8.2_beta_8/linux/diann-1.8.1.8

#FASTA       = $(shell echo $(DIANNFASTA0) | awk '{print $$NF}')

DIANNDEFAULTS = --threads 64 --verbose 1 --qvalue 0.01 --matrices --predictor --met-excision --cut K*,R* --min-pep-len 6 --max-pep-len 30 --unimod4 --smart-profiling
DIANNGENSLIB = --gen-spec-lib

#DIANNINPUT  = $(shell find $$PWD -maxdepth 1 -type d -name "*.d" -o -name "*.mzML" | awk '{print "--f "$$1}' )
#DIANNTEMP   = $(shell pwd)/temp-$(shell date -I)/
#DIANNOUTPUT = $(shell pwd)/out-$(shell date -I)/

DIANNLIBS := --out-lib ${DIANNOUTPUT}/WU${WORKUNITID}_report-lib.tsv --out-lib-copy

list:
	ls -l $(RAW)
	echo $(MZML)

help:
	grep "^## " Makefile

convert: $(MZML)

.raw.mzML:
	#declare -F write_workunitlog && write_workunitlog $(EXTERNALJOBID) "converting to mzML $< ..."
	docker run -t --network none -w ${PWD} -v ${PWD}:${PWD} $(MSCONVERTOPTS) $< -o $@

unzip: 
	ls | grep .d.zip$$ | parallel -j 16 unzip 

diannopts:
	@echo DIANNOPTS = $(DIANNTMP)
	@echo "--------"
	@$(call echo, $(EXTERNALJOBID), $(DIANNTMP))

diann: diannopts
	mkdir -p ${DIANNTEMP}
	mkdir -p ${DIANNOUTPUT}
	nice -19 $(DIANN) $(DIANNCFG) $(DIANNINPUT) $(DIANNLIBS)\
    --temp ${DIANNTEMP} \
    --out ${DIANNOUTPUT}/WU${WORKUNITID}_report.tsv \
    $(DIANNFASTA) | tee diann.log.txt 


.clean-mzML:
	find . -type f -name "*.mzML" -exec $(RM) -v {} \;

clean: .clean-mzML
