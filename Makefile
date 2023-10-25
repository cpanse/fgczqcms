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

.SUFFIXES: .mzML .raw .csv .fasta .fas .zip .xmllint .diann
.PHONY : dep clean .clean.raw .clean.zip test check .check-diann .check-FASTA 

NOW               = $(shell date +"%Y-%m-%dT%H:%M:%S%:z")
RAW               = $(shell find . -type f -name "*.raw")
MZML              = $(RAW:.raw=.mzML)
XMLLINT           = $(MZML:.mzML=.xmllint)
OUTPUT            = $(MZML:.mzML=.diann)

MSCONVERTOPTS     = chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert --mzML --64 --zlib --filter "peakPicking true 1-"

## DIANN specific
DIANN         = /usr/diann/1.8.2_beta_8/linux/diann-1.8.1.8
DIANNFASTA    = --fasta-search --fasta /srv/www/htdocs/FASTA/fgcz_uniprot-proteome_UP000005640.fasta
DIANNDEFAULTS = --threads 64 --verbose 1 --qvalue 0.01 --matrices --predictor --met-excision --cut K*,R* --min-pep-len 6 --max-pep-len 30 --unimod4 --smart-profiling
DIANNTEMP     = $(PWD)/temp


list:
	ls -l $(RAW)
	echo $(MZMLCHECK)
	echo $(XMLLINT)
	echo $(OUTPUT)

help:
	grep "^## " Makefile

convert: $(MZML)
diann: $(OUTPUT)

.raw.mzML:
	docker run -t --network none -w ${PWD} -v ${PWD}:${PWD} $(MSCONVERTOPTS) $< -o $(shell dirname $<)

unzip: 
	ls | grep .d.zip$$ | parallel -j 16 unzip 




statistics:
	@sort -t';' -k4 ../input.txt  > ../input.txt.s
	(find . -type f -name "diaqc_report.stats.tsv" | head -n 1| parallel head -n 1 | tr "\t" ";" | sed 's/File.Name;/File.Name;Md5;Time;Size;/'; \
	find . -type f -name "diaqc_report.stats.tsv" -exec tail -n 1 {} \; | sed -e 's/dump\///' -e 's/\.mzML/\.raw/' | tr "\t" ";" | sort | join -t";" -1 4 -2 1 ../input.txt.s -) \
	| sort -k3 -g -t';'

stat: statistics

.mzML.diann:
	mkdir -p $@ $(DIANNTEMP) && nice -19 $(DIANN) $(DIANNDEFAULTS) --f $< $(DIANNFASTA) --temp $< --temp $(DIANNTEMP) --out $@/diaqc_report.tsv || mv $@ $@.$(NOW).broken

.mzML.xmllint:
	xmllint --noout $< > $@ 

.clean-xmlllint:
	$(RM) $(XMLLINT)

.clean-mzML:
	find . -type f -name "*.mzML" -exec $(RM) -v {} \;
	find . -type d -name "*.mzML" -exec rmdir -v {} \;

.clean-diann:
	$(RM) -rfv $(OUTPUT)

check: $(XMLLINT)
clean: .clean-mzML .clean-xmllint .clean-diann
