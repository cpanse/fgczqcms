# README

The fgczqcms R package contains shell scripts and a shiny dashboard for handling
mass spec qc workflows (dia, dda) for Orbitrap and timsTof devices.

data are preprocessed as follow:

## pre-processing 

### autoQC01 

- compute host: `cpanse@fgcz-c-072:/scratch/FGCZ-MS-QC-DASHBOARD_A331/qc`
running 

```
while true;
do 
  clear; date;
  make -f Makefile.autoQC01 ;
  sleep 60;
done
```

for 

- [x] Thermo Orbitrap (QEXACTIVE_1|...)

### autoQC03dda
- compute host: running @fgcz-r-033

```
# autoQC03dda AUTOQC4L TIMSTOF
37 * * * * (/scratch/cpanse/timstof/bin/fgcz_timsTOF_autoQC4L.bash   2>&1)  >> /home/cpanse/data/runme_autoQC4L_TIMSTOF_1.log 
17 * * * * (/scratch/cpanse/timstof/bin/fgcz_TIMSTOFFLEX_1_autoQC03dda.bash 2>&1) >> /home/cpanse/data/runme_autoQC4L_TIMSTOFFLEX_1.log 

*/29    *       *       *       *     ( nice -19 /home/cpanse/__projects/2019/20190225-comet/runme_autoQC4L_comet.bash ) 2>/tmp/autoQC4L.2.log > /tmp/autoQC4L.1.log
```

TODO(cp): this script requires refactoring: `/home/cpanse/__projects/2019/20190225-comet/runme_autoQC4L_comet.bash`

- [x] Thermo Orbitrap (QEXACTIVE_1|...)
- [x] Bruker timsTOF 

see also:
- https://gitlab.bfabric.org/cpanse/timstof.git

### autoQC03dia
- compute host: `cpanse@fgcz-c-072:/scratch/FGCZ-MS-QC-DASHBOARD_A331/qc`

```
while true;
do 
  clear;
  date;
  make -f Makefile.DIANN stageinput \
    && make -f Makefile.DIANN convert \
    && make -f Makefile.DIANN diann \
    && make -f Makefile.DIANN output ;
    sleep 60;
done
```

- [x] Thermo Orbitrap (QEXACTIVE_1|...)
- [x] Bruker timsTOF  (TIMSTOF_1|TIMSTOFFLEX_1)


## links

- https://shiny-ms.fgcz.uzh.ch/fgczmsqc-dashboard/
- https://fgcz-bfabric.uzh.ch/bfabric/application/show.html?id=331&tab=details

## contact

Email: cp@fgcz.ethz.ch
