#!/bin/bash

#SBATCH --partition=prx
#SBATCH --nodelist=fgcz-r-033
#SBATCH --job-name=autoQC04
#SBATCH --workdir=/scratch/cpanse/comet/
#SBATCH --mem-per-cpu=1G
#SBATCH --dependency=singleton


# Christian Panse <cp@fgcz.ethz.ch>, 2019-01 

set -o pipefail

scriptname=$(basename $0)
lock="/tmp/${scriptname}"
exec 200>$lock
flock -n 200 || { echo "$0 is already running"; exit 1; }

export scratch=/scratch/cpanse/comet

EMAIL=cp@fgcz.ethz.ch

function report(){
  if [ $1 -eq 0 ];
  then
    #  mail ${EMAIL} -s "$0 succ" < /dev/null
      echo "succ"
  else
    mail ${EMAIL} -s "$0 failed" < /dev/null
  fi


  exit $1
}

trap "{ report $? < /dev/null; }" EXIT

function write_comet_para_pasef() {
  cat <<EOF
peptide_mass_tolerance = 50.0
peptide_mass_units = 2

use_A_ions = 0
use_B_ions = 1
use_C_ions = 0
use_X_ions = 0
use_Y_ions = 1
use_Z_ions = 0
use_NL_ions = 1

fragment_bin_tol = 0.02
fragment_bin_offset = 0.0

precursor_charge = 2 6
EOF
}

function write_comet_para_highres_by() {
  cat <<EOF
# comet_version 2018.01 rev. X
# Comet MS/MS search engine parameters file.
# Everything following the '#' symbol is treated as a comment.
database_name = /srv/www/htdocs/FASTA/fgcz_uniprot-proteome_UP000005640.fasta
decoy_search = 1
num_threads = 32
peptide_mass_tolerance = 10.0
peptide_mass_units = 2
mass_type_parent = 1
mass_type_fragment = 1
precursor_tolerance_type = 1
isotope_error = 0
search_enzyme_number = 1
num_enzyme_termini = 2
allowed_missed_cleavage = 2
variable_mod01 = 15.9949 M 0 3 -1 0 0
variable_mod02 = 0.0 X 0 3 -1 0 0
variable_mod03 = 0.0 X 0 3 -1 0 0
variable_mod04 = 0.0 X 0 3 -1 0 0
variable_mod05 = 0.0 X 0 3 -1 0 0
variable_mod06 = 0.0 X 0 3 -1 0 0
variable_mod07 = 0.0 X 0 3 -1 0 0
variable_mod08 = 0.0 X 0 3 -1 0 0
variable_mod09 = 0.0 X 0 3 -1 0 0
max_variable_mods_in_peptide = 5
require_variable_mod = 0
# low res
# fragment_bin_tol = 1.0005
# fragment_bin_offset = 0.4
## For high-res MS/MS data
fragment_bin_tol = 0.02
fragment_bin_offset = 0.0
theoretical_fragment_ions = 1
use_A_ions = 0
use_B_ions = 1
use_C_ions = 0
use_X_ions = 0
use_Y_ions = 1
use_Z_ions = 0
use_NL_ions = 1
output_sqtstream = 0
output_sqtfile = 0
output_txtfile = 1
output_pepxmlfile = 1
output_percolatorfile = 0
output_outfiles = 0
print_expect_score = 1
num_output_lines = 1
show_fragment_ions = 0
sample_enzyme_number = 1
scan_range = 0
precursor_charge = 0
override_charge = 0
ms_level = 2
activation_method = ALL
digest_mass_range = 600.0
num_results = 100
skip_researching = 1
max_fragment_charge = 3
max_precursor_charge = 6
nucleotide_reading_frame = 0
clip_nterm_methionine = 1
spectrum_batch_size = 0
decoy_prefix = REV_
output_suffix = .comet
mass_offsets = 0.0
minimum_peaks = 10
minimum_intensity = 0
remove_precursor_peak = 0
remove_precursor_tolerance = 1.5
clear_mz_range = 0.0
add_Cterm_peptide = 0.0
add_Nterm_peptide = 0.0
add_Cterm_protein = 0.0
add_Nterm_protein = 0.0
add_G_glycine = 0.0000
add_A_alanine = 0.0000
add_S_serine = 0.0000
add_P_proline = 0.0000
add_V_valine = 0.0000
add_T_threonine = 0.0000
add_C_cysteine = 57.021464
add_L_leucine = 0.0000
add_I_isoleucine = 0.0000
add_N_asparagine = 0.0000
add_D_aspartic_acid = 0.0000
add_Q_glutamine = 0.0000
add_K_lysine = 0.0000
add_E_glutamic_acid = 0.0000
add_M_methionine = 0.0000
add_O_ornithine = 0.0000
add_H_histidine = 0.0000
add_F_phenylalanine = 0.0000
add_U_selenocysteine = 0.0000
add_R_arginine = 0.0000
add_Y_tyrosine = 0.0000
add_W_tryptophan = 0.0000
add_B_user_amino_acid = 0.0000
add_J_user_amino_acid = 0.0000
add_X_user_amino_acid = 0.0000
add_Z_user_amino_acid = 0.0000
#
# COMET_ENZYME_INFO _must_ be at the end of this parameters file
#
[COMET_ENZYME_INFO]
0.  No_enzyme              0      -           -
1.  Trypsin                1      KR          P
2.  Trypsin/P              1      KR          -
3.  Lys_C                  1      K           P
4.  Lys_N                  0      K           -
5.  Arg_C                  1      R           P
6.  Asp_N                  0      D           -
7.  CNBr                   1      M           -
8.  Glu_C                  1      DE          P
9.  PepsinA                1      FL          P
10. Chymotrypsin           1      FWYL        P
EOF

}

function write_comet_para_lowres_by() {
  cat <<EOF
# comet_version 2018.01 rev. X
# Comet MS/MS search engine parameters file.
# Everything following the '#' symbol is treated as a comment.
database_name = /srv/www/htdocs/FASTA/fgcz_uniprot-proteome_UP000005640.fasta
decoy_search = 1
num_threads = 32
peptide_mass_tolerance = 10.0
peptide_mass_units = 2
mass_type_parent = 1
mass_type_fragment = 1
precursor_tolerance_type = 1
isotope_error = 0
search_enzyme_number = 1
num_enzyme_termini = 2
allowed_missed_cleavage = 2
variable_mod01 = 15.9949 M 0 3 -1 0 0
variable_mod02 = 0.0 X 0 3 -1 0 0
variable_mod03 = 0.0 X 0 3 -1 0 0
variable_mod04 = 0.0 X 0 3 -1 0 0
variable_mod05 = 0.0 X 0 3 -1 0 0
variable_mod06 = 0.0 X 0 3 -1 0 0
variable_mod07 = 0.0 X 0 3 -1 0 0
variable_mod08 = 0.0 X 0 3 -1 0 0
variable_mod09 = 0.0 X 0 3 -1 0 0
max_variable_mods_in_peptide = 5
require_variable_mod = 0
## low res
fragment_bin_tol = 1.0005
fragment_bin_offset = 0.4
## For high-res MS/MS data
theoretical_fragment_ions = 1
use_A_ions = 0
use_B_ions = 1
use_C_ions = 0
use_X_ions = 0
use_Y_ions = 1
use_Z_ions = 0
use_NL_ions = 1
output_sqtstream = 0
output_sqtfile = 0
output_txtfile = 1
output_pepxmlfile = 0
output_percolatorfile = 0
output_outfiles = 0
print_expect_score = 1
num_output_lines = 1
show_fragment_ions = 0
sample_enzyme_number = 1
scan_range = 0
precursor_charge = 0
override_charge = 0
ms_level = 2
activation_method = ALL
digest_mass_range = 600.0
num_results = 100
skip_researching = 1
max_fragment_charge = 3
max_precursor_charge = 6
nucleotide_reading_frame = 0
clip_nterm_methionine = 1
spectrum_batch_size = 0
decoy_prefix = REV_
output_suffix = .comet
mass_offsets = 0.0
minimum_peaks = 10
minimum_intensity = 0
remove_precursor_peak = 0
remove_precursor_tolerance = 1.5
clear_mz_range = 0.0
add_Cterm_peptide = 0.0
add_Nterm_peptide = 0.0
add_Cterm_protein = 0.0
add_Nterm_protein = 0.0
add_G_glycine = 0.0000
add_A_alanine = 0.0000
add_S_serine = 0.0000
add_P_proline = 0.0000
add_V_valine = 0.0000
add_T_threonine = 0.0000
add_C_cysteine = 57.021464
add_L_leucine = 0.0000
add_I_isoleucine = 0.0000
add_N_asparagine = 0.0000
add_D_aspartic_acid = 0.0000
add_Q_glutamine = 0.0000
add_K_lysine = 0.0000
add_E_glutamic_acid = 0.0000
add_M_methionine = 0.0000
add_O_ornithine = 0.0000
add_H_histidine = 0.0000
add_F_phenylalanine = 0.0000
add_U_selenocysteine = 0.0000
add_R_arginine = 0.0000
add_Y_tyrosine = 0.0000
add_W_tryptophan = 0.0000
add_B_user_amino_acid = 0.0000
add_J_user_amino_acid = 0.0000
add_X_user_amino_acid = 0.0000
add_Z_user_amino_acid = 0.0000
#
# COMET_ENZYME_INFO _must_ be at the end of this parameters file
#
[COMET_ENZYME_INFO]
0.  No_enzyme              0      -           -
1.  Trypsin                1      KR          P
2.  Trypsin/P              1      KR          -
3.  Lys_C                  1      K           P
4.  Lys_N                  0      K           -
5.  Arg_C                  1      R           P
6.  Asp_N                  0      D           -
7.  CNBr                   1      M           -
8.  Glu_C                  1      DE          P
9.  PepsinA                1      FL          P
10. Chymotrypsin           1      FWYL        P
EOF

}

# ET
function write_comet_para_lowres_bcyz() {
  cat <<EOF
# comet_version 2018.01 rev. X
# Comet MS/MS search engine parameters file.
# Everything following the '#' symbol is treated as a comment.
database_name = /srv/www/htdocs/FASTA/fgcz_uniprot-proteome_UP000005640.fasta
decoy_search = 1
num_threads = 32
peptide_mass_tolerance = 10.0
peptide_mass_units = 2
mass_type_parent = 1
mass_type_fragment = 1
precursor_tolerance_type = 1
isotope_error = 0
search_enzyme_number = 1
num_enzyme_termini = 2
allowed_missed_cleavage = 2
variable_mod01 = 15.9949 M 0 3 -1 0 0
variable_mod02 = 0.0 X 0 3 -1 0 0
variable_mod03 = 0.0 X 0 3 -1 0 0
variable_mod04 = 0.0 X 0 3 -1 0 0
variable_mod05 = 0.0 X 0 3 -1 0 0
variable_mod06 = 0.0 X 0 3 -1 0 0
variable_mod07 = 0.0 X 0 3 -1 0 0
variable_mod08 = 0.0 X 0 3 -1 0 0
variable_mod09 = 0.0 X 0 3 -1 0 0
max_variable_mods_in_peptide = 5
require_variable_mod = 0
## low res
fragment_bin_tol = 1.0005
fragment_bin_offset = 0.4
## For high-res MS/MS data
theoretical_fragment_ions = 1
use_A_ions = 0
use_B_ions = 1
use_C_ions = 1
use_X_ions = 0
use_Y_ions = 1
use_Z_ions = 1
use_NL_ions = 1
output_sqtstream = 0
output_sqtfile = 0
output_txtfile = 1
output_pepxmlfile = 0
output_percolatorfile = 0
output_outfiles = 0
print_expect_score = 1
num_output_lines = 1
show_fragment_ions = 0
sample_enzyme_number = 1
scan_range = 0
precursor_charge = 0
override_charge = 0
ms_level = 2
activation_method = ALL
digest_mass_range = 600.0
num_results = 100
skip_researching = 1
max_fragment_charge = 3
max_precursor_charge = 6
nucleotide_reading_frame = 0
clip_nterm_methionine = 1
spectrum_batch_size = 0
decoy_prefix = REV_
output_suffix = .comet
mass_offsets = 0.0
minimum_peaks = 10
minimum_intensity = 0
remove_precursor_peak = 0
remove_precursor_tolerance = 1.5
clear_mz_range = 0.0
add_Cterm_peptide = 0.0
add_Nterm_peptide = 0.0
add_Cterm_protein = 0.0
add_Nterm_protein = 0.0
add_G_glycine = 0.0000
add_A_alanine = 0.0000
add_S_serine = 0.0000
add_P_proline = 0.0000
add_V_valine = 0.0000
add_T_threonine = 0.0000
add_C_cysteine = 57.021464
add_L_leucine = 0.0000
add_I_isoleucine = 0.0000
add_N_asparagine = 0.0000
add_D_aspartic_acid = 0.0000
add_Q_glutamine = 0.0000
add_K_lysine = 0.0000
add_E_glutamic_acid = 0.0000
add_M_methionine = 0.0000
add_O_ornithine = 0.0000
add_H_histidine = 0.0000
add_F_phenylalanine = 0.0000
add_U_selenocysteine = 0.0000
add_R_arginine = 0.0000
add_Y_tyrosine = 0.0000
add_W_tryptophan = 0.0000
add_B_user_amino_acid = 0.0000
add_J_user_amino_acid = 0.0000
add_X_user_amino_acid = 0.0000
add_Z_user_amino_acid = 0.0000
#
# COMET_ENZYME_INFO _must_ be at the end of this parameters file
#
[COMET_ENZYME_INFO]
0.  No_enzyme              0      -           -
1.  Trypsin                1      KR          P
2.  Trypsin/P              1      KR          -
3.  Lys_C                  1      K           P
4.  Lys_N                  0      K           -
5.  Arg_C                  1      R           P
6.  Asp_N                  0      D           -
7.  CNBr                   1      M           -
8.  Glu_C                  1      DE          P
9.  PepsinA                1      FL          P
10. Chymotrypsin           1      FWYL        P
EOF

}


function run_comet() {

  mkdir -p $scratch || exit 1

  local instrument=`echo $1 | sed -e 's/.*\(EXPLORIS_[12]\|QEXACTIVE_[12]\|LUMOS_[12]\|QEXACTIVEHF_[2]\|FUSION_[2]\).*/\1/'`
  local raw="${instrument}.`basename $1`"

  echo $raw | grep -f ${scratch}/blacklist.txt > /dev/null
  [ $? -eq 0 ] && { echo "${raw} found in blacklist"; return 0; }

  echo "${scratch}/${raw}"
  sleep 2
  [ -L ${scratch}/${raw} -o -f ${scratch}/${raw} ] || cp -av $1 ${scratch}/${raw} || exit 1


  cd $scratch || exit 1

  local mgf=`basename ${raw} raw`
  # echo "consider ${mgf}"

  ls | egrep "${mgf}.+mgf$" >/dev/null
  if [ $? -eq 1 ];
  then
    echo "generate MGF ... "

    echo "source('/home/cpanse/__projects/2019/20190225-comet/top5.R'); mgf_split('${scratch}/$raw')" \
      | R --no-save
   # nice -19 Rscript --vanilla ~cpanse/__checkouts/prx/R/mgf.R ${scratch}/${raw} $PWD/ 
  fi

  ls \
    | egrep "${mgf}.+mgf$" \
    | while read f;
    do
      cometresult=`basename ${f} .mgf`.comet.txt
      if [ ! -f $cometresult ];
      then
        COMETPARA=${PWD}/comet_lowres_by.para

        echo ${f} | egrep "ET.+.lowres.mgf$" >/dev/null \
        && COMETPARA=${PWD}/comet_lowres_bcyz.para

        echo ${f} | egrep "HCD.highres.mgf$" > /dev/null \
        && COMETPARA=${PWD}/comet_highres_by.para


      
set -x
        [ -f ${COMETPARA} ] || exit 1;
        cp ${COMETPARA} ${f}.comet.para
        nice -19 /usr/local/bin/comet.exe -P${COMETPARA} ${f}
set +x
      fi
    done
}

write_comet_para_highres_by > ${scratch}/comet_highres_by.para
write_comet_para_lowres_by > ${scratch}/comet_lowres_by.para
write_comet_para_lowres_bcyz > ${scratch}/comet_lowres_bcyz.para

export -f run_comet

MYPATTERN=".+(EXPLORIS_[12]|QEXACTIVE_[12]|LUMOS_[12]|QEXACTIVEHF_[2]|FUSION_[2]).+20[1234567][0-9].+(autoQC4L|autoQC03dda).+raw$"

ssh r35 "egrep '${MYPATTERN}' /srv/www/htdocs/Data2San/sync_LOGS/pfiles.txt" \
 | egrep "/202[0-9].*\.raw" \
 | tee pfiles_autoQC4L.txt \
 | awk -F';' '{print "/srv/www/htdocs/"$NF}' \
 | tail -n 10 \
 | parallel -P 4 run_comet 

sleep 1
test $? -eq 0 && R -q --no-site-file --no-save < /home/cpanse/src/gitlab.bfabric.org/proteomics/qc/inst/scripts/comet_autoQC03dda-post.R
exit $?
