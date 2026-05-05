#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
curl -L --fail --retry 3 -o 'E-MTAB-14704.idf.txt' 'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/704/E-MTAB-14704/Files/E-MTAB-14704.idf.txt'
curl -L --fail --retry 3 -o 'E-MTAB-14704.sdrf.txt' 'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/704/E-MTAB-14704/Files/E-MTAB-14704.sdrf.txt'
curl -L --fail --retry 3 -o 'P_SKMel29_AdDNp73beta_1.CEL' 'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/704/E-MTAB-14704/Files/P_SKMel29_AdDNp73beta_1.CEL'
curl -L --fail --retry 3 -o 'P_SKMel29_AdDNp73beta_2.CEL' 'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/704/E-MTAB-14704/Files/P_SKMel29_AdDNp73beta_2.CEL'
curl -L --fail --retry 3 -o 'P_SKMel29_AdDNp73beta_3.CEL' 'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/704/E-MTAB-14704/Files/P_SKMel29_AdDNp73beta_3.CEL'
curl -L --fail --retry 3 -o 'P_SKMel29_AdGFP_1.CEL' 'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/704/E-MTAB-14704/Files/P_SKMel29_AdGFP_1.CEL'
curl -L --fail --retry 3 -o 'P_SKMel29_AdGFP_2.CEL' 'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/704/E-MTAB-14704/Files/P_SKMel29_AdGFP_2.CEL'
curl -L --fail --retry 3 -o 'P_SKMel29_AdGFP_3.CEL' 'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/704/E-MTAB-14704/Files/P_SKMel29_AdGFP_3.CEL'
curl -L --fail --retry 3 -o 'P_SKMel29_AdTAp73alpha_1.CEL' 'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/704/E-MTAB-14704/Files/P_SKMel29_AdTAp73alpha_1.CEL'
curl -L --fail --retry 3 -o 'P_SKMel29_AdTAp73alpha_2.CEL' 'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/704/E-MTAB-14704/Files/P_SKMel29_AdTAp73alpha_2.CEL'
curl -L --fail --retry 3 -o 'P_SKMel29_AdTAp73alpha_3.CEL' 'https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/704/E-MTAB-14704/Files/P_SKMel29_AdTAp73alpha_3.CEL'
