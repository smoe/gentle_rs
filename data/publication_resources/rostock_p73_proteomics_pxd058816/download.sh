#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
curl -L --fail --retry 3 -o '190607_HFx_UL_BMFZ_Ctrl_1.msf' 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058816/190607_HFx_UL_BMFZ_Ctrl_1.msf'
curl -L --fail --retry 3 -o '190607_HFx_UL_BMFZ_Ctrl_1.raw' 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058816/190607_HFx_UL_BMFZ_Ctrl_1.raw'
curl -L --fail --retry 3 -o '190607_HFx_UL_BMFZ_Ctrl_2.msf' 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058816/190607_HFx_UL_BMFZ_Ctrl_2.msf'
curl -L --fail --retry 3 -o '190607_HFx_UL_BMFZ_Ctrl_2.raw' 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058816/190607_HFx_UL_BMFZ_Ctrl_2.raw'
curl -L --fail --retry 3 -o '190607_HFx_UL_BMFZ_TAp73_1.msf' 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058816/190607_HFx_UL_BMFZ_TAp73_1.msf'
curl -L --fail --retry 3 -o '190607_HFx_UL_BMFZ_TAp73_1.raw' 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058816/190607_HFx_UL_BMFZ_TAp73_1.raw'
curl -L --fail --retry 3 -o '190607_HFx_UL_BMFZ_TAp73_2.msf' 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058816/190607_HFx_UL_BMFZ_TAp73_2.msf'
curl -L --fail --retry 3 -o '190607_HFx_UL_BMFZ_TAp73_2.raw' 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058816/190607_HFx_UL_BMFZ_TAp73_2.raw'
curl -L --fail --retry 3 -o 'Human_uniprot_01_2017.fasta' 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058816/Human_uniprot_01_2017.fasta'
curl -L --fail --retry 3 -o 'checksum.txt' 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/08/PXD058816/checksum.txt'
