# `test_files/fixtures/sequencing_confirmation`

Public benchmark shortlist and provenance plan for sequencing-confirmation
fixtures.

Status:

- This directory currently contains provenance notes only.
- No sequencing-confirmation fixture payloads are committed here yet.
- The first committed pack should stay small and deterministic, then layer real
  public read evidence on top for manual/exploratory benchmarking.

## Priority Order

1. `U13852` / `pGEX-3X` derived synthetic read pack for deterministic CI
2. `PRJNA1066256` / `SRR27605537` as a real public plasmid-read benchmark
3. Biopython `Tests/Abi/*.ab1` fixtures for ABI/AB1 parser intake

## 1. Deterministic phase-1 baseline: `U13852` / `pGEX-3X`

- Origin:
  - ENA accession `U13852`, `pGEX-3X cloning vector, complete sequence`
  - already represented in-repo as `test_files/pGEX-3X.embl`
- Why this is the top choice:
  - public and compact
  - expected construct sequence is exact and already part of GENtle parser
    coverage
  - ideal for tiny committed fixtures that lock the confirmation verdict logic
    before raw trace intake begins
- Planned fixture role:
  - perfect forward read
  - reverse-complement read
  - truncated read
  - contradicted junction or indel read
- Deterministic retrieval:
  - fetch ENA EMBL export for `U13852`
  - example source URL:
    `https://www.ebi.ac.uk/ena/browser/api/embl/U13852`
- Planned deterministic self-creation:
  - derive short synthetic read sequences from fixed coordinate windows on the
    `U13852` sequence
  - freeze the exact windows and mutation rules in this manifest when the
    synthetic reads are committed
- Planned GENtle usage:
  - `src/engine/tests.rs` confirmation-classification tests
  - `src/engine_shell/tests.rs` CLI/report export parity tests
  - future GUI specialist smoke coverage once confirmation inspection lands

## 2. Real public read benchmark: `PRJNA1066256` / `SRR27605537`

- Origin:
  - ENA study `PRJNA1066256`, `Recombination Plasmid Sequencing`
  - study description:
    `Whole plasmid nanopore sequencing of plasmids containing both the Bxb1 SSR and it's attP/attB sites.`
  - representative run:
    - `SRR27605537`
    - experiment title:
      `PromethION sequencing: Whole plasmid sequencing of PV521 purified from Syn61-deltaT`
- Why this is the best current real-read candidate:
  - explicitly whole-plasmid sequencing
  - small enough read count (`694`) to be practical for manual/exploratory
    confirmation benchmarking
  - closer to the intended construct-verification story than generic RNA or
    genome read sets
- Deterministic retrieval:
  - study metadata:
    `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1066256&result=study&fields=study_accession,study_alias,study_title,study_description,center_name,first_public,last_updated&format=json`
  - run metadata:
    `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR27605537&result=read_run&fields=run_accession,study_accession,study_title,experiment_title,sample_accession,instrument_platform,library_strategy,read_count,base_count,fastq_ftp,submitted_ftp&format=json`
  - FASTQ:
    `ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR276/037/SRR27605537/SRR27605537_1.fastq.gz`
- Current limitation:
  - this is not yet a committed CI fixture because the exact expected plasmid
    reference sequence for `PV521` has not been frozen in-repo
  - treat it as the first manual/external benchmark after the phase-1 synthetic
    pack is stable
- Planned GENtle usage:
  - exploratory confirmation runs against a real plasmid-read set
  - future read-downsampling experiments for compact public benchmark
    derivatives if redistribution/provenance remains acceptable

## 3. ABI/AB1 parser fixtures: Biopython `Tests/Abi`

- Origin:
  - Biopython repository, commit
    `d76a5581ca3317ec850d72450c458ec4e34e2583`
  - source tree:
    `Tests/Abi/`
- Recommended positive fixtures:
  - `310.ab1`
  - `3100.ab1`
  - `3730.ab1`
  - `A6_1-DB3.ab1`
  - `no_smpl1.ab1`
- Useful edge fixtures:
  - `empty.ab1`
  - `fake.ab1`
  - `nonascii_encoding.ab1`
  - `test.fsa`
- Why these are suitable:
  - already used as public ABI parser fixtures in a mature bioinformatics
    library
  - small and purpose-built for format handling rather than broad biology
    semantics
  - good nostalgia bridge to `gentle-m` `ABItype.cpp`/`SCFtype.cpp` work
- Deterministic retrieval:
  - raw file pattern pinned to commit:
    `https://raw.githubusercontent.com/biopython/biopython/d76a5581ca3317ec850d72450c458ec4e34e2583/Tests/Abi/<FILE>`
  - example:
    `https://raw.githubusercontent.com/biopython/biopython/d76a5581ca3317ec850d72450c458ec4e34e2583/Tests/Abi/310.ab1`
- Planned GENtle usage:
  - `src/genbank.rs` or new trace-parser module tests once ABI/AB1 import work
    starts
  - negative/edge parsing regressions before trace-aware confirmation is added
- Licensing note:
  - preserve upstream attribution when copying individual fixtures from
    Biopython into GENtle
  - recheck redistribution comfort at commit time, even though the fixtures are
    public in the upstream repository

## Not The Right Primary Benchmark

- `test_files/mapping/SRR32957124_10000.fasta.gz` remains useful for RNA-read
  report plumbing, but it is not an expected-construct confirmation dataset and
  should not be the anchor benchmark for this workflow.

## Adoption Plan

1. Commit the tiny `U13852`-derived synthetic pack first and freeze its exact
   derivation rules here.
2. Keep `PRJNA1066256` as the first real public benchmark for exploratory
   confirmation runs outside the default CI fixture set.
3. Bring in the smallest useful subset of Biopython ABI fixtures only when
   raw-trace import begins.
