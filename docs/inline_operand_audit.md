# Inline ASCII Operand Audit

Date: 2026-05-14

This audit follows the architecture invariant in `docs/architecture.md` and
`docs/decisions.md`: low-latency sequence inspection should be state-optional
when an operation only needs sequence letters, optional topology, and an
optional local span.

Scope:

- `src/engine.rs` `Operation` variants, executed through `src/engine/ops/`.
- Inline-capable protocol support in `src/engine/protocol.rs` and relevant
  shared records in `crates/gentle-protocol/src/lib.rs`.
- This document tracks the audit state; follow-up implementation PRs should
  stay tiny and one operation at a time.

## Protocol Records

| Record | Current operand shape | Inline ASCII today | Classification |
| --- | --- | --- | --- |
| `SequenceScanTarget` | `SeqId { seq_id, span_* }` or `InlineSequence { sequence_text, topology, id_hint, span_* }` | yes | inline-ok |
| `RestrictionSiteScanReport` | report carries `target_kind`, `target_label`, source/scan spans, topology, rows | yes | inline-ok result carrier |
| TFBS scan/score reports | reports are built from `SequenceScanTarget` in the engine and carry target labels/spans | yes | inline-ok result carrier |
| `SequenceFeatureQuery` | `seq_id` plus feature/range/qualifier filters | no | inherently-stateful: queries persisted annotations |
| `FeatureExpertTarget` | feature/site/panel selector paired with `RenderFeatureExpertSvg.seq_id` | no | inherently-stateful: needs existing features, panels, or projections |

## Bucket 1: Inline-Ok Today

These operations already accept inline ASCII through `SequenceScanTarget`.

| Operation | Current operand shape | Inline ASCII today |
| --- | --- | --- |
| `RenderTfbsScoreTracksSvg` | `target: SequenceScanTarget`, motifs, score/render options, path | yes |
| `FindRestrictionSites` | `target: SequenceScanTarget`, enzymes, cut-geometry options, optional output path | yes |
| `SummarizeTfbsScoreTracks` | `target: SequenceScanTarget`, motifs, score options, optional output path | yes |
| `SummarizeTfbsTrackSimilarity` | `target: SequenceScanTarget`, anchor/candidate motifs, ranking options, optional output path | yes |
| `ScanTfbsHits` | `target: SequenceScanTarget`, motifs, thresholds, optional output path | yes |
| `AlignSequences` | `query`/`target: SequenceScanTarget`, legacy `query_seq_id`/`target_seq_id` accepted, scoring options | yes |

## Bucket 2: Should Accept Inline But Does Not Yet

These operations appear to be inspect/readout operations over sequence letters
where persisted state is a convenience rather than a biological requirement.

| Priority | Operation | Current operand shape | Smallest source change | Estimated source diff |
| --- | --- | --- | --- | --- |
| 1 | `RenderSequenceSvg` | `seq_id`, `mode`, `path` | Replace `seq_id` with `target: SequenceScanTarget`; for inline input render an unannotated temporary `DNAsequence` with default display state and use the target label in messages. | 70-120 LOC |
| 2 | `RenderRnaStructureSvg` | `seq_id`, `path` | Replace `seq_id` with an inline-capable target, reusing the RNA renderer against a temporary sequence; confirm inline `U`/RNA text normalization before exposing the shell route. | 60-100 LOC |
| 3 | `RenderTfbsScoreTrackCorrelationSvg` | `seq_id`, start/end span, motifs, score/correlation options, path | Replace `seq_id` plus start/end with `target: SequenceScanTarget`; the handler already wraps the old fields into `SequenceScanTarget::SeqId` before calling the inline-capable scorer. | 30-60 LOC |
| 4 | `ComputeDotplot` | `seq_id`, optional `reference_seq_id`, query/reference spans, dotplot options, store id | Replace query and optional reference sequence ids/spans with `SequenceScanTarget` operands; omit reference annotation for inline references and keep stored dotplot output explicit. | 120-220 LOC |
| 5 | `ComputeFlexibilityTrack` | `seq_id`, optional span, model/window options, store id | Replace `seq_id` plus span with `target: SequenceScanTarget`; decide whether inline results are returned report-only or stored as detached tracks keyed by `id_hint`. | 90-160 LOC |

Notes:

- `FindRestrictionSites`, `ScanTfbsHits`, and TFBS score summaries already cover
  the prompt's high-frequency restriction/motif examples.
- The current enum does not expose a dedicated stateless GC-content, ORF-scan,
  or simple translate-to-protein operation. Existing related operations either
  mutate project state (`DeriveProteinSequences`, `Reverse*`, `ExtractRegion`)
  or require annotations/reports. Those are outside this `seq_id_or_inline`
  audit unless a later PR adds a new pure readout operation.

## Bucket 3: Inherently Stateful Or Not A Sequence-Letter Operand

These operations either mutate stored state, require containers/lineage/reports,
use persisted annotations or analysis stores, operate on external catalogs/files,
or do not inspect caller-supplied sequence letters.

| Operation | Current operand shape | Inline ASCII today | Classification |
| --- | --- | --- | --- |
| `LoadFile` | `path`, optional `as_id` | N/A | inherently-stateful import |
| `SaveFile` | `seq_id`, `path`, format | no | inherently-stateful sequence export |
| `RenderDotplotSvg` | `seq_id`, `dotplot_id`, render options, path | no | inherently-stateful: needs stored dotplot analysis |
| `RenderFeatureExpertSvg` | `seq_id`, `FeatureExpertTarget`, path | no | inherently-stateful: needs features/sites/panels/projections |
| `RenderIsoformArchitectureSvg` | `seq_id`, `panel_id`, optional expression TSV, path | no | inherently-stateful: needs imported isoform panel |
| `RenderLineageSvg` | path | N/A | inherently-stateful lineage view |
| `RenderPoolGelSvg` | sequence/container/arrangement ids, ladders, gel conditions, path | no | inherently-stateful pool/container view |
| `RenderProteinGelSvg` | `report_id`, ladders, path | no | inherently-stateful report view |
| `RenderProteinGelReportsSvg` | `report_ids`, ladders, path | no | inherently-stateful report view |
| `RenderProteaseDigestGelSvg` | optional protein `seq_id` or `report_id`, proteases, path | no | inherently-stateful protein/report gel view |
| `RenderProtein2dGelSvg` | `report_id`, ladders, path | no | inherently-stateful report view |
| `RenderProtocolCartoonSvg` | protocol kind, path | N/A | no sequence operand |
| `RenderProtocolCartoonTemplateSvg` | template path, output path | N/A | external template render |
| `ValidateProtocolCartoonTemplate` | template path | N/A | external template validation |
| `RenderProtocolCartoonTemplateWithBindingsSvg` | template path, bindings path, output path | N/A | external template render |
| `ExportProtocolCartoonTemplateJson` | protocol kind, path | N/A | no sequence operand |
| `ApplyGibsonAssemblyPlan` | plan JSON | N/A | plan/state operation |
| `CreateArrangementSerial` | container ids, arrangement metadata, ladders | N/A | inherently-stateful arrangement operation |
| `SetArrangementLadders` | arrangement id, ladders | N/A | inherently-stateful arrangement operation |
| `SetContainerDeclaredContentsExclusive` | container id, flag | N/A | inherently-stateful container operation |
| `CreateRackFromArrangement` | arrangement id, rack/profile metadata | N/A | inherently-stateful rack operation |
| `PlaceArrangementOnRack` | arrangement id, rack id | N/A | inherently-stateful rack operation |
| `MoveRackPlacement` | rack id and coordinates | N/A | inherently-stateful rack operation |
| `MoveRackSamples` | rack id and coordinates | N/A | inherently-stateful rack operation |
| `MoveRackArrangementBlocks` | rack id, arrangement ids, coordinate | N/A | inherently-stateful rack operation |
| `SetRackProfile` | rack id, profile | N/A | inherently-stateful rack operation |
| `ApplyRackTemplate` | rack id, template | N/A | inherently-stateful rack operation |
| `SetRackFillDirection` | rack id, fill direction | N/A | inherently-stateful rack operation |
| `SetRackProfileCustom` | rack id, row/column dimensions | N/A | inherently-stateful rack operation |
| `SetRackBlockedCoordinates` | rack id, blocked coordinates | N/A | inherently-stateful rack operation |
| `ExportRackLabelsSvg` | rack id, optional arrangement id, preset, path | N/A | inherently-stateful rack export |
| `ExportRackFabricationSvg` | rack id, template, path | N/A | inherently-stateful rack export |
| `ExportRackIsometricSvg` | rack id, template, path | N/A | inherently-stateful rack export |
| `ExportRackOpenScad` | rack id, template, path | N/A | inherently-stateful rack export |
| `ExportRackCarrierLabelsSvg` | rack id, optional arrangement id, template/preset, path | N/A | inherently-stateful rack export |
| `ExportRackSimulationJson` | rack id, template, path | N/A | inherently-stateful rack export |
| `ExportDnaLadders` | path, optional name filter | N/A | catalog export |
| `ExportRnaLadders` | path, optional name filter | N/A | catalog export |
| `ExportPool` | sequence ids, pool metadata, path | no | inherently-stateful sequence-pool export |
| `ExportProcessRunBundle` | path, optional run id | N/A | inherently-stateful run export |
| `ExportLabAssistantInstructions` | path, optional run id/title/audience | N/A | inherently-stateful run export |
| `PrepareGenome` | genome id, catalog/cache paths, timeout | N/A | external genome cache operation |
| `ExtractGenomeRegion` | genome id, chromosome/range, output id, annotation options | N/A | external genome extraction into state |
| `ExtractGenomeGene` | genome id, gene query, output id, annotation options | N/A | external genome extraction into state |
| `ExtractGenomePromoterSlice` | genome id, gene/transcript query, window, output id | N/A | external genome extraction into state |
| `ExtendGenomeAnchor` | stored `seq_id`, side/length, genome cache options | no | inherently-stateful: needs stored genome anchor metadata |
| `VerifyGenomeAnchor` | stored `seq_id`, genome cache options | no | inherently-stateful: needs stored genome anchor metadata |
| `ImportGenomeBedTrack` | `seq_id`, track path/options | no | inherently-stateful annotation import |
| `ImportGenomeBigWigTrack` | `seq_id`, track path/options | no | inherently-stateful annotation import |
| `ImportGenomeVcfTrack` | `seq_id`, track path/options | no | inherently-stateful annotation import |
| `ProjectMicroarrayTrack` | `seq_id`, manifest path, contrast/filter options | no | inherently-stateful annotation/evidence projection |
| `ProjectGenomeInterval` | source/target genome ids, projection path, interval | N/A | external genome projection |
| `ListCutRunDatasets` | optional filter/catalog path | N/A | catalog query |
| `ShowCutRunDatasetStatus` | dataset id, catalog/cache paths | N/A | catalog/cache query |
| `PrepareCutRunDataset` | dataset id, catalog/cache paths | N/A | external dataset preparation |
| `ProjectCutRunDataset` | `seq_id`, dataset id, projection options | no | inherently-stateful evidence projection |
| `InterpretCutRunReads` | `seq_id`, read/dataset inputs, align config, report/checkpoint ids | no | inherently-stateful read interpretation |
| `ListCutRunReadReports` | optional `seq_id` | no | inherently-stateful report listing |
| `ShowCutRunReadReport` | report id | N/A | inherently-stateful report lookup |
| `ExportCutRunReadCoverage` | report id, coverage kind, path | N/A | inherently-stateful report export |
| `InspectCutRunRegulatorySupport` | `seq_id`, dataset/read report ids, promoter window options | no | inherently-stateful evidence/feature inspection |
| `ImportIsoformPanel` | `seq_id`, panel path/id, strict flag | no | inherently-stateful panel import |
| `ImportUniprotSwissProt` | file path, entry id | N/A | external entry import |
| `FetchUniprotSwissProt` | query, entry id | N/A | external entry fetch |
| `FetchEnsemblGene` | query/species, entry id | N/A | external entry fetch |
| `FetchEnsemblRegion` | species/chromosome/range, output id | N/A | external sequence fetch into state |
| `FetchEnsemblProtein` | query, entry id | N/A | external entry fetch |
| `FetchGenBankAccession` | accession, optional id | N/A | external sequence fetch into state |
| `FetchDbSnpRegion` | rs id, genome id/cache, output id, annotation options | N/A | external sequence fetch into state |
| `FetchUniprotLinkedGenBank` | UniProt entry/accession, output id | N/A | external sequence fetch into state |
| `ImportUniprotEntrySequence` | entry id, output id | N/A | persisted entry import |
| `ImportEnsemblGeneSequence` | entry id, output id | N/A | persisted entry import |
| `ImportEnsemblProteinSequence` | entry id, output id | N/A | persisted entry import |
| `ProjectUniprotToGenome` | `seq_id`, entry id, projection/transcript ids | no | inherently-stateful projection |
| `QueryProteinResidueGenomicCoordinates` | `seq_id`, transcript/residue range | no | inherently-stateful annotation/projection query |
| `AuditUniprotProjectionConsistency` | projection/report/entry ids | N/A | inherently-stateful projection audit |
| `AuditUniprotProjectionParity` | projection/report/entry ids | N/A | inherently-stateful projection audit |
| `ImportBlastHitsTrack` | `seq_id`, BLAST hit records, track options | no | inherently-stateful annotation import |
| `DigestContainer` | container id, enzymes, output prefix | N/A | inherently-stateful container mutation |
| `MergeContainersById` | container ids, output prefix | N/A | inherently-stateful container mutation |
| `LigationContainer` | container id, protocol/output options | N/A | inherently-stateful container mutation |
| `FilterContainerByMolecularWeight` | container id, size filter, output options | N/A | inherently-stateful container mutation |
| `Digest` | input `SeqId`, enzymes, output prefix | no | inherently-stateful: materializes fragments/lineage |
| `QueryRepeatAnnotations` | genome id, RepeatMasker path, filter, path | N/A | external annotation index query |
| `QueryRepeatOverlaps` | `seq_id`, repeat index path, span, limit, path | no | inherently-stateful sequence-coordinate query |
| `MaterializeRepeatFeatures` | `seq_id`, repeat index path, materialization options | no | inherently-stateful feature mutation |
| `BuildRepeatEnvironmentCohort` | genome id, repeat path, flank/filter options | N/A | external genome cohort construction |
| `SummarizeWindowCohortTfbs` | cohort report, motifs, score/cache options | N/A | report/cohort summarization |
| `Ligation` | input `SeqId`s, protocol/output options | no | inherently-stateful: materializes sequence/lineage |
| `MergeContainers` | input `SeqId`s, output prefix | no | inherently-stateful: materializes sequence/container |
| `Pcr` | template `SeqId`, primer strings, output options | no | inherently-stateful: materializes PCR product |
| `PcrAdvanced` | template `SeqId`, primer specs, output options | no | inherently-stateful: materializes PCR product |
| `PcrMutagenesis` | template `SeqId`, primer specs, mutations, output options | no | inherently-stateful: materializes edited product |
| `DesignPrimerPairs` | template `SeqId`, ROI, constraints, report id | no | inherently-stateful report tied to template coordinates |
| `DesignInsertionPrimerPairs` | template `SeqId`, insertion intent, constraints, report id | no | inherently-stateful report tied to template coordinates |
| `ExportPrimerDesignReport` | report id, path | N/A | inherently-stateful report export |
| `AssessPrimerPairSpecificity` | primer/report inputs, target genome/cache, policy, path | N/A | external genome specificity report |
| `PrepareRestrictionCloningPcrHandoff` | template/vector `SeqId`s, primer report, enzymes | no | inherently-stateful handoff over persisted designs |
| `PcrOverlapExtensionMutagenesis` | template `SeqId`, edit span, insert sequence, constraints | no | inherently-stateful: materializes edited product |
| `DesignQpcrAssays` | template `SeqId`, ROI, constraints, optional transcript targeting | no | inherently-stateful report tied to template/transcripts |
| `TestCdnaPcr` | `seq_id`, source feature id, primers, transcript/report options | no | inherently-stateful transcript-feature assay |
| `TestCdnaQpcr` | `seq_id`, source feature id, primers/probe, transcript/report options | no | inherently-stateful transcript-feature assay |
| `BuildTranscriptQpcrPanel` | `seq_id`, source feature id, qPCR report id | no | inherently-stateful transcript panel |
| `TestCdnaQpcrFasta` | cDNA FASTA paths, primers/probe, output paths | N/A | external FASTA assay check |
| `DeriveTranscriptSequences` | `seq_id`, feature ids, scope, output prefix | no | inherently-stateful annotation-derived materialization |
| `PlanExonSkippedIsoform` | `seq_id`, transcript feature id, criteria, plan id | no | inherently-stateful splicing plan |
| `MaterializeExonSkippedIsoform` | plan id, candidate ids, output options | N/A | inherently-stateful plan materialization |
| `DeriveProteinSequences` | `seq_id`, feature ids/query, scope, output/report ids | no | inherently-stateful annotation-derived translation |
| `ReverseTranslateProteinSequence` | protein `seq_id`, output id, speed/translation options | no | inherently-stateful protein-to-DNA materialization |
| `ProteaseDigestProteinSequence` | protein `seq_id`, proteases, output/materialization options | no | inherently-stateful protein report/materialization |
| `BuildProteinToDnaHandoffReasoning` | DNA/protein `SeqId`s, projection/feature ids, ranking options | no | inherently-stateful reasoning graph |
| `ComputeDotplotOverlay` | owner/reference `SeqId`s, overlay query specs, store id | no | inherently-stateful dotplot overlay |
| `DeriveSplicingReferences` | `seq_id`, span, seed feature id, scope/output prefix | no | inherently-stateful splicing-feature materialization |
| `ImportSequencingTrace` | trace path/id, optional `seq_id` | N/A | external trace import |
| `ListSequencingTraces` | optional `seq_id` | no | inherently-stateful trace listing |
| `ShowSequencingTrace` | trace id | N/A | inherently-stateful trace lookup |
| `SuggestSequencingPrimers` | expected `SeqId`, primer `SeqId`s, report options | no | inherently-stateful primer/sequence inspection |
| `ConfirmConstructReads` | expected/baseline/read `SeqId`s, traces, targets, alignment options | no | inherently-stateful confirmation report |
| `ListSequencingConfirmationReports` | optional expected `SeqId` | no | inherently-stateful report listing |
| `ShowSequencingConfirmationReport` | report id | N/A | inherently-stateful report lookup |
| `ExportSequencingConfirmationReport` | report id, path | N/A | inherently-stateful report export |
| `ExportSequencingConfirmationSupportTsv` | report id, path | N/A | inherently-stateful report export |
| `ReadAcquireStatus` | manifest/cache/work dirs | N/A | external read-acquisition state |
| `ReadAcquirePrepare` | manifest/cache/work dirs, acquisition options | N/A | external read acquisition |
| `ReadAcquireInspect` | SRA accession/cache/work dirs | N/A | external read acquisition |
| `ReadAcquireCancel` | SRA accession/cache/work dirs | N/A | external read acquisition |
| `InterpretRnaReads` | `seq_id`, seed feature id, read input, profiles/filters, report/checkpoint ids | no | inherently-stateful RNA-read interpretation |
| `AlignRnaReadReport` | report id, selection, align override | N/A | inherently-stateful report refinement |
| `PreflightRnaReadIsoforms` | `seq_id`, seed feature id, control FASTAs, filters | no | inherently-stateful isoform preflight |
| `ListRnaReadReports` | optional `seq_id` | no | inherently-stateful report listing |
| `ShowRnaReadReport` | report id | N/A | inherently-stateful report lookup |
| `SummarizeRnaReadGeneSupport` | report id, gene ids, selection/rule, path | N/A | inherently-stateful report summarization |
| `InspectRnaReadGeneSupport` | report id, gene ids, selection/rule/filter, path | N/A | inherently-stateful report summarization |
| `RunRnaReadBatchMap` | manifest path, `seq_id`, seed feature id, batch/read options | no | inherently-stateful batch interpretation |
| `SummarizeTfbsRegion` | `seq_id`, focus/context spans, occurrence filters, path | no | inherently-stateful: summarizes persisted TFBS features |
| `SummarizeMultiGenePromoterTfbs` | genome id, gene queries, motifs, windows/cache/path | N/A | external genome/gene summary |
| `RenderMultiGenePromoterTfbsSvg` | genome id, gene queries, motifs, windows/cache/path | N/A | external genome/gene render |
| `SummarizeJasparEntries` | motifs, random-sequence presentation settings, path | N/A | catalog/motif summary |
| `ResolveTfQueries` | TF query strings, path | N/A | catalog query |
| `BenchmarkJasparRegistry` | random sequence settings, path | N/A | catalog benchmark |
| `ListJasparCatalog` | filter/limit/remote metadata flags/path | N/A | catalog query |
| `SyncJasparRemoteMetadata` | motifs/filter/limit/path | N/A | catalog sync |
| `InspectJasparEntry` | motif, random-sequence settings, remote metadata flags/path | N/A | catalog entry inspection |
| `AnnotatePromoterWindows` | input `SeqId`, gene/transcript labels, window/collapse options | no | inherently-stateful feature mutation |
| `SummarizeVariantPromoterContext` | input `SeqId`, variant/gene/transcript labels, windows/path | no | inherently-stateful variant/feature summary |
| `SummarizeAlternativePromoterComparison` | input `SeqId`, gene/transcript labels, windows/path | no | inherently-stateful promoter-feature summary |
| `SummarizePromoterEvidenceMatrix` | input `SeqId`, gene/transcript labels, windows/path | no | inherently-stateful feature/evidence summary |
| `SummarizeIsoformPromoterComparison` | input `SeqId`, gene/transcript labels, windows/path | no | inherently-stateful isoform/feature summary |
| `SummarizePromoterExpressionEvidence` | input `SeqId`, gene/transcript labels, expression rows, path | no | inherently-stateful promoter/expression summary |
| `ExportPromoterArtifactManifest` | input `SeqId`, gene label, artifact entries, path | no | inherently-stateful manifest export |
| `SuggestPromoterReporterFragments` | input `SeqId`, variant/gene/transcript labels, window/candidate options | no | inherently-stateful feature-aware fragment suggestion |
| `MaterializeVariantAllele` | input `SeqId`, variant selector, allele choice, output id | no | inherently-stateful sequence mutation |
| `ExportRnaReadReport` | report id, path | N/A | inherently-stateful report export |
| `ExportRnaReadHitsFasta` | report id, selection/subset, path | N/A | inherently-stateful report export |
| `ExportRnaReadSampleSheet` | optional `seq_id`, report ids, gene ids, path | no | inherently-stateful report/sample export |
| `ExportRnaReadTargetQuality` | report id, gene ids, rule, path | N/A | inherently-stateful report export |
| `ExportRnaReadExonPathsTsv` | report id, selection/subset, path | N/A | inherently-stateful report export |
| `ExportRnaReadExonAbundanceTsv` | report id, selection/subset, path | N/A | inherently-stateful report export |
| `ExportRnaReadScoreDensitySvg` | report id, density options, path | N/A | inherently-stateful report render |
| `ExportRnaReadAlignmentsTsv` | report id, selection/limit/subset, path | N/A | inherently-stateful report export |
| `ExportRnaReadAlignmentDotplotSvg` | report id, selection/max points, path | N/A | inherently-stateful report render |
| `MaterializeRnaReadHitSequences` | report id, selection/subset, output prefix | N/A | inherently-stateful report materialization |
| `ExtractRegion` | input `SeqId`, range, output id | no | inherently-stateful: materializes extracted sequence |
| `ExtractAnchoredRegion` | input `SeqId`, anchor/direction/constraints/primers/output options | no | inherently-stateful candidate materialization |
| `GenerateCandidateSetBetweenAnchors` | set name, `seq_id`, anchors, length/step/limit | no | inherently-stateful candidate-set creation |
| `SelectCandidate` | input `SeqId`, criterion, output id | no | inherently-stateful candidate materialization |
| `FilterByMolecularWeight` | input `SeqId`s, size/error filters, output options | no | inherently-stateful sequence filtering/materialization |
| `FilterByDesignConstraints` | input `SeqId`s, quality filters, output options | no | inherently-stateful sequence filtering/materialization |
| `GenerateCandidateSet` | set name, `seq_id`, feature/distance filters, length/step/limit | no | inherently-stateful candidate-set creation |
| `DeleteCandidateSet` | set name | N/A | inherently-stateful candidate-set mutation |
| `UpsertGuideSet` | guide set id, guide candidates | N/A | inherently-stateful guide-set mutation |
| `DeleteGuideSet` | guide set id | N/A | inherently-stateful guide-set mutation |
| `FilterGuidesPractical` | guide set id, filter config, output guide set id | N/A | inherently-stateful guide-set mutation |
| `GenerateGuideOligos` | guide set id, template id, output oligo set id/options | N/A | inherently-stateful guide-set mutation |
| `ExportGuideOligos` | guide/oligo set ids, format/plate/path | N/A | inherently-stateful guide export |
| `ExportGuideProtocolText` | guide/oligo set ids, path, QC flag | N/A | inherently-stateful guide export |
| `ExportFeaturesBed` | `SequenceFeatureQuery`, coordinate/export options, path | no | inherently-stateful annotation export |
| `InspectSequenceContextView` | `seq_id`, mode/viewport/classes/coordinate options | no | inherently-stateful annotation/display inspection |
| `ExportSequenceContextBundle` | `seq_id`, mode/viewport/export options, output dir | no | inherently-stateful annotation/display export |
| `ScoreCandidateSetExpression` | set name, metric, expression | N/A | inherently-stateful candidate-set scoring |
| `ScoreCandidateSetDistance` | set name, metric, feature filters | N/A | inherently-stateful candidate-set scoring |
| `FilterCandidateSet` | input/output set names, metric filters | N/A | inherently-stateful candidate-set filtering |
| `CandidateSetOp` | operator, left/right/output set names | N/A | inherently-stateful candidate-set mutation |
| `ScoreCandidateSetWeightedObjective` | set name, metric, objective terms | N/A | inherently-stateful candidate-set scoring |
| `TopKCandidateSet` | input/output set names, metric/k/tie policy | N/A | inherently-stateful candidate-set filtering |
| `ParetoFrontierCandidateSet` | input/output set names, objectives/options | N/A | inherently-stateful candidate-set filtering |
| `UpsertWorkflowMacroTemplate` | template metadata/ports/script | N/A | inherently-stateful macro catalog mutation |
| `DeleteWorkflowMacroTemplate` | template name | N/A | inherently-stateful macro catalog mutation |
| `UpsertCandidateMacroTemplate` | template metadata/script | N/A | inherently-stateful macro catalog mutation |
| `DeleteCandidateMacroTemplate` | template name | N/A | inherently-stateful macro catalog mutation |
| `Reverse` | input `SeqId`, output id | no | inherently-stateful: materializes derived sequence |
| `Complement` | input `SeqId`, output id | no | inherently-stateful: materializes derived sequence |
| `ReverseComplement` | input `SeqId`, output id | no | inherently-stateful: materializes derived sequence |
| `Branch` | input `SeqId`, output id | no | inherently-stateful lineage branch |
| `SetDisplayVisibility` | display target, visible flag | N/A | inherently-stateful display setting |
| `SetLinearViewport` | viewport start/span | N/A | inherently-stateful display setting |
| `SetTopology` | `seq_id`, circular flag | no | inherently-stateful sequence mutation |
| `RecomputeFeatures` | `seq_id` | no | inherently-stateful feature mutation |
| `SetParameter` | parameter name/value | N/A | inherently-stateful parameter setting |
| `AnnotateTfbs` | `seq_id`, motifs, thresholds, clear/max options | no | inherently-stateful feature mutation |
