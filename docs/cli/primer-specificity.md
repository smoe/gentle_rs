# Reference-Bound Primer Specificity

<!-- manual-counterparts:start -->
Related GUI reference: [Primer, qPCR, and Sequence-Specific RT Design](../gui/primer-qpcr-and-sequence-specific-rt-design.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

## Standalone Reference-Bound Primer Specificity

[Tutorial 04.09](../tutorial/04-09_reference_bound_primer_specificity.md) explains
why genomic carryover, transcriptome cross-amplification and unavailable evidence
must not be collapsed. Assess one saved pair or explicit full oligos with declared
annealing/tail boundaries against 1-8 explicit already prepared references.

```sh
gentle_cli --state study.gentle.json primers specificity-multi-handoff @request.json NEW_BUNDLE_DIR
# Inspect the handoff. A trusted external runner executes its structured commands
# and fills execution_manifest.json with actual exits and retained-byte hashes.
gentle_cli --state study.gentle.json primers specificity-multi-import NEW_BUNDLE_DIR/handoff.json NEW_BUNDLE_DIR/execution_manifest.json
gentle_cli --state study.gentle.json primers specificity-multi-list
gentle_cli --state study.gentle.json primers specificity-multi-show SUMMARY_ID
```

Required unavailable references fail before files; optional unavailable references
are unassessed. Preparation runs no BLAST and never prepares/downloads resources.
Import is an approved mutation; it validates commands, raw bytes, pair/policy,
source snapshot, complete searches and current reference/annotation fingerprints.
The two dimensions report `pass|fail|incomplete|not_requested|not_required` separately.
An imported summary is not automatic panel-readiness evidence or order approval.
Show/list inspect immutable history without database probes. Refer to the
[protocol](../protocol.md#standalone-multi-reference-primer-specificity) for schemas,
MCP tools, applicability and non-claims; the same commands work in GUI Shell.

For bounded import diagnostics only, `GENTLE_PRIMER_IMPORT_TRACE=1` emits
fixed-phase wall-clock JSON lines on stderr (at most 512 events per import).
It is off by default and does not add timings to stdout, reports or state;
completion events are not biological pass claims. Capture stderr to a local
file and follow the [privilege-free diagnostic procedure](../glen_post_blast_import_phase_trace.md).

### Inspect Retained Panel Specificity

```sh
gentle_cli --state study.gentle.json primers transcript-assay-specificity-matrix PANEL_REPORT_ID
```

This read-only command displays every retained selection receipt separately by
panel design, exact reference and policy, including active/superseded history,
per-assay binding diagnostics and unassessed current-design dimensions. It does
not search, probe databases, write evidence, choose a reference or return an
aggregate panel pass. Same-oligo standalone summaries remain unlinked citations,
not panel evidence. GUI Shell and MCP `transcript_assay_specificity_matrix` use
the same result; see [tutorial 04.09](../tutorial/04-09_reference_bound_primer_specificity.md#inspect-the-panels-retained-evidence-matrix).


- Shared shell:
  - Commands:
    - Primer specificity confirmation notes (`primers specificity`):
      - returns and persists `gentle.primer_specificity_report.v4`; repeated
        assessment of the same biological inputs and policy reuses one stable
        content-derived `report_id`
      - `--full-alignment disabled|best-effort|required` controls end-to-end
        realignment of complete primers after BLAST locates candidate loci.
        `best-effort` is the default; `required` makes unavailable alignment
        evidence incomplete rather than falling back to a partial local HSP
      - `--reviewed-off-target-allowlist JSON_ARRAY_OR_@FILE` accepts exact reviewed
        products bound to target space, subject, optional coordinates,
        reviewer, reason, and evidence hash. An allowance remains visibly
        off-target and cannot repair an incomplete search
      - `primers specificity-alignment-html REPORT_ID OUTPUT.html` projects
        the stored complete alignments without rerunning BLAST
      - `collections run primer-specificity` is the engine-owned `map`
        continuation for a persisted gene-set resolution or a set of project
        sequences. It calls this same assessment once per resolved primer
        report and returns `gentle.collection_operation.v1` with typed
        per-member outcomes and produced specificity-report ids
      - logical gene-set members require explicit repeated
        `--member-report MEMBER_ID=PRIMER_REPORT_ID` bindings. A gene symbol
        does not identify one assay and is never used to guess a report.
        Project-sequence members may omit a binding only when exactly one
        stored primer-design report targets that sequence; zero or multiple
        matches are retained as member failures
      - a member outcome of `succeeded` means that GENtle completed and
        persisted the child assessment. The child report may still have
        `summary.status = fail|incomplete|not_assessed`; those biological
        verdicts are preserved and summarized as aggregate warnings rather
        than mislabeled as execution errors
      - `collections run restriction-scan` maps the same scanner used by
        `features restriction-scan` over project sequences or explicitly
        sequence-bound gene-set members. Empty enzyme input uses configured
        preferred enzymes and then the built-in default set; unknown names
        fail closed. The returned `gentle.collection_restriction_site_scan.v1`
        wrapper includes successful child reports and aggregate counts because
        restriction scans are not otherwise stored by report id.
      - `collections run construct-reasoning-inspection` assembles a fresh
        construct-reasoning snapshot for each project sequence or explicitly
        sequence-bound gene-set member and returns the recommended inspection
        actions in `gentle.collection_construct_reasoning_inspection.v1`.
        Snapshots are built on an engine clone: `graph_persisted=false`, generic
        `produced_report_ids` stays empty, and the live project's graph store is
        not changed. A failed member remains a typed per-member failure rather
        than aborting successful inspections.
      - `collections run tfbs-scan` maps the same non-mutating scanner used by
        `features tfbs-scan` over project sequences or explicitly
        sequence-bound gene-set members. Motifs are always explicit and use
        the same token expansion and threshold rules as the direct command.
        The returned `gentle.collection_tfbs_hit_scan.v1` wrapper embeds each
        successful child report and labels aggregate counts as retained hits.
        `aggregate_counts_complete=false` identifies failed members, capped
        child scans, or effective motifs that were not scanned; omitted
        `--max-hits` means unlimited per-member computation.
      - `digest` exposes the existing mutating `Digest` operation through the
        shared shell. Enzyme names are explicit and unknown names reject the
        complete request before any fragment is created.
      - `collections run digest` maps the same topology-aware digest over
        project sequences or explicitly sequence-bound gene-set members.
        Preview is the default and returns `gentle.collection_digest.v1`
        without mutation. Apply requires the preview's
        `plan_fingerprint_sha256`; the lock covers source sequence snapshots,
        effective enzymes, output-prefix policy, and sequentially reserved
        fragment ids, so a changed source or namespace rejects before any
        insertion. Successful apply creates individual fragment sequences with
        exact source lineage and does not create a pool/container. `--path`
        writes the portable report JSON only.
      - `collections run export-pool` combines the declared exhaustive members
        of one physical container into the same `gentle.pool.v1` artifact as
        direct `export-pool`. It returns a
        `gentle.collection_pool_export.v1` wrapper with source-container
        provenance, the canonical membership fingerprint, and one common
        report link for all successful members. Non-exclusive containers are
        rejected before the artifact is written because `gentle.pool.v1`
        cannot represent unknown/unlisted physical contents. In-silico
        `selection` containers are also rejected; physical `singleton` and
        `pool` containers remain eligible.
      - the report follows GENtle's computational-artifact contract with
        `op_id`, `run_id`, sequence links, external inputs/database
        fingerprints, request/effective-setting summaries, a reopen hint, and
        supported export kinds
      - independent characterization rows use
        `pass|fail|incomplete|not_run`; design provenance is cited when the
        assessed pair came from a persisted primer-design report, while raw or
        commercial/literature primer strings explicitly remain `not_run`
      - cited design pairs are independently fingerprinted against the
        normalized oligos actually assessed. A drifted pair makes
        `design_provenance` fail rather than trusting only the report id
      - explicit target-space products make intended isoform coverage and
        junction/genomic-carryover interpretation separate four-state
        dimensions; absent transcript products or unresolved genomic geometry
        remain visible as `incomplete` or `not_run`
      - intended isoform coverage compares expected and observed transcript
        products independently of the off-target verdict, so a transcriptome
        specificity failure does not erase successful intended coverage
      - this operation does not silently rerun Primer3/oligo-pair QC, variant
        masking, or repeat/low-complexity masking. Those dimensions remain
        `not_run` unless their own evidence exists
      - uses local BLAST+ through the prepared reference-genome index path; no
        NCBI remote submission is performed
      - saved report mode consumes `REPORT_ID --pair-rank N` from a persisted
        `PrimerDesignReport`
      - explicit mode consumes `--forward SEQ --reverse SEQ`
      - saved-design 5' non-annealing tails are excluded from BLAST queries but
        retained in report provenance
      - output includes BLAST preflight/invocation metadata, per-primer hits,
        forward/reverse products, forward/forward and reverse/reverse warning
        products, intended/unintended classification, and pass/fail summary
      - short-primer searches are exhaustive with respect to database subjects:
        GENtle validates the database sequence count with `blastdbcmd` and sets
        `-max_target_seqs` above that count. `max_hits_per_primer` remains a
        post-search review threshold and does not truncate the search
      - `search_completeness` records the validated database count, required
        subject limit, observed command limit, and `complete|incomplete`
        status. An unproven search can never produce a specificity pass
      - per-hit primer coverage is calculated from that HSP's inclusive
        `qstart..qend` span and primer length; BLAST's subject-aggregated
        `qcovs` value is retained only in raw BLAST provenance. Query bases
        outside a partial HSP count as effective mismatches together with
        aligned mismatches when products are ranked and screened
      - inward-facing products are paired by subject coordinates and strands,
        independently of whether the forward-role primer is the leftmost hit;
        minus-strand intended products are therefore treated symmetrically
      - intended genomic products are matched by prepared-FASTA subject and
        genomic interval, never by equality to a cDNA amplicon length;
        junction-spanning primers may therefore have no contiguous intended
        genomic product
      - transcript-aware assays use a transcript-set intended target. Every
        declared transcript may contribute one expected product, so a clean
        assay targeting many isoforms is not forced into a one-product model
      - Ensembl transcript stable ids match an optional final numeric version
        suffix (`ENST00000465287` matches `ENST00000465287.1`) when recognizing
        and counting intended products. The full observed subject id remains
        in the report. This rule is deliberately restricted to Ensembl
        transcript ids; arbitrary dotted ids and genomic/gene accessions are
        not stripped
      - wrapped BLAST subject identifiers are normalized against the prepared
        FASTA index, repeated warnings are aggregated, and genomic-DNA versus
        transcriptome/cDNA assessments remain separate in the report
      - `--readiness-max-amplicon-bp N` controls the ordinary pass/fail window,
        while `--exploratory-max-amplicon-bp N` retains longer compatible
        products as `long_product_warning` evidence. The legacy
        `--max-target-amplicon-bp N` sets both ceilings
      - whole-panel planning inherits the panel's declared maximum allowed
        amplicon length as its readiness ceiling unless the caller overrides
        it; the report records the effective ceiling source as
        `transcript_assay_panel_allowed_range`, `explicit_override`, or
        `legacy_policy_max`
      - `--report-detail compact` is the default for newly created reports. It
        retains accepted primer sites and intended, terminal-compatible, or
        failing products plus raw/retained counts. `--report-detail full`
        keeps every normalized HSP and candidate product inline. Imported
        handoffs also cite the raw BLAST TSV paths, byte sizes, and digests
      - candidate pairing is indexed by subject, strand, and coordinate; the
        report records `pairing_candidate_comparison_count` so large inputs can
        be regression-tested without elapsed-time thresholds
      - `--path OUTPUT.json` writes the same structured report returned on
        stdout; the persisted copy can later be inspected or exported with
        `primers show-report REPORT_ID` and `primers export-report`
      - use `primers specificity-plan` when a scheduler, ClawBio skill, or
        another wrapper should own the BLAST processes. The returned
        `gentle.primer_specificity_handoff.v1` stores authoritative
        `program`/`args[]`, query FASTAs, expected `-out` files, accepted exit
        codes, an `all_commands_success` completion policy, and `primers
        specificity-import` follow-up metadata; planning itself starts no
        BLAST search process
      - `primers specificity-multi-handoff REQUEST_JSON_OR_@FILE OUTPUT_DIR`
        prepares one pair against 1-8 explicitly selected genomic/transcriptome
        references, with required/optional flags and optional reference-specific
        intended-target mappings. Required indexes must already be validated;
        unavailable optional indexes are unassessed, not zero-hit searches.
        Choose a new output directory with an existing parent. It does not
        execute BLAST, prepare resources or attach panel/order-readiness gates.
      - after both planned commands exit successfully, `primers
        specificity-import HANDOFF.json` parses their TSVs and applies the same
        specificity interpretation used by the inline route. Empty output is
        allowed as a completed no-hit result, so the wrapper must use process
        exit rather than file size as its completion signal
      - regenerating the same deterministic handoff clears its old declared
        output TSVs, preventing a fresh run from importing stale results
      - planned handoffs bind the inspected BLAST database content fingerprint
        and index kind, and each command carries the same proven subject limit;
        import/finalization rejects replacement content at the same prefix.
        Legacy or modified commands lacking sufficient subject coverage yield
        `incomplete`, not `pass`
      - the plan/import convenience commands dispatch the shared
        `PreparePrimerPairSpecificityHandoff` and
        `ImportPrimerPairSpecificityHandoff` operations. The same operation
        JSON is accepted by CLI `op`/`workflow`, MCP `op`, JavaScript
        `apply_operation`, and Lua `apply_operation`
      - for a persisted transcript-assay panel, `primers
        transcript-assay-specificity-plan` emits one aggregate handoff plus a
        process-manifest template covering every selected assay and both primer
        searches. The outer scheduler runs only the declared `program` and
        `args[]`; it must return one row per command with its exit code and the
        exact output byte length/hash, even when a command fails
      - always call `primers transcript-assay-specificity-finalize` after the
        scheduler finishes. A successful command with an empty TSV is a valid
        completed no-hit result. Complete biological evidence can yield
        `specificity_fail`; complete searches whose intended target geometry
        cannot be assessed yield `not_assessed`; missing, duplicate, failed,
        stale, or provenance-mismatched execution evidence yields `incomplete`
      - finalization binds the current panel, assay ids/ranks, annealing
        sequences, policy, prepared genome/BLAST database, handoff schema, and
        output identities. Panels persist their design-time genome anchor so
        later finalization does not depend on mutable live provenance; old
        panels without an anchor remain readable and may resolve to
        `not_assessed`. Only a complete all-assay `pass` is accepted. Complete
        `pass`, `specificity_fail`, and `not_assessed` outcomes atomically
        persist every per-assay report by assay and target kind; an
        `incomplete` outcome is returned/exportable but does not attach a
        partial panel result. The acceptance record separately lists
        biologically failing, not-assessed, incomplete, and execution-failed
        assay ids, so those states and `not_run` remain distinguishable while
        genomic plus transcriptome evidence can coexist
      - after a complete `specificity_fail`, use `primers
        transcript-assay-specificity-redesign REQUEST_JSON_OR_@FILE [--path
        OUTPUT.json]`. The request binds the panel, acceptance bytes, and exact
        original design operation. The read-only report retains passing assays
        and gives each failed required assay either content-bound substitute
        candidates or an explicit constrained-infeasibility disposition;
        substitutes must traverse the same cDNA and genomic gates
      - `primers build-transcript-assay-cdna-similarity-map
        REQUEST_JSON_OR_@FILE [--path OUTPUT.json]` produces the advisory map
        independently when a caller wants to review it before design
      - `primers specificity-alignment-html REPORT_ID OUTPUT.html` renders the
        full-primer alignments already stored in a specificity report and does
        not rerun BLAST

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Primer Design and Handoff](primers.md) | [Next: Transcript Assay Panels](transcript-assays.md)
<!-- manual-page-nav:end -->
