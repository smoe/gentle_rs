# TSS Collection Outer-Agent Experiment

Date: 2026-09-30

Reviewed source: `f70be167c3f60e673987ef7a5fde517cdf6c88a3`

Fixture: public synthetic tutorial 08.15 starter

The planner smoke used ClawBio source
`629a50350e6781ed47378a5de1f586e55edc09f4` with an isolated temporary skill
root; unrelated untracked files in the long-lived ClawBio checkout were not
read or changed.

## Question

What does the outer ClawBio/OpenClaw layer add beyond direct GENtle use, and
can its tutorial-relevant mapping be verified without running GENtle's inner
agent?

## Baseline

The merged source already ships a generic `gentle-cloning` runner and a focused
PCR descriptor. Its 59 generic intents contained no TSS, transcription-start,
or tutorial 08.15--08.17 route. The generic runner could nevertheless execute
an exact TSS shared-shell command supplied by a caller.

That distinction is intentional:

- GENtle owns command parsing, annotation geometry, preview identity,
  derivation, validation, state, typed results and biological warnings.
- The generic wrapper owns exact execution plus a topic-neutral receipt.
- A focused descriptor owns natural-language discovery, route identity and the
  decision to require proposal-then-approved execution.

## Measurement

Two runs started from the same starter bytes
`sha256:efd4626f662d5b752b973196bd497c5fc99c3999668dbeed2f530f3a33191b67`.

### Read-only inventory

Direct `gentle_cli` and the generic wrapper ran:

```text
promoters tss-inventory @docs/examples/assets/tss_tutorial/inventory.request.json
```

Both emitted exactly the same stdout bytes:

```text
sha256:ad0f436184ad26854f7856426158f2c667aa490fd675dcf97b506a0b47a85f99
```

Both left project bytes unchanged. The wrapper additionally retained the
normalized request, input-file hash, exact command, exact runtime binary hash,
before/after state hashes, exit status, native output hashes, report and replay
environment. It did not change or reinterpret GENtle's result.

### Generic structured materialization

A direct structured wrapper request ran:

```text
promoters tss-materialize @docs/examples/assets/tss_tutorial/materialize.request.json
```

It immediately created the three selected windows, emitted stdout
`sha256:ff078b14858a853d77b083c8024051c355d8635608c8746b5692bd594d95abc5`,
and changed the project hash. It emitted no execution proposal. This is the
documented direct-use boundary, not a bypass: direct structured GENtle calls do
not acquire conversational authority rules merely because they use the wrapper.

### Focused delegated materialization

The new `gentle-tss-collection` descriptor selected the same command and
delegated it with route identity `materialize_reviewed_starts`. The first call:

- returned `approval_required`;
- ran only runtime preflight, not the scientific command;
- retained the identical before/after project hash;
- bound the descriptor, catalog, route, plan step, resolved request file,
  executable, project state and exact command;
- produced proposal digest
  `sha256:96031715acbdadf84d7dbd315fec403c455974b9a110d753e9b7e741ffd17605`.

A second synthetic-test invocation attested that exact digest and executed the
stored proposal. It created the same three IDs and emitted the same materialized
result bytes as the generic run. The two resulting project files differed only
in operation timestamps (`created_at_unix_ms` and `recorded_at_unix_ms`), not in
the scientific payload.

## Result

The outer layer adds three concrete things:

1. **Discovery:** it maps bounded natural-language terms to one reviewed shared
   command instead of requiring the caller to know GENtle's grammar.
2. **Authority separation:** mutating conversational routes become a bound
   proposal and require approval of an exact digest; direct structured API use
   remains available and is not misrepresented as conversational approval.
3. **Portable evidence:** every run binds request, inputs, runtime, state,
   outcome and replay material without duplicating biological computation.

It does **not** add a TSS algorithm, improve the scientific result, prove that a
promoter is active, authenticate the approver, or establish that a GUI host
opened a window.

## Inner-Agent Conclusion

The inner agent is not required to verify the tutorial-relevant mapping covered
here. The outer test can verify:

- every route template against tutorial 08.15's machine-readable parity cases;
- shared-shell parsing through the same GENtle command grammar;
- read-only versus mutating classification;
- descriptor/catalog/delegate compatibility;
- proposal non-execution, state binding and exact approved replay.

The inner agent remains useful for model-response quality, conversational voice,
context selection and live GUI-host adoption. Those are different acceptance
questions and were not claimed by this experiment.

## Bounded Scope

This first descriptor covers tutorial 08.15's six cases. Tutorials 08.16 and
08.17 contribute another 13 outer-agent cases for report attachment, local
factor scoring, SVG export and public TP73 handling. They should be added only
as further descriptor routes over existing shared commands, preserving their
artifact-write and biological-selection approval requirements. No new wrapper
mode is needed.

Not tested here: a live natural-language router, inner model, native GUI host,
Windows/macOS, private data, real-locus biological acceptance, or approver
authentication.
