# Delegated Execution Approval

Confirmation-gated delegated skills use a two-invocation contract. Natural
language may select and fill a route, but it cannot directly execute that
route.

## 1. Produce A Proposal

Run the normal delegated `gentle.clawbio_skill_request.v1`. The generic runner
performs only deterministic, read-only preflights and returns status
`approval_required`. It writes
`reproducibility/execution_proposal.json` with schema
`gentle.clawbio_execution_proposal.v1`.

The proposal binds:

- the normalized request and exact argument vector;
- the delegated skill, version, descriptor, route, and plan-step hashes;
- the generic runner, catalog, and intent-descriptor file hashes that enforce
  the delegation contract;
- the resolved execution directory and state/output/input paths;
- selected pair/report/material and specificity-space slots;
- caller-supplied biological assumptions;
- the GENtle runtime version and executable-file hashes; known local Cargo
  launchers are reduced to the built `gentle_cli` binary before proposal
  creation, so later execution cannot silently rebuild changed sources;
- material PATH, cache, asset-root, and external-tool environment overrides,
  including executable hashes when an override resolves to a local file;
- project-state and scientific-input hashes before execution;
- a Primer backend resolved to `internal` or a path/version/hash-pinned
  `primer3` executable.

`proposal_digest` is the SHA-256 of the canonical `approval_basis`. The
proposal timestamp and its storage location are deliberately outside that
digest, so repeating the same preflight context produces the same approval
identity even though each proposal file records its own creation time.

No scientific command runs during this invocation.

For OCI execution, use an immutable image reference such as
`ghcr.io/owner/image@sha256:...`; mutable tags are rejected for
confirmation-gated routes. A local Apptainer/Singularity image is bound by its
file hash. Custom launcher scripts can only be bound as far as their visible
command and file arguments allow, so prefer a direct `gentle_cli` executable,
the included local-checkout launcher, a hash-bound local image, or a
digest-addressed OCI image.

## 2. Approve The Exact Digest

The caller control plane must show the proposal to the user and obtain approval
for its exact digest. It then supplies a separate request such as:

```json
{
  "schema": "gentle.clawbio_approved_execution_request.v1",
  "proposal_path": "/absolute/path/execution_proposal.json",
  "approval": {
    "schema": "gentle.clawbio_execution_approval.v1",
    "approval_scope": "execute_exact_proposal",
    "proposal_digest": "sha256:...",
    "approval_id": "caller-issued-id",
    "approved_by": "caller-attested-reviewer",
    "approved_at_utc": "2026-08-03T12:00:00Z"
  }
}
```

GENtle validates the assertion but does not decide who may approve. Identity,
authentication, chat permissions, and the proof that a human actually reviewed
the digest belong to the caller/OpenClaw control plane.

## 3. Execute Without Rerouting

The runner reloads the normalized request from the approved proposal. It never
reconstructs the operation from the original prose. Before the main command it
rechecks the proposal digest, route/descriptor contract, exact command,
working directory, runtime version/files, project state, scientific inputs,
resolved paths, material execution environment, and any pinned Primer3
executable. Any drift fails closed.
For the included local-checkout launcher and Cargo fallback, that exact command
uses the already-built, hash-bound `gentle_cli` executable rather than invoking
Cargo again.

Read-only list/show/preflight routes may remain automatic. A delegated route
that selects biological material, selects a pair/backend, mutates project
state, or writes/replaces an artifact must use this proposal contract.

Direct structured use of `gentle-cloning` without delegated-skill metadata is
unchanged. This gate controls conversational delegation; it does not add a
second scientific implementation or restrict ordinary `gentle_cli` use.

## Multi-Stage Scientific Planning

A deterministic planner and the operations it recommends are two approval
subjects. For gene isoform-assay studies, the first proposal must bind the
complete normalized planner request, including effective defaults, policy and
evidence hashes, prior-plan reference, observations, and override. Planning
then emits exact ordered operation payloads and one batch digest. A second
proposal binds and executes those stored payloads without regeneration. Either
approval authorizes only its exact computation; neither validates biology.

Selecting declared presentation blocks or a named profile from an immutable
canonical GENtle report is not another scientific decision. It may run without
a new scientific approval, while ordinary file-write/overwrite safety and a
content-bound execution receipt remain required.
