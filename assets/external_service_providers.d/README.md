# Suppliers Publishing Readable Protocols

Scope revised by the owner on 2026-10-19: **global suppliers with public,
machine-readable protocol documents**, not European companies only. Readable
HTML and text-bearing PDFs qualify; the format is explicit. This first
source-backed directory is **non-exhaustive**, not a guarantee that every
eligible supplier has been found. Product pages, order forms, promotional
brochures or a country office alone do not qualify as protocol evidence.

## Runtime And Filtering

`protocol_suppliers.json` extends the existing provider catalog. All twelve
new records are `catalog_only`: no GENtle quote adapter, direct submission,
status tracking or automatic selection is implemented. The base GeneArt and
metabion records retain their existing handoff contracts and are not overridden.
They need not satisfy this directory's protocol criterion to remain inspectable.

`country_code` is the ISO 3166-1 alpha-2 country of the **named company or
operation**, with `country_source_url` for review. For international suppliers,
this is not every country where they operate, a statement of ultimate ownership,
manufacturing origin, shipping eligibility or commercial preference. Unknown
legacy country metadata is absent, not a guessed country. The doctor validates
the uppercase two-letter shape and source URL, not corporate domicile.

Each `protocols[]` reference carries a title, official URL, format
(`html` or `text_pdf`), `checked_on` date and scope notes. These are readable
documents, **not JSON/XML workflow definitions or executable instructions**.
No protocol is downloaded, parsed or executed when listing the catalog.
A dated reference is not a byte-pinned copy; review the current document and
retain its revision/hash separately before relying on it in a real experiment.
Older manuals can remain readable after a product has been discontinued.

## Official Sources (Checked 2026-10-19)

| Supplier / ID | Country | Role | Published protocol |
| --- | --- | --- | --- |
| Addgene / `addgene` | [US](https://www.addgene.org/about/) | Nonprofit plasmid repository, not a commercial manufacturer | [Plasmid cloning by restriction enzyme digest](https://www.addgene.org/protocols/subcloning/) (`html`); [product/service context](https://www.addgene.org/) |
| Agilent Technologies / `agilent` | [US](https://www.agilent.com/en/contact-us/page) | Manufacturer and product brands | [pBluescript II XR predigested vector manual](https://www.agilent.com/Library/usermanuals/Public/212240-12_pBluescriptIIXR_Manual.pdf) (`text_pdf`); [product/service context](https://www.agilent.com/en/product/mutagenesis-cloning/cloning-vectors-kits) |
| Applied Biological Materials (abm) / `applied_biological_materials` | [CA](https://www.abmgood.com/contact) | Manufacturer and construct-service supplier | [Pro Ligation-Free Cloning Kit quick-start protocol](https://www.abmgood.com/document/p/r/pro_ligation-free_cloning_kit_manual_-_20210326.pdf) (`text_pdf`); [product/service context](https://www.abmgood.com/pshuttle-vector-a003.html) |
| GenScript / `genscript` | [US](https://www.genscript.com/contact.html?src=footer) | International custom synthesis and cloning supplier; named US operation | [GenBuilder cloning kit protocol](https://www.genscript.com/gsfiles/techfiles/cloning_kit_protocol_GenBuilder.pdf) (`text_pdf`); [product/service context](https://www.genscript.com/gene_synthesis.html) |
| IBA Lifesciences / `iba_lifesciences` | [DE](https://www.iba-lifesciences.com/downloadcenter/) | Manufacturer | [StarGate direct transfer cloning protocol](https://www.iba-lifesciences.com/media/84/b0/92/1723034008/Protocol_StarGate-Direct-Transfer-Cloning.pdf?ts=1733731878) (`text_pdf`); [product/service context](https://www.iba-lifesciences.com/applications/cloning/cloning/) |
| Integrated DNA Technologies / `idt` | [US](https://www.eu.idtdna.com/page/contact-us) | International synthesis supplier; named US operation | [Blunt-end cloning protocol and strategies](https://www.eu.idtdna.com/page/support-and-education/decoded-plus/cloning-strategies-part-3-blunt-end-cloning/) (`html`); [product/service context](https://www.idtdna.com/pages/products/genes-and-gene-fragments/custom-gene-synthesis) |
| Merck / Sigma-Aldrich / Novagen / `merck_novagen` | [DE](https://digitalhub.merckgroup.com/legal/imprint) | Manufacturer and product brands | [pET system manual TB055](https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/product/documents/147/299/tb055.pdf) (`text_pdf`); [product/service context](https://www.sigmaaldrich.com/US/en/product/mm/70777?emdredirect=1) |
| New England Biolabs / `neb` | [US](https://www.neb.com/en-us/about-neb) | Manufacturer | [Restriction endonuclease reaction protocol](https://www.neb.com/en-nz/protocols/2012/12/07/optimizing-restriction-endonuclease-reactions) (`html`); [product/service context](https://www.neb.com/products/n3041-puc19-vector) |
| NZYtech / `nzytech` | [PT](https://www.nzytech.com/media/dds/brochurescertificates/mb053_ifu_en_v2301.pdf) | Manufacturer | [NZY-A PCR cloning instructions for use](https://www.nzytech.com/media/dds/brochurescertificates/mb053_ifu_en_v2301.pdf) (`text_pdf`); [product/service context](https://www.nzytech.com/media/dds/brochurescertificates/mb282_mb319-mb330_ifu_en_v2401.pdf) |
| OriGene Technologies / `origene` | [US](https://www.origene.com/about-us/contact-us) | Manufacturer | [Entry vector cloning protocol](https://www.origene.com/support/learning-resources/protocols/entry-vector-cloning-protocol) (`html`); [product/service context](https://www.origene.com/products/vectors/) |
| Promega / `promega` | [US](https://www.promega.com/support/contact-promega/) | Manufacturer; requested worked example | [pGL4 luciferase reporter vectors technical manual TM259](https://worldwide.promega.com/-/media/files/resources/protocols/technical-manuals/0/pgl4-luciferase-reporter-vectors-protocol.pdf?rev=395cb76a0eed4807ae0bbd6f67b7bb77&sc_lang=en) (`text_pdf`); [product/service context](https://www.promega.com/products/luciferase-assays/genetic-reporter-vectors-and-cell-lines/promoterless-firefly-luciferase-basic-vectors/?catNum=E6651) |
| Twist Bioscience / `twist_bioscience` | [US](https://www.twistbioscience.com/company/contact-us) | Custom DNA synthesis supplier | [Gene fragment cloning guidelines](https://www.twistbioscience.com/content/dam/twistbioscience/resources/2024-07/FRM-001257_CloningGuidelines_REV4%20singles.pdf) (`text_pdf`); [product/service context](https://www.twistbioscience.com/products/genes/resources) |

Promega is the owner's requested worked example, not a ranking preference.
The existing helper entry for pGL4.10[luc2] identifies E6651, AY738222.1,
4,242 bp and circular topology; the commercial sequence is not bundled.
The technical manual covers multiple vectors, which are not interchangeable.

## Provenance And Maintenance

Origin: hand-curated factual references to official supplier publications.
No vendor sequence, artwork, reaction instructions, pricing or advertising
performance claims are copied. Recreate by reviewing the named official
protocol and company/country reference, confirming readable document text, then
updating the JSON review date and scope notes. Keep stable IDs and distinguish
manufacturers, service providers and nonprofit distributors. Do not infer that
all product lines have protocols or that all protocols remain current.

Add any further qualifying supplier regardless of country, using these fields
and an actual protocol reference; do not populate the directory from reseller
lists or invent a protocol URL. A supplier without verified protocol evidence
is unassessed for this directory, not biologically or commercially unsuitable.
API-documentation flags describe GENtle's audited adapter contract, not whether
a supplier offers an API.

Use: shared provider loader/doctor, preflight and readiness guards,
catalog/round-trip tests and [tutorial 09.04](../../docs/tutorial/09-04_cloning_time_cost_and_protocol_suppliers.md).
Local listing is offline. Opening an official website is a separate explicit
network action; no private sequences, credentials, quotes or orders are sent.
