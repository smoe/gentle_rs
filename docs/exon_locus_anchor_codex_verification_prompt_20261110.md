# Codex prompt: current exon-locus anchor verification

Use this prompt with Codex after Steffen's latest push. Resolve and report the
exact `origin/main` SHA before running anything; do not infer that a nearby
checkout or this fixture branch is the requested source.

> Glen, bitte nach Steffens Push den aktuellen `main` mit Anchor-Korrektur
> `5aed0ae9` gemäß
> [aktuellem Handoff](https://github.com/smoe/gentle_rs/blob/main/docs/glen_next_verification_handoff.md)
> gezielt prüfen:
>
> - **Anchor-Regression:** Query `[65,77)`, Referenz `[165,177)`, Anchor **171**,
>   weiterhin außerhalb des Fensters ab 200; JSON/TSV vollständig, kein
>   ungültiges Glyph.
> - **Offene Exon-Locus-Prüfungen:** aktuelle Locked-/Default-/Headless-Checks
>   und noch fehlende Core-/Service-Regressions; nichtleere Testauswahl
>   sicherstellen.
> - **Echte Ausführung:** JS-/Lua-VMs, MCP sowie quellenidentischer
>   CLI-Rehearsal mit LF/CRLF, allen fünf Exportformaten, deterministischen
>   Ausgaben, Pixel-Guard und No-clobber. Skips zählen nicht als Abnahme.
> - Native GUI nur bei verfügbarer Umgebung; Windows über passende normale CI,
>   kein zusätzlicher Dispatch.
>
> Die abgeschlossenen Renderer-, `7fa6b15a`-, Reader-/BLAST-/Barcode-Kampagnen
> nicht wiederholen. Claudes **32/32** bleiben ein korrigierter
> Scratch-Overlay-Lauf; die geschätzten 360 Reads sind keine gemessene
> Kapazitätsgrenze.
>
> Eigenes Target und freien Speicher vorab prüfen. Bitte SHA,
> Binary-/Lock-/Input-Hashes, Plattform/Toolchain, Befehle, Testzahlen, Skips
> und manifestierte Evidenz zurückgeben. Keine Quelländerungen, privaten
> Datenläufe, Pushes oder Veröffentlichungen.

The real TP53/TP73 cDNA fixture delivered on the accompanying Glen branch is a
separate, owner-authorized repository asset. Its presence does not expand this
verification prompt into a private-data analysis or authorize Codex to modify
source, commit, push, publish, or repeat already closed campaigns.
