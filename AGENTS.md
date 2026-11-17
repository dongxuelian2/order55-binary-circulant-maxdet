# Repository instructions

## Mandatory legacy research audit

- Review **every existing tracked research document**, not only the files that contain Chinese text. Audit definitions, quantifiers, hypotheses, edge cases, lemma dependencies, completeness of proofs, calculations, code-to-claim consistency, bibliographic accuracy, and all statements of novelty or scope. Resolve ambiguities and supply necessary missing definitions, arguments, examples, and references.
- Translate Chinese-language mathematics into accurate English without changing its meaning, preserving equations, labels, hypotheses, and all proof steps. If a purported proof has an actual gap, repair it rigorously if possible; otherwise identify the precise unresolved step and correctly classify the result as incomplete. Never silently upgrade evidence or invent conclusions.
- Organize corrections into focused per-repository PRs. Review the full diff for residual Chinese prose, broken mathematical links and citations, and unintended proof changes. Run existing relevant checks and provide proof-level review; automated language scans and tests alone are insufficient.
- Keep all future research files unambiguous and in English, including comments and manuscript prose. Avoid duplicative progress reports or unnecessary engineering scaffolding.

This is a research repository. Keep proof statements, mathematical explanations, and program documentation precise, verifiable, and easy to navigate. Preserve the established directory layout and proof assumptions when editing.

## Mandatory English-language policy

- English is the sole language for repository-authored explanatory prose. Write and maintain `README.md`, `AGENTS.md`, `PROOF.md`, research notes, Markdown documentation, manuscript exposition, source-code comments/docstrings, and pull-request descriptions in English. Use English for new or revised public-facing messages and instructions.
- Audit **all tracked Markdown files**, including nested directories and archives, for Chinese-language explanatory text. Translate existing Chinese proofs and documentation into precise, idiomatic mathematical English in dedicated PRs. Preserve definitions, hypotheses, quantifiers, equations, proof dependencies, numbering, cross-references, bibliographic details, code behavior, and the exact strength of every claim. Never silently shorten or reinterpret a proof during translation.
- Keep mathematical notation, source identifiers, and faithful quotations or official titles in their original form when required; provide an English explanation around any unavoidable original-language quotation. Do not leave mixed-language paragraphs, headings, tables, or code comments.
- Before completing a documentation PR, inspect every changed file for residual Chinese prose and broken links/references. Review mathematical equivalence against the original and run existing relevant checks. A character scan is a useful audit aid, not a replacement for reading or proof review.
- This requirement applies to all future contributions and edits, not only the present translation cleanup. Do not create new Chinese prose or alternate Chinese documentation copies. Keep English clear, direct, and technically specific; avoid defensive boilerplate and unsupported claims.
