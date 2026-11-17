# Repository instructions

## Zero-defensive-prose rule — mandatory

Apply this rule to all repository-authored English prose: README, mathematical research notes, `PROOF.md`, LaTeX exposition, documentation, pull-request descriptions, and code comments. Enforce it when auditing historical material and when writing anything new.

**Prohibited rhetorical constructions (including grammatical variants):**
- Chinese examples to recognize and eliminate in translated source material: `不是……`, `并不是……`, `不构成……`, `还没有做……`, `尚未……`, `而不是……`, `并未……`, `不意味着……`, `不代表……`, `不能据此……`, `这并不证明……`, `我们不声称……`, `不保证……`, `仅仅是……`.
- English equivalents to eliminate from authored explanations: `not a ...`, `is not ...`, `does not ...`, `do not ...`, `has not ...`, `have not ...`, `not yet ...`, `rather than ...`, `instead of ...`, `does not constitute ...`, `does not establish ...`, `does not imply ...`, `should not be interpreted as ...`, `we do not claim ...`, `no claim is made ...`, `merely ...`, `only a preliminary ...`, `beyond the scope ...`, `not intended to ...`.
- Delete disclaimers, self-exculpation, premature qualification, repetitive status caveats, and defensive contrasts. Do not relocate them to a different paragraph or replace them with synonyms serving the same rhetorical purpose. Zero such constructions are acceptable in explanatory prose.

**Required replacement style:** Lead with the actual theorem, deliverable, tested behavior, evidence, applicability conditions, exact scope, and concrete open mathematical question, using affirmative statements. Examples:
- `We prove the classification for n >= N under hypotheses H1 and H2.`
- `Verified cases: n = 2,...,7. General-case proof target: ...`
- `CI checks: 16 passed. Proof: Proposition 4. Research priority: external review.`
- `Algorithm: FPT for disjoint ordinary cardinality quotas. Next target: robust mixed-status source-SCC accessibility.`

**Mathematical truth takes precedence over stylistic cleanup.** Preserve essential logical negation in theorems, inequalities, counterexamples, constraints, source quotations, and executable code where changing it would alter correctness. In ordinary explanation, express boundaries by positive quantified scope and precise next objectives. If a claim is unsupported, correct or remove it; never invent evidence to avoid a negative sentence.

Before opening or updating a PR, review every added or modified paragraph for prohibited rhetoric and its semantic equivalents. Record the direct result and evidence; keep mathematical limitations in compact, exact statements of domain and open problems.

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
