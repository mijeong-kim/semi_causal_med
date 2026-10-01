# Snapshot Verification

The retained numerical outputs were last validated locally on 2026-09-30 using
R 4.6.1 on macOS (arm64). Later document-only checks are recorded below.
These are checks of the supplied JKSS snapshot, not a new full Monte Carlo run.

## Executed Checks

`make smoke validate` completed successfully from the new `github_JKSS` root
on 2026-09-29.
It also passed again after replacing the 13 large CSV files with compressed
copies. No manual extraction was performed for that validation run.

On 2026-09-30, `JKSS_REUSE_RESULTS=1 Rscript run_all.R` passed again after
the supplement reorganization. The main, comparator, paired-fit, calibrated-power,
variance-sensitivity, confounding, full-curve and application checks used the
retained full outputs; no full Monte Carlo rerun was performed.

| Check | Result |
| --- | --- |
| Required R packages | Available in the validation environment |
| Lossless compression, gzip reading and fresh-CSV precedence | Passed |
| All five estimators and both total-effect decompositions | Passed |
| Effect-map derivatives and baseline-mean propagation | Passed |
| JOBS II data, baseline terms and application diagnostics | Passed |
| Huber coefficients, covariance, shift invariance and calibration | Passed |
| Independent calibrated-power output validation | Passed |
| Main simulation and application output validation | Passed |
| Comparator seeds, summaries and representative OLS refits | Passed |
| Variance-misspecification output validation | Passed |
| Confounding-sensitivity algebra, gradients and output | Passed |
| Full-curve population identities, covariance and output | Passed |

The main validation checked 16,000 attempted method fits (15,807 valid).
The comparator validation checked 20,000 attempted method fits (19,940 valid),
all 4,000 data-set seeds and 100 summary cells. The largest discrepancy in the
16 reconstructed comparator OLS fits was 4.89e-15. Variance sensitivity checked
14,000 attempts (13,792 valid); confounding sensitivity checked 6,000 attempts
(5,997 valid), with the full-curve analysis using those same data sets.

Two supplementary Huber checks emitted `quantreg` warnings that an intermediate
quantile-regression solution may be nonunique. All assertions passed; warnings
were not treated as evidence of a failed check.

## Environment

| Package | Loaded version |
| --- | --- |
| Matrix | 1.7.5 |
| rootSolve | 1.8.2.4 |
| sn | 2.1.3 |
| mediation | 4.5.1 |
| medflex | 0.6.11 |
| quantreg | 6.1 |

For this local check, `sn` was loaded from a pre-existing auxiliary R library
through `R_LIBS`; it was not installed in the default library. No library or
machine-specific path is bundled. Other users should install the dependencies
listed in `README.md`, then run `Rscript R/check_dependencies.R`.

## File Checks

- All 39 R files, including the runner and compression helpers, parsed.
- All 42 distinct LaTeX input/figure paths resolved across 33 TeX files.
- The citation file parsed as YAML and identifies an unpublished manuscript.
- The included `jobs` object exactly matched `mediation` version 4.5.1.
- The original snapshot passed a byte-identical comparison of 92 scientific
  source, study-data, numerical-output and manuscript/PDF files. The compact
  revision preserves the scientific inputs and numeric values; 18 CSV-reading
  entry points now use a shared plain/gzip reader. The 2026-09-29 theory update
  leaves the score equations and all reported finite-sample tuning values
  unchanged while expressing numerical safeguards as convergent sequences.
- Every compressed file was restored and checked for MD5 equality before its
  redundant CSV was removed. All 13 uncompressed hashes also match the original
  records in `Latex/JKSS/results`; see `RESULT_COMPRESSION.csv`.
- The revised theory and algorithm pages of both supplied PDFs were rendered
  and visually checked.
- The supplement now follows the manuscript in ten sections. All 16 tables
  and six figures are retained; main-text cross-references match the new
  numbering. Revised simulation, calibration, sensitivity and application pages
  were rendered and checked. LaTeX reported no undefined references or box warnings.
- All 47 CSV/RDS numerical files in the primary manuscript's results directory
  remained byte-identical across the document rebuild.
- Figure 4 was redrawn from the 25 retained JOBS II estimates with a separate
  legend panel and wider legend-column spacing. The standalone figure and its
  manuscript page were rendered and checked; estimates and intervals were unchanged.
- Supplement Section S10 was shortened to a repository link, archive provenance
  and essential reproduction commands. The rebuilt PDF's GitHub hyperlink was
  checked, and the final page was rendered and visually inspected.
- The main-text contribution statements and interpretive prose were revised for
  positive or neutral framing. Theoretical scope, assumptions, numerical results
  and adverse findings were retained; shared supplementary prose was synchronized.
  Both PDFs were recompiled and revised pages were rendered for layout checks.
- Repeated interpretation was consolidated in the main article, with the
  supplement retaining derivations, computational details and complete tables.
  Duplicate variance-sensitivity and full-curve findings paragraphs were removed
  from the supplement. Numbered equations, propositions, numerical files and all
  figure/table assets were preserved; both PDFs and cross-references were checked.
- No cover letters, build logs, local R libraries or temporary checkpoints
  were included. The retained data and PDFs are not excluded by `.gitignore`.

`FILE_MANIFEST.csv` records relative paths, byte sizes and MD5 hashes of the
delivered snapshot, excluding the manifest itself. It is an integrity inventory,
not a package lockfile or cryptographic authenticity claim. Regenerating tables,
figures or validation reports can change hashes; the inventory describes the
current local submission snapshot only. It does not certify a remote GitHub upload.

## Rebuild Scope

The complete Monte Carlo experiments and application analyses were not refitted.
All retained-output validators were rerun, however, and
`JKSS_REUSE_RESULTS=1 Rscript run_all.R` recompiled both PDFs and rebuilt the
reproducibility ZIP from the stored full results. No GitHub upload or release
was performed by this local verification step.

## Repository-only Distribution Update (2026-09-30)

The article and supplement now direct readers to the public GitHub repository
for code, data and retained numerical output; Online Resource 1 remains the
supplementary PDF. References to a separately supplied Online Resource 2 were
removed from the active TeX files and table generators. The runner's existing
ZIP output is described as a local backup, not a submitted supplement.

Both PDFs were rebuilt without final LaTeX warnings. All 12 pages whose text
changed were rendered and visually checked, and the GitHub hyperlinks in both
PDFs were verified. All 47 primary CSV/RDS result files remained byte-identical;
all 35 primary and 39 upload-snapshot R files parsed. Analysis code changed only
in repository wording and the completion message, with no changes to estimators,
simulation settings or stored numerical results. No analyses were refitted and
the full numerical validators were not rerun for this wording-only update.

The final PDFs and this snapshot's file manifest were refreshed. Existing local
ZIP backups were left unchanged. No GitHub upload or release was performed;
the latest files must still be uploaded and checked against the submission.

## Supplementary Material Naming Update (2026-09-30)

At the author's request, the current article's cross-references and the
supplement's title now use "Supplementary Material" without a resource number.
The README and upload instructions use the same name. The submission checklist
and upload guide flag that JKSS's current instructions use "Online Resource"
for attached supplementary files, so the final submission naming needs review.

Both PDFs were recompiled without final LaTeX warnings and retain 32 and 33
pages, respectively. All changed pages were rendered for layout review; the
primary and upload PDFs have identical extracted text. The old resource label
is absent from the active TeX files and both PDFs. The 47 primary CSV/RDS result
files remain byte-identical, and no analysis scripts, results or figures were
changed. The final PDFs and file manifest were refreshed locally; no GitHub
upload, release, simulation rerun or backup ZIP rebuild was performed.

## Equation Numbering Audit (2026-10-01)

Unreferenced equation numbers were removed from nine main-article displays
and one supplementary display. Their mathematical content is unchanged; only
the numbering environments, unused labels and redundant numbering commands
were edited. Normalizing those commands reproduces the previous source text.

All remaining 15 main-article and six supplementary equation numbers have a
direct reference or are included in a cited equation range. In particular,
the kernel-density display keeps its number because it belongs to the cited
range (7)-(10). The audit covered both documents and their included TeX files.
Both PDFs rebuilt without final LaTeX warnings or unresolved references, and
all ten changed pages were rendered and visually checked. Page counts remain
32 and 33. The final PDFs, upload sources and manifest were refreshed locally;
the analysis scripts, study data, numerical results and figure assets were not
changed. No analyses were rerun and no GitHub upload was performed.

## AI Disclosure and Submission Naming (2026-10-01)

The methods section now discloses OpenAI Codex assistance with manuscript
drafting and revision, mathematical derivations and exposition, and development
and debugging of analysis code. The statement records the author's responsibility
without asserting that the author's independent verification has already been
completed. The submission checklist retains that author review as an open item.

This update supersedes the September 30 unnumbered-supplement naming choice.
The article introduces Online Resource 1 (Supplementary Material), uses Online
Resource 1 in subsequent explicit file citations, and provides its description
under Supplementary information. The supplementary title includes the resource
number and identifies the target journal. Sections S1-S10, all 16 tables, all six
supplementary figures and the equation labels are unchanged. Code and numerical
materials remain distributed through GitHub; no Online Resource 2 is submitted.

Both PDFs were rebuilt without unresolved references or box warnings and retain
32 main-article and 33 supplementary pages. The 35 pages with changed extracted
text were rendered and visually checked. Primary and upload PDFs have identical
extracted text; displayed mathematics and source labels were checked against the
pre-edit sources. All 47 retained CSV/RDS numerical files are unchanged.

The rebuild runners now also create `output/pdf/ESM_1.pdf` as a byte-identical
submission copy of `JKSS_supplement.pdf`. Both runners parsed successfully and
their new copy steps were tested separately. Existing supplementary PDF links
are preserved; the two filenames do not represent two different supplements.
The README, upload instructions, final PDFs and file manifest were refreshed.
No simulation, full numerical-validation workflow, backup ZIP rebuild or GitHub
upload was performed for this editorial update.
