# data/meta

- The pulsar key is `name` (as in glitches.csv: B names where they exist, otherwise J names).
  `jname` has no leading "J". The 21 rows are the pulsars with >= 5 usable glitches.
- source_audit.csv is prefilled with the 49 primary references (and how many of the 21 pulsars'
  glitches each supplies). Several "New - Our measurement" variants are the catalogue's own
  measurements (Basu et al. style); merge them when you audit. The top 8 references cover ~80% of glitches.
- Every value needs a source. When two end dates are defensible, use the earlier one.

## Provenance of source_audit.csv

- The 49 reference rows and their `n_glitches` counts were already present when this project's files were first
  committed (2026-10-31), with every other column empty. The records do not show how that list and those counts
  were produced, and no count of catalogue entries reproduces the counts (see the caption of Table 8 in the paper,
  which describes `n_glitches` as a bookkeeping count whose derivation was not recorded).
- The audit columns of 12 rows (the eight largest sources, and the catalogue maintainers' own measurements, which
  appear under four spellings) were filled on 2026-10-31 from an audit table supplied in the author's instructions
  for that step. The records do not show how that table was compiled. The other 37 rows were not audited.
