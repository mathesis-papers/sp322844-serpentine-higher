# Serpentine matrices: the two conjectures in OEIS A322844

Revision 2, November 5, 2026.

**Haobo Ma**, ChronoAI Pte Ltd, Singapore; **Wenlin Zhang**, National University of Singapore; ChronoAI Pte Ltd, Singapore.

License: Creative Commons Attribution 4.0 International (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/.

## Results and sources

The note extends Stefano Spezia's determinant result recorded on December 1, 2018 in [OEIS A322277](https://oeis.org/A322277). It proves both assertions in the conjectures block of [OEIS A322844, Comments, revision 68](https://oeis.org/A322844). The second-trace sequence dates to December 28, 2018; the higher-trace assertion was added on January 6, 2019, in revision 5. Both appear together under the conjectures heading in revision 14, January 31, 2019. The [entry history](https://oeis.org/history?seq=A322844) supplies those locators.

For the one-based matrix whose odd rows increase and even rows decrease, the rank is two for every size `n >= 2`. Its characteristic polynomial is `x^(n-2)(x^2-T_n*x-a_n)`, where

```
T_n = n*(n^2 + n%2)/2
a_n = n^2*(3*(1+n^2)-2*(2+n^2)*(n%2))/12
```

The second exterior trace is `-a_n`; every exterior power of order greater than two is the zero map. The note also proves the spectrum, minimal polynomial, ordinary-power recurrence, a coefficient formula for the permanent, row/column/diagonal identities, diagonal aggregates, full and parity generating functions, and seven reciprocal sums across A322844, A002415 and A037270. It rederives the related `i+j` and `i-j` matrix formulas recorded by Cloitre and Somos in 2002. Explicit corrections distinguish the conventional exponential series from shifted-factorial weights and correct particular formulas in the linked entries. Existing determinant, matrix and sum formulas retain their attribution; no priority claim is made.

## Files and building

- `note.tex` contains the complete argument and an embedded bibliography.
- `verify.py` is a standard-library Python 3 checker; `checks.json` is its deterministic output.
- `metadata.json` records the byline, license, problem identifiers and mathematical scope.
- `manifest.json` records exact file hashes in the source archive and the accompanying PDF hash.

Build with **Tectonic 0.17.0**:

```sh
tectonic note.tex
```

Alternatively, use a LaTeX distribution containing AMS packages, `geometry`, `booktabs`, `tabularx`, `longtable`, `lmodern` and `hyperref`:

```sh
pdflatex note.tex
pdflatex note.tex
```

The bibliography is embedded; BibTeX is unnecessary. Run the checks from the extracted source directory without Python optimization:

```sh
python3 verify.py > checks-reproduced.json
cmp checks.json checks-reproduced.json
```

## Mathematical coverage

All assertions in the note have written proofs. Exact finite checks cover every principal minor through size 10 (2,047 minors), full characteristic polynomials through size 16, factorization/rank/sums through size 32, permanents through size 9, the ordinary-power recurrence through size 16, both linked matrix families through size 16, generating coefficients through degree 64, and reciprocal partial fractions through index 64. Numerical reciprocal sums use 4,096 terms and explicit tail bounds; they supplement the analytic proof.

The pinned Lean sources use the same matrix, with indices shifted by one:

- [SpeziaSerpentineSecondExteriorTrace.lean](https://github.com/the-omega-institute/trureturing/blob/d9bb8a406783c5222849b798f3c2ab738f653d75/D5/S3/Quantum/Matrix/SpeziaSerpentineSecondExteriorTrace.lean): `result` expresses the absolute second characteristic-coefficient exterior trace for `n > 1`; `coefficient_two` expresses the signed coefficient; `charpoly_factor` factors the characteristic polynomial through a two-by-two product for `n >= 2`.
- [SpeziaSerpentineHigherExteriorTrace.lean](https://github.com/the-omega-institute/trureturing/blob/d9bb8a406783c5222849b798f3c2ab738f653d75/D5/S3/Quantum/Matrix/SpeziaSerpentineHigherExteriorTrace.lean): `result` expresses the vanishing characteristic-coefficient exterior trace for every natural size and every order `k > 2`, including `k > n`.

The note proves the principal-minor bridge to actual exterior-space traces. The zero-map assertion, spectral/minimal-polynomial results, permanent formula, scalar generating functions, reciprocal sums and linked sum identities have their written proofs in the note; they are distinct from those Lean declarations. The coverage table identifies source assertions, earlier recorded results, corrections and consequences not posed in the sources. Unrelated sequence interpretations and bivariate series outside the triangular matrix-sum domain are not covered.
