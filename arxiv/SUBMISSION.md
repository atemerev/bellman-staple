# arXiv submission package

**File to upload:** `staple-arxiv.tar.gz` (24 KB)

Contents:
```
staple.tex                      main source (self-contained: inline TikZ/pgfplots
                                figures, inline bibliography — no .bbl, no image files)
anc/README.txt                  describes the ancillary files
anc/tight_staple_proof.py       verifies Theorem 3 (exact algebraic signs)
anc/interval_certificate.py     verifies the 116 Sturm certificates of Theorem 1
anc/interval_construction.json  the 51 rational coefficients of S(u)
```

The `anc/` directory follows arXiv's ancillary-file convention; the files are
listed automatically on the abstract page. The source is pure ASCII and
compiles under pdflatex (a `\pdfoutput=1` hint, guarded so XeTeX ignores it,
is on line 2 for arXiv's autodetection).

## Metadata for the submission form

**Title:**
An improved upper bound for Bellman's lost-in-a-forest problem in isosceles triangles

**Authors:** Alexander Temerev (University of Geneva)

**Abstract:** (plain text, math in $...$)

We prove a new upper bound for the escape length of isosceles triangles in
Bellman's lost-in-a-forest problem, on an entire interval of base angles
inside the unsolved range. We construct an explicit piecewise-polynomial
family of three-segment paths $S(u)$, $u=\tan(\beta/2)$, with
isosceles-trapezoidal convex hull ("staples"), and prove that $S(u)$ escapes
the triangle $T(\beta)$ with unit legs and base angle $\beta$, and is
strictly shorter than every escaping square path, for all
$\beta\in[31.08^\circ,42.16^\circ]$; on $[31.50^\circ,42.16^\circ]$ it is
also strictly shorter than the diameter, the scaled Zalgaller caliper, and
the Besicovitch--Movshovich zee. This interval strictly contains
$(32.36^\circ,41.34^\circ)$, the range on which Ward's square path was the
shortest known candidate, so the square path is removed from the frontier
everywhere, and the trapezoidal improvement Ward reported numerically in
2008 becomes a theorem, with the enlarged range he predicted certified at
both ends. The proofs eliminate the continuum of positions and headings
exactly: escape from a triangle reduces, via the planar case of a
containment theorem of Lutwak, to the containment of a disk in an explicit
Minkowski sum, and every resulting inequality is a rational-polynomial sign
condition verified by exact Sturm certificates. At the type example, the
golden gnomon (apex $108^\circ$), the bound reads
$\ell(T)\le(3751558+2\sqrt{20693661817348})/10^7=1.2849616\ldots$ against
Ward's $1.2934740\ldots$, with the closed-form bracket
$\frac{3}{10}\sqrt{25-5\sqrt5}\le\ell(T)$. Optimality is not claimed.

**Primary category:** math.MG (Metric Geometry)
**Cross-list:** cs.CG (Computational Geometry)
**MSC classes:** 52A10 (Primary), 52B55, 52A38 (Secondary)

**Comments field:**
12 pages, 6 figures, 2 tables. Ancillary files contain exact verification
code (Python/SymPy); code and data also at
https://github.com/atemerev/bellman-staple

**License:** author's choice; CC BY 4.0 recommended if you want the widest
reuse, otherwise the arXiv non-exclusive license is the minimal option.

## Pre-submission checklist

- [x] Source compiles from the exact tarball contents in a clean directory
- [x] Pure ASCII source; no external image files; no .bbl needed
- [x] `\pdfoutput=1` hint present (XeTeX-guarded)
- [x] Ancillary scripts re-run from the packaged copies:
      `tight_staple_proof.py` → ALL CERTIFICATES PASS;
      `interval_certificate.py` → 116/116 PASS
- [ ] Final proofread of the compiled PDF by the author
