# arXiv submission package

**File to upload:** `staple-arxiv.tar.gz` (41 KB)

Contents:
```
staple.tex                      main source (self-contained: inline TikZ/pgfplots
                                figures, inline bibliography — no .bbl, no image files)
anc/README.txt                  describes the ancillary files
anc/tight_staple_proof.py       verifies the gnomon theorem (Thm 5.1, exact signs)
anc/interval_certificate.py     verifies the 116 Sturm certificates of Thm 4.1
anc/interval_construction.json  the 51 rational coefficients of S(u)
anc/arc_certificate.py          constructs + verifies the line-arc witness (Thm 6.1)
```

The `anc/` directory follows arXiv's ancillary-file convention; the files are
listed automatically on the abstract page. The source is pure ASCII and
compiles under pdflatex (a `\pdfoutput=1` hint, guarded so XeTeX ignores it,
is on line 2 for arXiv's autodetection).

## Metadata for the submission form

**Title:**
The three regimes of escape paths for isosceles triangles in Bellman's lost-in-a-forest problem

**Authors:** Alexander Temerev (University of Geneva), Alessio Doria
**(TODO: add Doria's affiliation + email in staple.tex before submitting;**
**there is a marked %% TODO line in the preamble)**

**Abstract:** (plain text, math in $...$)

We organize the best-known escape paths for isosceles triangles $T(\beta)$
(unit legs, base angle $\beta$) in Bellman's lost-in-a-forest problem into
three regimes, and determine the character of both regime boundaries. On
$[45^\circ,60^\circ]$ the Besicovitch--Movshovich zee is optimal;
numerically it remains the frontier down to $\beta_\times=42.287^\circ$.
Below $\beta_\times$ the frontier is a single connected branch of line--arc
paths that traverse the boundary of their convex hull minus a chord:
near-polygonal trapezoidal "staples" at the top, progressively rounded by
arcs whose curvature radius is $\sin\beta$ (the width of the triangle),
terminating on the scaled Zalgaller caliper at $\beta^*\approx27.6^\circ$
--- the "Zalgalloid" family Ward speculated about in 2008. The boundaries
differ in kind. At $\beta_\times$ two combinatorially distinct branches
cross transversally (slopes $-0.016$ and $+0.013$ per degree): each
persists as a local optimum beyond the crossing, the zee resists rounding,
and the upper-bound envelope has a corner there --- its maximum, so the
hardest isosceles forest in the surveyed range sits at the phase boundary,
with $\ell(T)\le1.38920$. At $\beta^*$ nothing crosses: the branch merges
tangentially into the caliper, bar length and caliper gap vanishing
together. Inside the middle regime we prove two anchors. An explicit
piecewise-polynomial family of polygonal staples $S(u)$, $u=\tan(\beta/2)$,
escapes and is strictly shorter than every escaping square path for
$\beta\in[31.08^\circ,42.16^\circ]$, and shorter than the diameter, the
scaled caliper, and the zee on $[31.50^\circ,42.16^\circ]$ --- so Ward's
square path is removed from the frontier everywhere, and the trapezoidal
improvement he reported numerically in 2008 becomes a theorem. The proofs
eliminate the continuum of positions and headings exactly: escape from a
triangle reduces, via the planar case of a containment theorem of Lutwak,
to the containment of a disk in an explicit Minkowski sum, and every
resulting inequality is a rational-polynomial sign condition verified by
exact Sturm certificates. At the type example, the golden gnomon
($\beta=36^\circ$), an exact line--arc certificate gives
$\ell(T)\le1.2826799$, and the critical line--arc path is characterized
algebraically: window tangency forces the curvature radius $\sin\beta$, the
sines of its arc angles are roots of explicit quartics over $\Q(\sqrt5)$,
and the critical length $1.2826760\ldots$ has the same closed shape as
Zalgaller's strip constant. The branch ends above in an exact rectangle
path of length $\sqrt2$ at $\beta=45^\circ$, tied with the diameter.
Optimality is claimed nowhere below $45^\circ$.

**Primary category:** math.MG (Metric Geometry)
**Cross-list:** cs.CG (Computational Geometry)
**MSC classes:** 52A10 (Primary), 52B55, 52A38 (Secondary)

**Comments field:**
19 pages, 11 figures, 2 tables. Ancillary files contain exact verification
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
      `interval_certificate.py` → 116/116 PASS;
      `arc_certificate.py` → ALL WINDOW CERTIFICATES PASS (18 windows)
- [ ] Final proofread of the compiled PDF by the author
- [ ] Doria's affiliation and email added (%% TODO in staple.tex preamble)
