# arXiv submission package

**File to upload:** `staple-arxiv.tar.gz`

Contents:
```
staple.tex                      main source (inline bibliography; two PNG figures)
results/fig_regimes.png         Figure 1: the zee-to-caliper frontier
results/fig_zalgalloid_morph.png Figure 2: the Zalgalloid continuation
anc/README.txt                  describes the ancillary files
anc/tight_staple_proof.py       verifies the golden-gnomon polygonal certificate
anc/interval_certificate.py     verifies the 116 Sturm certificates of Thm 2.1
anc/interval_construction.json  the 51 rational coefficients of S(u)
anc/arc_certificate.py          constructs + verifies the line-arc anchor (Prop 4.1)
```

The `anc/` directory follows arXiv's ancillary-file convention. The source
compiles under pdflatex (a `\pdfoutput=1` hint, guarded so XeTeX ignores it,
is on line 2 for arXiv's autodetection).

## Metadata for the submission form

**Title:**
From the equilateral triangle to the flat limit: the zee-to-caliper frontier in Bellman's forest problem

**Authors:** Alexander Temerev (University of Geneva), Alessio Doria
**(TODO: add Doria's affiliation + email in staple.tex before submitting;**
**there is a marked %% TODO line in the preamble)**

**Abstract:** (plain text, math in $...$)

In Bellman's lost-in-a-forest problem (Bellman, 1956), a hiker lost in a
forest of known shape -- but at an unknown position and heading -- seeks the
shortest path that is guaranteed to reach the boundary.  For a convex forest
this is the shortest escape path: a path no congruent copy of which fits in
the forest's interior.  We study the isosceles triangles $T_\alpha$ with unit
legs and apex angle $60^\circ\le\alpha<180^\circ$.  At the equilateral end
the Besicovitch--Movshovich zee is optimal through $\alpha=90^\circ$; at the
flat end the natural construction is Zalgaller's strip caliper.  We map the
unresolved upper-bound frontier between them.

Our exact result removes Ward's proposed square phase.  We construct an
explicit piecewise-polynomial family of trapezoidal staples which escapes
$T_\alpha$ and is shorter than every escaping square path for
$95.6734^\circ\le\alpha\le117.8562^\circ$.  On
$95.6734^\circ\le\alpha\le117.0062^\circ$ it also beats Ward's other three
comparators: the zee, the diameter, and the scaled caliper.  Escape is reduced
to containment of a disk in one Minkowski sum, and the resulting $116$
polynomial inequalities are certified exactly by Sturm's theorem.

Numerical continuation supplies the global synthesis.  The best paths found
after the zee lie on one line--arc Zalgalloid branch: it crosses the zee
transversally at $\alpha_\times\approx95.43^\circ$ and meets the caliper
tangentially near $\alpha^*\approx124.8^\circ$, thereafter being the caliper.
Thus the current candidate upper envelope from the equilateral triangle to the
flat limit has a single nonsmooth transition, rather than distinct square,
trapezoid, Zalgalloid, and caliper phases.  We also give an exact line--arc
anchor at $\alpha=108^\circ$.  Optimality beyond $\alpha=90^\circ$ remains
open.

**Primary category:** math.MG (Metric Geometry)
**Cross-list:** cs.CG (Computational Geometry)
**MSC classes:** 52A10 (Primary), 52B55, 52A38 (Secondary)

**Comments field:**
7 pages, 2 figures. Ancillary files contain exact verification code
(Python/SymPy); code and data also at
https://github.com/atemerev/bellman-staple

**License:** author's choice; CC BY 4.0 recommended if you want the widest
reuse, otherwise the arXiv non-exclusive license is the minimal option.

## Pre-submission checklist

- [x] Source compiles from the exact tarball contents in a clean directory
- [x] Ancillary scripts re-run from the packaged copies:
      `tight_staple_proof.py` -> ALL CERTIFICATES PASS;
      `interval_certificate.py` -> 116/116 PASS;
      `arc_certificate.py` -> ALL WINDOW CERTIFICATES PASS (18 windows)
- [ ] Final proofread of the compiled PDF by the author
- [ ] Doria's affiliation and email added (%% TODO in staple.tex preamble)
