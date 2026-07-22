Ancillary files for "From the equilateral triangle to the flat limit: the
zee-to-caliper frontier in Bellman's forest problem" (A. Temerev, A. Doria).

  tight_staple_proof.py       Verifies the golden-gnomon polygonal staple:
                              reconstructs the Minkowski sum W in the field
                              Q(sqrt5, sqrt(10-2*sqrt5)), determines the exact
                              algebraic signs of the twelve edge certificates
                              of Table 1, and checks the length comparisons.
                              Expected output: "ALL CERTIFICATES PASS".
                              Runtime: ~2 minutes.

  interval_construction.json  The 51 rational coefficients of the piecewise-
                              polynomial staple family S(u) of the interval
                              theorem (Theorem 2.1), together with the fixed
                              support-test vertex choices used in the escape
                              certificates.

  interval_certificate.py     Verifies all 116 Sturm certificates of the
                              interval theorem (Theorem 2.1 and Corollary
                              2.2): re-derives every polynomial from the
                              construction data and the exact triangle
                              geometry, then checks strict positivity on each
                              rational interval by endpoint signs and exact
                              root counts.
                              Expected output: "116/116 certificates PASS".
                              Runtime: ~10-15 minutes.

  arc_certificate.py          Constructs the rational line-arc witness of the
                              exact line-arc anchor (Proposition 4.1: bar + tangent
                              segments + circular arcs with fully rational
                              data) and verifies its 18 direction-window
                              certificates and the length enclosure
                              ell(T) <= 1.2826798604... < 1.282680.
                              Expected output: ALL WINDOW CERTIFICATES PASS.
                              Runtime: ~3 minutes.

Requirements: Python 3.9+, sympy, numpy, scipy.
Code and data are also maintained at
https://github.com/atemerev/bellman-staple
