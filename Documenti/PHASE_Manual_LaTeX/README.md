# PHASE 7 manual source

This is the editable LaTeX source of the PHASE 7.2.5 user manual. The original
Overleaf layout, typography, colors, title-page artwork and scientific diagrams
have been retained. `main.tex` is the root document; upload this entire folder
as a ZIP to Overleaf and select `main.tex` as the main file. Compile with pdfLaTeX
and BibTeX.

The obsolete screenshots of the separate PHASE 5/6 `.mlapp` interfaces are no
longer referenced by the document. Those image assets remain in `images/` for
source provenance. Five early PHASE 7 interface screenshots supplied from a
macOS session illustrate the hub and the three sections. A PHASE 7 PSI output
image supplied during the macOS test is included as an operational example,
not as independent scientific validation.

The built PDF is at the repository root as `PHASE_Manual.pdf`.

The 7.2.5 update covers PSI suitability, CSK ground AOIs/DEM normalization, PSI
sub-areas, Radar grid & PS, map zoom, Automatic/Custom models, covariance and
residual controls, query points and current outputs. Historical screenshots are
identified as interface examples. No private test inputs are bundled.

Local build: run `latexmk -pdf -interaction=nonstopmode
-outdir=/tmp/phase-manual-build main.tex` from this folder with pdfLaTeX/BibTeX.
Copy the verified PDF to the repository root. Release assets include that PDF
and a versioned source ZIP for Overleaf.
