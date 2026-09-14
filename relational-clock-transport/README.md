# Relational Clock Transport Across Interaction Change

**A Connected Jacobi Model — Scott A. Jordan**  
ORCID: [0009-0005-7444-3386](https://orcid.org/0009-0005-7444-3386)

Manuscript and supporting materials prepared for submission to *Foundations of Physics*, from the author's submission package dated 14 September 2026. This is an author manuscript, not a journal acceptance or published version of record.

- [Read the manuscript](manuscript.pdf)
- [LaTeX source](manuscript.tex) and [BibTeX references](references.bib)
- [Online Resource 1: reproducibility code](Online_Resource_1.py)
- Figures: [Fig. 1 PDF](Fig1.pdf), [Fig. 1 EPS](Fig1.eps), [Fig. 2 PDF](Fig2.pdf), [Fig. 2 EPS](Fig2.eps)
- [Verified numerical output](numerical_check.txt)

## Reproduce the figures and numerical results

From this directory, using Python 3.12:

```sh
python -m pip install -r requirements.txt
python Online_Resource_1.py
```

The script regenerates `Fig1.pdf`, `Fig1.eps`, `Fig2.pdf`, and `Fig2.eps` in its own directory, overwriting those figure files, and prints the conserved charges, endpoint transport ratio, and fixed-endpoint sensitivity check. All model inputs are in the script and manuscript; no external dataset is required. Poppler's `pdfinfo` is optional.

The GitHub script differs from the submission-package script only in its optional executable lookup: it uses Python's `shutil.which` so the check works on Windows as well as Unix-like systems. The numerical model and calculations are unchanged. Use this GitHub copy if you also want that portability fix in the journal's Online Resource 1 upload.

## Numerical verification

Tested with Python 3.12, NumPy 1.26.4, SciPy 1.13.1, and Matplotlib 3.9.0. The script regenerated both PDF and EPS figures and reproduced the manuscript's conserved charges and endpoint ratio. The finite-difference sensitivity was `0.3538851790` and the direct integral was `0.3538851789`, differing by about `1.84e-10`. These agree with the manuscript at its quoted approximate precision; last digits may vary with numerical libraries.

The verification environment required a test-only temporary-directory workaround for Windows sandbox permissions during EPS font export. This workaround is not part of the distributed script. The figure files committed here are the originals supplied in the submission package.

## Compile the manuscript

Use a LaTeX installation with pdfLaTeX and BibTeX. Obtain `sn-jnl.cls` and `sn-mathphys-num.bst` from the [official Springer Nature journal template](https://www.springernature.com/gp/authors/campaigns/latex-author-support), and place them beside `manuscript.tex`. The submission package uses the December 2024 template. These publisher template dependencies are not redistributed here.

```sh
pdflatex manuscript.tex
bibtex manuscript
pdflatex manuscript.tex
pdflatex manuscript.tex
```

Keep `references.bib` and both PDF figures in the same directory. The supplied manuscript PDF is the author's compiled submission copy; it was not rebuilt for this repository update.

The editorial cover letter remains in the author's private submission package. See the [main collection](../README.md) for other manuscripts and repository rights information.
