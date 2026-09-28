# Relational Clock Transport in Connected Jacobi Systems

**Complete Observables and Clock-Fiber Geometry — Scott A. Jordan**  
ORCID: [0009-0005-7444-3386](https://orcid.org/0009-0005-7444-3386)

Author manuscript from the Reports on Mathematical Physics submission package supplied on 27 September 2026. This description identifies the package, not journal acceptance or publication.

- [Manuscript PDF](manuscript.pdf)
- [LaTeX source](manuscript.tex), [BibTeX references](references.bib), and [compiled bibliography](manuscript.bbl)
- [Online Resource 1: reproducibility code](Online_Resource_1.py)
- Figures: [Fig. 1 PDF](Fig1.pdf), [Fig. 1 EPS](Fig1.eps), [Fig. 2 PDF](Fig2.pdf), [Fig. 2 EPS](Fig2.eps)
- [Dependencies](requirements.txt) and [numerical verification record supplied with the package](numerical_check.txt)

The manuscript studies compact phase-clock transport in connected Jacobi systems, including the independent transport factors for multiple clocks, clock-fiber geometry, stationary-action sensitivity, and induced clock inertias from a stable fast-mode model.

## Reproduce the numerical example

With Python 3.12, run from this directory:

```sh
python -m pip install -r requirements.txt
python Online_Resource_1.py
```

The script overwrites both PDF and EPS figures in this directory and prints the conserved charges, endpoint transport ratio, and sensitivity comparison. Its calculations reproduce the illustrative two-clock numerical example; they do not numerically verify every analytical result in the revised paper. All inputs are provided in the script and manuscript, with no external dataset needed.

The requirements file pins NumPy 1.26.4, SciPy 1.13.1, and Matplotlib 3.9.0. The supplied numerical verification record describes a separate run with newer versions. Both files are preserved as supplied. The script is unchanged in content from the previously tested GitHub supplement, apart from line endings, and includes the portable optional `pdfinfo` lookup.

## Compile the manuscript

This revision uses the standard LaTeX `article` class and the packages listed in its preamble. Keep the bibliography and PDF figures beside the source. With pdfLaTeX and BibTeX installed:

```sh
pdflatex manuscript.tex
bibtex manuscript
pdflatex manuscript.tex
pdflatex manuscript.tex
```

Alternatively, run `latexmk -pdf manuscript.tex`. The earlier Springer template files are not required. The uploaded PDF is the separately supplied `Jordan_Relational_Clock_Transport_ROMP.pdf`; its extracted text matches the package PDF on all 20 pages. The PDF was not rebuilt during this upload.

## Version history

This revision replaces the earlier manuscript titled *Relational Clock Transport Across Interaction Change: A Connected Jacobi Model* in the main collection. The [14 September 2026 version and supporting files](https://github.com/sjordannc/progression-research-papers/tree/7c2759b316db6e0314e22edb13082dc4d31c83a0/relational-clock-transport) remain available at their permanent link.

The editorial cover letter is not part of the public collection. See the [main collection](../README.md) for related papers and repository rights information.
