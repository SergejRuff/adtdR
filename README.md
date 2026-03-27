# R implementation of ADTD

This repository provides an **R implementation of Adaptive Digital Tissue Deconvolution (ADTD)** based on the Python implementation in **Deconomix**:

- Deconomix: https://gitlab.gwdg.de/MedBioinf/MedicalDataScience/Deconomix
- Deconomix bioRxiv preprint: https://www.biorxiv.org/content/10.1101/2024.11.28.625894v1

The implementation here aims to reproduce the core ADTD workflow in R while following the structure and behavior of the original **Deconomix/ADTD** codebase as closely as possible.

## Background

ADTD is a digital tissue deconvolution approach that estimates:

- **cellular compositions**
- **cell-type-specific reference adaptations / gene regulation**
- **hidden background contributions**

The original ADTD method was described in:

**Görtler, F., Mensching-Buhr, M., Skaar, Ø., Schrod, S., Sterr, T., Schäfer, A., ... & Altenbuchinger, M. (2024). _Adaptive digital tissue deconvolution_. Bioinformatics, 40(Supplement_1), i100-i109.**

The broader Python framework used as the basis for this R port was described in:

**Mensching-Buhr, M., Sterr, T., Seifert, N., Völkl, D., Tauschke, J., Rayford, A., ... & Altenbuchinger, M. (2024). _Deconvolution of omics data in python with deconomix–cellular compositions, cell-type specific gene regulation, and background contributions_. bioRxiv, 2024-11.**

## Purpose of this repository

This repository is intended for users who would like to run ADTD-style analyses in **R** instead of Python.  
It includes R implementations of the main algorithmic components inspired by the original ADTD implementation.

## Notes

- This is an **R reimplementation**, not the original Deconomix/ADTD package.
- The original Python implementation and methodology remain the primary reference.
- For methodological details, please consult the original publications listed above.
- There is no gurantee given that this package gets new updates or will be maintained in the future.
- The original authors of the python implementation had no involvement in the R implementation. If you find bugs/issues please contact me (package maintainer) and not the original authors of ADTD.
- some functions such as prop_metrics and sim_bulk_seurat_sce are new implementations not found in the python original.

## References

1. Görtler, F., Mensching-Buhr, M., Skaar, Ø., Schrod, S., Sterr, T., Schäfer, A., ... & Altenbuchinger, M. (2024). *Adaptive digital tissue deconvolution*. **Bioinformatics**, 40(Supplement_1), i100-i109.
2. Mensching-Buhr, M., Sterr, T., Seifert, N., Völkl, D., Tauschke, J., Rayford, A., ... & Altenbuchinger, M. (2024). *Deconvolution of omics data in python with deconomix–cellular compositions, cell-type specific gene regulation, and background contributions*. **bioRxiv**.