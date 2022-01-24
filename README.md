# TiRP

T cells acquire a regulatory phenotype when their T cell receptors (TCRs) experience an intermediate-high affinity interaction with a self-peptide presented on MHC. Using TCR sequences from FACS-sorted human cells, we identified TCR features that shape affinity to these self-peptide-MHC complexes, finding that 1) CDR3β hydrophobicity and 2) certain TRBV genes promote Treg fate. We developed a scoring system for **TCR-intrinsic regulatory potential (TiRP)** and found that within the tumor microenvironment clones exhibiting Treg-Tconv plasticity had higher TiRP than expanded clones maintaining the Tconv phenotype. To elucidate drivers of these predictive TCR features, we examined the two elements of the Treg TCR ligand separately: the self-peptide via murine Tregs, and the human MHC II molecule via human memory Tconvs. These analyses revealed that CDR3β hydrophobicity promotes reactivity to self-peptides, while TRBV gene usage shapes the TCR’s general propensity for MHC II-restricted activation.

A demonstration of TiRP calculation is provided in **TiRP_demo.ipynb.**

To see how TCR features are defined, check out **define_TCRfeatures_demo.ipynb.**

Jupyter notebooks to reproduce all main figures are provided in **display_items.**. For clarity and brevity, functions are defined separately in **utils.R.**