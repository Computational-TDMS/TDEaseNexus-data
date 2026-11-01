# TDEaseNexus soybean study

This workspace contains the soybean Analysis-mode session and its data.

- `session.html`: recorded session for viewing in a browser.
- `.tdease/`: session and tool-call records.
- `soybean/`: three mzML files, TopPIC/TopFD outputs for ten samples, the FASTA database, and analysis results under `analysis/`.
- `zhou2022-source-search/`: published-study TopPIC results and settings for comparison.
- `toppicr-comparison/`: data for the TopPICR comparison in Figure S2 and Table S4.

The original dataset is available from MassIVE under accession MSV000089984.

## Data sources

Soybean data: Zhou et al., *Discovery top-down proteomics in symbiotic soybean root nodules*, Frontiers in Analytical Science (2022), [doi:10.3389/frans.2022.1012707](https://doi.org/10.3389/frans.2022.1012707). Original data: [MassIVE MSV000089984](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=6ac410bc49ec41f1a91f49e063943cce). The FASTA contains UniProt sequences and contaminant entries; UniProt data attribution: [UniProt Consortium](https://www.uniprot.org/).

The mzML files are stored with Git LFS. Install Git LFS before cloning to obtain the data files.
