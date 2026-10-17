# Reactive oxygen species inhibit the differentiation of monocyte-derived dendritic cells through a PPARγ-dependent mechanism

Analysis code for the RNA-seq data in Figures 3C–E and 4A–D, 4F–J of:

> Mazzitelli IG, *et al.* Reactive oxygen species inhibit the differentiation of monocyte-derived dendritic cells through a PPARγ-dependent mechanism. *Cell Death Differ* (2026). DOI: [to be added]

Human monocytes from 4 donors were cultured with GM-CSF, GM-CSF+IL-4 or GM-CSF+NAC (plus freshly isolated monocytes), and analysed by bulk RNA-seq.

## Data

Raw sequencing data are available at GEO under accession **[GSEXXXXXX]**.

Reads were pseudoaligned with kallisto (version [X.X.X]) against [reference transcriptome/index]. The script starts from the kallisto abundance files.

## How to run

1. Clone this repository and open `paper-nac.Rproj` in RStudio (or set the working directory to the repository root).
2. Place the input files as follows:
   ```
   data/
   ├── kallisto/                    16 kallisto abundance files (<donor>_<condition>.tsv,
   │                                e.g. D1_GM.tsv, D1_GM+IL4.tsv, D1_GM+NAC.tsv)
   └── transcripts_to_genes.csv     transcript-to-gene table for tximport
   ```
3. Run:
   ```r
   source("run_all_figures.R")
   ```

The script writes one PDF per panel to `figures/` and the intermediate tables (DESeq2 results, TF activity, GSEA, ORA) to `results/`. It requires internet access to download the CollecTRI network (decoupleR).

## Panels generated

| Panel | Content | Method |
|---|---|---|
| 3C, 3D | Cytokines and chemokines, GM-CSF+NAC vs GM-CSF+IL-4 | DESeq2 |
| 3E | Cytokines and chemokines, GM-CSF+NAC vs GM-CSF | DESeq2 |
| 4A | PCA | Donor-wise z-scores of normalized counts |
| 4B | Top 350 variable genes | Donor-wise z-scores, hierarchical clustering |
| 4C | Volcano plot, GM-CSF+NAC vs GM-CSF+IL-4 | DESeq2 |
| 4D | Transcription factor activity | decoupleR (ULM, CollecTRI), paired t-test |
| 4F, 4H | Hallmark GSEA | fgsea, paired t-test ranking |
| 4G | Hallmark over-representation | clusterProfiler |
| 4I | NRF2-associated genes | Donor-matched log2FC |
| 4J | Cholesterol/sterol genes | DESeq2 |

Panels not listed were generated outside R (Western blot, flow cytometry, qPCR, ELISA).

## Software

R and package versions used are recorded in `results/sessionInfo.txt` after running the script. Main packages: tximport, DESeq2, decoupleR, fgsea, msigdbr, clusterProfiler, ComplexHeatmap, ggplot2.

## Contact

Ignacio Gabriel Mazzitelli — INBIRS (UBA–CONICET), Buenos Aires, Argentina
