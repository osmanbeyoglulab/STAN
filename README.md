# STAN, a computational framework for inferring spatially informed transcription factor activity across cellular contexts

STAN (Spatially informed Transcription factor Activity Network) is a linear mixed-effects computational method that predicts spot-specific, spatially informed TF activities by integrating curated TF-target gene priors, mRNA expression, spatial coordinates, and morphological features from corresponding imaging data. 

## Acquiring datasets:
- Visium Spatial Gene Expression datasets: [10x Genomics lymph node](https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-1-0), [10x Genomics glioblastoma](https://www.10xgenomics.com/datasets/human-glioblastoma-whole-transcriptome-analysis-1-standard-1-2-0), [Ravi et al. glioblastoma](https://doi.org/10.5061/dryad.h70rxwdmj), [Wu et al. breast cancer](https://zenodo.org/record/4739739)
- Visium CytAssist Spatial Gene and Protein Expression datasets: [Tonsil](https://www.10xgenomics.com/datasets/gene-protein-expression-library-of-human-tonsil-cytassist-ffpe-2-standard), [Tensil (with addons)](https://www.10xgenomics.com/datasets/visium-cytassist-gene-and-protein-expression-library-of-human-tonsil-with-add-on-antibodies-h-e-6-5-mm-ffpe-2-standard), [Breast cancer](https://www.10xgenomics.com/datasets/fresh-frozen-visium-on-cytassist-human-breast-cancer-probe-based-whole-transcriptome-profiling-2-standard)
- scRNA-seq datasets: [Wu et al. breast cancer](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078)

## Methods for inferring TF activity
- `notebook01_run_stan.ipynb`: the proposed method, including reading datasets and inference pipeline
- `notebook02_run_ridge.ipynb`: the Ridge regression model
- `notebook03_run_decoupler.ipynb`: the decoupleR method

## Analysis of STAN
- `notebook04_run_cv.ipynb` and `notebook04_anaysis_cv.ipynb`: to analyze the impact of the parameters in the regression
- `notebook05_run_image.ipynb` and `notebook05_analysis_image.ipynb`: to analyze the impact of including morphological features
 
## Applications of STAN
- Identifying germinal centers and cell-type-specific TFs in the human lymph node: `analysis_lymphnode.ipynb`
- Identifying similar/different TFs associated with pathological regions in breast cancer: `analysis_breast_part1.ipynb` (a series notebooks)
- Linking ligands and receptors to TFs in glioblastoma: `analysis_glioblastoma.ipynb`
- Linking surface proteins to TFs in CytAssist datasets: `analysis_cytassist.ipynb`

![image](https://github.com/osmanbeyoglulab/Tutorials-on-ISMB-2024/blob/main/hands-on_tutorial/session-1/resources_stan/stan.png?raw=true)

## Supporting resources
A gene set resource comprising TF–target gene priors from [hTFtarget](https://guolab.wchscu.cn/hTFtarget/#!/) will be obtained, and TFs identified in the Human Transcription Factor database [(humantfs)](https://www.cell.com/cell/fulltext/S0092-8674(18)30106-5) will be retained to generate the TF–target gene prior matrix. 

## Citation
L Zhang, A Sagan, B Qin4, H Wang, E Kim, B Hu, HU Osmanbeyoglu (2024). STAN, a computational framework for inferring spatially informed transcription factor activity across cellular contexts. _bioRxiv_, bioRxiv (2024): 2024.06.26.600782. 

doi: <https://doi.org/10.1101/2024.06.26.600782>
