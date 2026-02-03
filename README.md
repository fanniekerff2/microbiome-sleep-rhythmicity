# microbiome-sleep-rhythmicity

This is the code repository associated with the manuscript [The gut microbiota and sleep in infants: a focus on diurnal rhythmicity patterns](https://www.biorxiv.org/content/10.1101/2025.07.11.664337v1).

Please note that the publication includes results from 20 infants (163 samples). However, because the parents of 3 infants declined to make their data public, the publicly available dataset associated with this repository (**PRJEB104111**) contains only the samples from the remaining 17 infants (130 samples).

## Setup

To reproduce the analysis, after cloning the repository (`git clone https://github.com/fanniekerff2/microbiome-sleep-rhythmicity.git`), we strongly suggest you set up a conda environment with the required dependencies.

### Setup the project environment: `qiime2-amplicon-2024.10`

Start by installing [mamba](https://github.com/mamba-org/mamba) in the base environment (if not installed yet). Then install the project environment (`qiime2-amplicon-2024.10`), either using the `environment.yml` file in the directory (option 1) or via the official QIIME 2 distribution file and plugin repositories (option 2).

```shell
conda install mamba -n base -c conda-forge
```

#### Option 1: installation using the `environment.yml` file

Install the project environment via the file.

```shell
mamba env create -f environment.yml  # if on Apple Silicon chip following [1] prefix with: CONDA_SUBDIR=osx-64
conda activate qiime2-amplicon-2024.10
# if on Apple Silicon chip following [1] also run: conda config --env --set subdir osx-64
# qiime info # can be run to check for a correct installation

pip install https://github.com/caporaso-lab/q2-boots/archive/refs/tags/2024.10.beta.zip
# qiime boots --help # can be run to check for a correct installation

pip install q2_kmerizer@git+https://github.com/bokulich-lab/q2-kmerizer.git@main
# qiime kmerizer --help # can be run to check for a correct installation
```
[1] Mac users with Apple Silicon chips (M1, M2, etc) must configure the installation of QIIME 2 plugins in Rosetta 2 emulation mode (as ARM builds are not yet available). For more information on this, see [QIIME2 guidelines](https://docs.qiime2.org/2024.10/install/native/).

#### Option 2: installation using the official QIIME 2 distribution file and plugin repositories

Using the correct installation instructions required for your machine ([QIIME2 guidelines](https://docs.qiime2.org/2024.10/install/native/)), create a fresh installation of QIIME2 Amplicon version 2024.10. For example, below are the instructions for a macOS machine with an Apple Silicon chip.

```shell
CONDA_SUBDIR=osx-64 conda env create \
  --name qiime2-amplicon-2024.10 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2024.10/amplicon/released/qiime2-amplicon-macos-latest-conda.yml
conda activate qiime2-amplicon-2024.10
conda config --env --set subdir osx-64
# qiime info # can be run to check for a correct installation
```

Then, install the required dependencies for this project:

```shell
conda install -c conda-forge scikit-learn seaborn joypy statsmodels
pip install openpyxl

pip install q2_kmerizer@git+https://github.com/bokulich-lab/q2-kmerizer.git@main
# qiime kmerizer --help # can be run to check for a correct installation

mamba install -y \ # if on Apple Silicon chip following [1] prefix with: CONDA_SUBDIR=osx-64
   -c https://packages.qiime2.org/qiime2/2024.10/metagenome/released/ \
   -c conda-forge -c bioconda -c defaults \
   q2-fondue
# qiime fondue --help # can be run to check for a correct installation

pip install https://github.com/caporaso-lab/q2-boots/archive/refs/tags/2024.10.beta.zip
# qiime boots --help # can be run to check for a correct installation
```

#### Mandatory configuration for option 1 and 2

* Refresh the QIIME 2 CLI cache:
```shell
qiime dev refresh-cache
```
* Run the `vdb-config` tool and exit by pressing x (needed to initialize the wrapped SRA Toolkit for more information see [here](https://github.com/ncbi/sra-tools/wiki/05.-Toolkit-Configuration))
```shell
vdb-config -i
```
* In case you need to configure a proxy server (for e.g. on ETH HPC), run the following command:
```shell
vdb-config --proxy <your proxy URL> --proxy-disable no
```

### Setup the additional environment: `empress`

To prevent `empress`'s older dependencies from breaking your main QIIME 2 installation, create `empress` as a standalone environment.

```shell
CONDA_SUBDIR=osx-64 conda create -n empress python=3.9
conda activate empress
# if on Apple Silicon chip following [1] also run: conda config --env --set subdir osx-64
pip install cython "numpy>=1.12.0,<2.0" # or pip install cython "numpy >= 1.12.0" # required for empress installation
pip install iow --no-build-isolation
pip install empress
# empress --help # can be run to check for a correct installation
```

### Setup the additional environment: `qiime2-amplicon-2025.7`

Using the correct installation instructions required for your machine ([QIIME2 guidelines](https://docs.qiime2.org/2024.10/install/native/)), create a fresh installation of QIIME2 Amplicon version 2025.7. This more recent version of QIIME 2 is required to run [ANCOM-BC2](https://www.nature.com/articles/s41592-023-02092-7) within QIIME 2, a compositionally-aware linear regression model that allows testing for differentially abundant features across sample groups while also implementing bias correction. For example, below are the instructions for a macOS machine with an Apple Silicon chip.

```shell
CONDA_SUBDIR=osx-64 conda env create \
  --name qiime2-amplicon-2025.7 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2025.7/amplicon/released/qiime2-amplicon-macos-latest-conda.yml
conda activate qiime2-amplicon-2025.7
conda config --env --set subdir osx-64
# qiime info # can be run to check for a correct installation
```

Then please run the mandatory configuration step above.

---

## Repos structure

```
microbiome-sleep-rhythmicity/
  ├── scripts/            <- python notebooks/scripts
  │   └── qiime2/         <- sbatch scripts
  ├── data
  │   ├── meta_data/      <- processed metadata files
  │   │   └──raw/         <- raw metadata files
  │   ├── seq_data/       <- raw sequences
  │   └── processed_data/ <- processed sequence data
  ├── figures/            <- vizualizations (.pdf)
  ├── environment.yaml
  ├── README.md
  └── .gitignore
```

---

## Data processing and analysis overview

This project is implemented through a series of **Python notebooks**, supported by modular **Bash scripts** for specialized tasks. All Bash scripts called by the notebooks are located in the `scripts/qiime2/` directory, while core utility functions are maintained in `functions_script.py` to ensure clean, reproducible code.

The project workflow is structured into the following stages:

### 1. Data pre-processing
This phase involves the end-to-end preparation of both sequencing and metadata.
* **Sequencing Data:** Includes fetching raw sequences, denoising, filtering, taxonomic annotation, and the computation of diversity measures (including k-mers).
* **Metadata Processing:** Involves integrating alpha diversity metrics, generating beta diversity tables, and calculating specialized longitudinal metrics such as babySQUID scores, microbial rhythmicity, and temporal volatility per age group.

*Associated Notebooks:* `1-seq-proc.ipynb`, `1.2-meta-proc.ipynb`

### 2. Exploratory data analysis (EDA)
In this stage, we examine the dataset's foundational characteristics. This includes evaluating sample distributions across timepoints and visualizing community structure using **Principal Coordinates Analysis (PCoA)** based on beta-diversity metrics.

*Associated Notebook:* `2-EDA.ipynb`

### 3. Primary analysis
The core analysis investigates the temporal dynamics of the gut microbiome and its relationship with infant sleep development:

* **Diurnal patterns:** Investigation of rhythmicity patterns within the samples, specifically focusing on the synchronization between gut microbiota and infant sleep rhythmicity patterns.
    * *Associated Notebook:* `3-rhythmicity.ipynb`
* **Microbial volatility & environmental factors:** Analysis of microbial temporal volatility across different infant ages and the subsequent effects of feeding history and sleep patterns on the gut microbiota diversity.
    * *Associated Notebook:* `4-factors.ipynb`
* **Microbiota-gut-sleep axis modulation:** Exploration of how the microbiota influences infant sleep quality (babySQUID) and sleep rhythmicity (CFI). This includes a **Random Forest model** utilizing microbial composition to predict babySQUID scores.
    * *Associated Notebooks:* `5-mb-sleep-CFI.ipynb`, `5.2-mb-sleep-babySQUID.ipynb`, `5.3-mb-sleep-RF.ipynb`
* **Gut melatonin exploration:** A targeted exploratory analysis investigating the relationship between gut melatonin levels, gut microbiota composition and sleep outcomes.
    * *Associated Notebook:* `6-melatonin.ipynb`

---

## Contact and license

If you use this project code and/or data, please [cite us](https://www.biorxiv.org/content/10.1101/2025.07.11.664337v1) accordingly. In case of questions or comments, feel free to raise an issue in this repository or to contact the lead author via e-mail.

---
