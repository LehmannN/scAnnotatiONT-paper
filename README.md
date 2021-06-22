# Improved annotation with scRNA-seq and Oxford Nanopore bulk long reads

Recovering significant 10xGenomics scRNA-seq signal through improved annotation with long reads.
Applied on chick embryo neuro-epithelial progenitors at 66h of development.

GitHub repository started on 2021-04-02.

## Repository content and structure

### Repository root files

Files found at the root of the repository are of general purpose:

* **_bookdown.yml** - Configuration `YAML` file for the notebook (which is writen with [bookdown](https://bookdown.org/))
* **_build.log** - Log file (stdout and stderr) generated when running the analysis
* **_build.sh** - Bash script to run the analysis
* **_deploy.sh** - Bash script to deploy the analysis on GitHub (e.g. `git push`)
* **index.Rmd** - R Markdown file used for the setup of the analysis (load libraries, define variables and paths...)
* **README.md** - This file
* **_workflowr.yml** - [Workflowr](https://jdblischak.github.io/workflowr/index.html) configuration file

### Repository structure

* **analysis** - [R Markdown](https://rmarkdown.rstudio.com/) analysis files (e.g. `Rmd` files)
* **data** - Data used for the analyses
* **docs** - Rendered analysis reports (e.g. `html` files) and figures generated by the notebook (e.g. `png` files)
* **env** - Code used to create the conda environment for this study (e.g. `bash` script) and the corresponding [`YAML` file](https://medium.com/@balance1150/how-to-build-a-conda-environment-through-a-yaml-file-db185acf5d22).
* **figures** - Figures created before the notebook creation and used as input (e.g. `png` or `pdf` files)
* **scripts** - Reusable code (e.g. functions)
* **output** - Output files (e.g. `tsv` or `csv` files)

### Analysis

Current analysis files include:

* **01-Impact-ref-annotation-scRNA.Rmd** - Here we explore the discrepancies between the references annotations (Ensembl and RefSeq) and their impact on common scRNA-seq analyses.
* **02-Incomplete-annotations-induce-signal-loss.Rmd** - We study here the loss in scRNA-seq signal (e.g. genes) due to significant deficiencies in the reference annotation, specially in the 3' UTRs annotations.
* **03-Approaches-to-improve-transcriptome-with-Long-Reads.Rmd** - We compare various tools dedicated to transcriptome reconstruction in bulk RNA-seq (StringTie2, scallop), a dedicated signal detection approach and broad 3' UTR extension and apply them to scRNA-seq data and ONT bulk long reads.
* **04-Impact-reannotation-scRNA.Rmd** - We assess the impact of our various reannotations on common scRNA-seq analyses.
* **05-Validation-of-novel-genes-with-scRNA.Rmd** - We evaluate the ability of our approach to identify novel genes and use scRNA-seq analyses as a filter to highlight genes of biological interest in chick embryo neuro-epithelial progenitors.
* **06-A-tool-and-pipeline-to-improve-annotation-for-scRNA.Rmd** - Description and recommendations to use our pipeline on other scRNA-seq data.
* **07-Session-info.Rmd** - Session info output.

### Data

Input data files include:

*TBA*


### Scripts

Current code files include:

*TBA*

### Output

Output files (all the files in this folder are created when running the analysis).
