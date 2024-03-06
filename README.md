# distiller-sm

## A modular Hi-C mapping pipeline for reproducible data analysis.

The `distiller` pipeline aims to provide the following functionality:

- Align the sequences of Hi-C molecules to the reference genome
- Parse .sam alignment and form files with Hi-C pairs
- Filter PCR duplicates
- Aggregate pairs into binned matrices of Hi-C interactions

### Installation

First, clone the repository:

`git clone https://github.com/open2c/distiller-sm.git`

The recommended way to get all the requirements is to create a conda environment using `workflow/envs/environment.yml`.
We recommend using mamba to handle creation and modification of conda environments, like this:
```
cd distiller-sm
mamba env create -f workflow/envs/environment.yml
```
Other than snakemake, it just installs sra-tools, needed to run the test.
Feel free to simply install snakemake in any way you like.

To check your installation, first run the workflow with a small test dataset.
```
conda activate distiller-sm
sh setup_test.sh
snakemake --use-conda --cores $Ncores --configfile config/config.yml
```
This will also create all required conda environments which will be reused in future runs of the same workflow.

### Customization

To setup a new project, modify the file `config/config.yml`.

Other than your .fastq files, you'll need the sequence of your reference genome of choice and the chrom.sizes file. You need to specify the path to the basename of the bwa index (typially same as the path to the genome sequence, potentially without the extension), but `distiller-sm` can generate the index for you if it doesn't exist in the provided location.

To start the workflow, run `snakemake --use-conda --cores $Ncores --configfile config/config.yml`.

#### Resources

To modify resources used by specific jobs, you should edit the file workflow/profiles.default/config.yml, or write your own workflow profile and point to its folder in the above command with --profile <path> argument.

See more information here: https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles