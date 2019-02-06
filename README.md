For the purpose of the webinar we will be working on the head node of an HPC cluster running Slurm. This is the most likely infrastructure that fellow bioinformaticians already find themselves using on a regular basis.

The execution of the Snakemake workflow will actually take place on the cluster head node with jobs being submitted to Slurm for queing and processing. From the head node, Snakemake will monitor the submitted jobs for their completion status and submit new jobs as dependent jobs complete sucessfully.

# Snakemake Installation

See also: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

The recommended installation route for Snakemake is through a conda environemnt. As such, you need Anaconda3, usually avaiable to you on your cluster via the module system.

```bash
# Load the Anaconda3 module on your cluster
# If it's unavailable contact the cluster sysadmin
module load \
  Anaconda3

# Create a new conda environment for snakemake v5.4.0 (currently the latest version)
# https://anaconda.org/search?q=snakemake
conda create \
  --name snakemake \
  --channel bioconda --channel conda-forge \
  --yes \
  snakemake=5.4.0
```

# Snakemake Usage

First, load the conda environment you set up:

```bash
# Load the Anaconda3 module on your cluster
# If it's unavailable contact the cluster sysadmin
module load \
  Anaconda3

# activate the conda environment you created above
source activate snakemake
```

# Example Workflow

```bash
# Clone this repo
git clone https://github.com/UofABioinformaticsHub/2019_EMBL-ABR_Snakemake_webinar

cd 2019_EMBL-ABR_Snakemake_webinar/snakemake-tutorial
```

Now you are ready to run Snakemake. Lets visalise the DAG of job dependencies:

```bash
snakemake \
  --dag \
  | dot -Tsvg \
  > dag.svg
```

See what rules and commands would be executed for a couple of differet targets:

```bash
# The "setup_data" target can be used to create a small reference and read data sets for testing purposes
snakemake \
  setup_data \
  --dryrun --printshellcmds

# Run the whole pipeline (target "all"), which would also generate the test data set
snakemake \
  --dryrun --printshellcmds
```

Provided in this repo is a slurm profile, cluster config files and a conda environment file. This allows Snakemake to interact with Slurm and have the relevant software available to each rule.

Running the pipeline, for real, can be accomplished using:

```bash
snakemake \
  --profile slurm \
  --cluster-config cluster-configs/default.yaml \
  --cluster-config cluster-configs/phoenix.yaml \
  --use-conda
```

Add a new accession and take a look at the DAG of job dependencies

```bash
# Uncomment the "Alsen" line
sed -i 's/^#  "(Alsen)",/  "\1",/' Snakefile

# Generate DAG of job dependencies
snakemake \
  --dag \
  | dot -Tsvg \
  > dag2.svg
```

Add in the remaining accessions:

```bash
# Uncomment all the accessions
sed -i -r 's/^#  "(.+)",/  "\1",/' Snakefile

# Generate DAG of job dependencies
snakemake \
  --dag \
  | dot -Tsvg \
  > dag3.svg
```

Wow, that's a big/messy DAG. Lets look at the simpler DAG of rules:

```bash
snakemake \
  --rulegraph \
  | dot -Tsvg \
  > rulegraph.svg
```

Lets just have a look at a count of jobs that would be executed:

```bash
# Get a count of all jobs that would be run
snakemake \
  --dryrun \
  --quiet
```

# Cleanup the Enviroment

Once we have finished, we will deactive the conda environment:

```bash
source deactivate
``` 
