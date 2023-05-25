# Nextflow De Novo Assembly Pipeline for Oxford Nanopore Data

This pipeline automates the process of generating a de novo assembly of long-read sequencing data produced by Oxford Nanopore Technology.

## Technical Considerations

### Minimum Read Length
The assembler used for this pipeline, [Flye](https://github.com/fenderglass/Flye), has a minimum read length limit of 1000bp. For most data generated using Nanopore, this should be fine. However, if the majority of your reads are less than this length, this pipeline will not be useful for your analysis. In my experience, the noisy/long read assembler [Miniasm](https://github.com/lh3/miniasm) allows parameters to be set for reads shorter than 1000bp (Find a helpful tutorial for this tool [here](https://faculty.washington.edu/sr320/?p=13602)).

### Medaka Model

Medaka requires information about the pore, sequencing device, and basecaller. This information is specified to the tool through a 'model', which is a string of text in the following format:
```
{pore}_{device}_{caller variant}_{caller version}
```

The pipeline requires a medaka model as input. To see a list of medaka models, use the command ```medaka tools list_models```.

As well, it is important to note, the models will *not* contain individual versions of Guppy. Thus, you should choose the version most closest to and less than the version you used (Ex: using Guppy 6.00, enter R***_min_hac_g507)

### Medaka Batch Size

A known issue with medaka is that it can use too much GPU memory and crash resulting in an incomplete assembly. If this happens, you can try entering the following command to allow for more GPU memory to be allotted:
```
export TF_FORCE_GPU_ALLOW_GROWTH=true
```
If this still does not solve the issue, you can reduce the batch size. The pipeline allows for the modification of this value through an optional argument ```--medakaBatchSize INT```. The value is defaultly 100.
## Installation

To install this pipeline, enter the following commands:
```
# Clone the repository
git clone https://github.com/rchapman2000/ont-de-novo-assembly.git

# Create a conda environment using the provided environment.yml file
conda env create -f environment.yml

# Activate the conda environment
conda activate ONT-DeNovoAssembly
```

## Updating the Pipeline
If you already have the pipeline installed, you can update it using the following commands:
```
# Navigate to your installation directory
cd ont-de-novo-assembly

# Use git to pull the latest update
git pull

# Activate the conda environment and use the environment.yml file to download updates
conda activate ONT-DeNovoAssembly
conda env update --file environment.yml --prune
```

## Usage
To run the pipeline, use the following command:
```
# You must either be in the same directory as the main.nf file or reference the file location.
nextflow run main.nf [OPTIONS] --input INPUT_DIR --output OUTPUT_DIR --model MEDAKA_MODEL
```

### Optional Arguments
The pipeline also supports the following optional arguments:

| Option | Type | Description |
|---|---|---|
| --trimONTAdapters | *None* | Enables ONT Adapter/Barcode trimming using Porechop [Default = off] |
| --minReadLen | *int* | If supplied, the pipeline will perform length filtering using Chopper excluding reads less than this size [Default = off] |
| --maxReadLen | *int* | If supplied, the pipeline will perform length filtering using Chopper excluding reads larger than this size [Default = off] |
| --medakaBatchSize | *int* | Medaka uses a lot of GPU memory, and if you're assembly is large enough it may cause Medaka to crash. Reducing the batch size will help solve this issue. [Default = 100] |
| --preGuppy5 | *None* | Flye handles data generated version of Guppy < 5.0 differently. Supply this parameter if your data was generated pre-Guppy 5.0 |
| --threads | *int* | The number of CPU threads that can be use to run pipeline tools in parallel |

To view the list of options from the command line, use the following command:
```
nextflow run main.nf --help
```