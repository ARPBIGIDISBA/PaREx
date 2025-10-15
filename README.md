
# *Pa*REx: an open‑source pipeline for the automated analysis of *Pseudomonas aeruginosa* resistomes from whole‑genome sequences

The *Pseudomonas aeruginosa* Resistome Explorer (*Pa*REx) is an open-source Python-based customizable pipeline that has been specifically designed for the automated analysis of *P. aeruginosa* resistomes from Illumina® paired-end reads. *Pa*REx uses different open-source bioinformatics tools, software and publicly available databases along with custom-built databases, scripts and tools and is composed by two main components the PaREx pipeline and the PaREx databases.

## Requirements

Parex requires:
- Git
- Python 3.8 or higher
- Other packages listed in `requirements.txt`

## Installation of *Pa*REx custom-built scripts and tools

Here are the steps to install *Pa*REx from the GitHub repository:

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/project.git  CAMBIAR POR: https://github.com/ARPBIG/PaREx.git
   cd parex
   ```

2. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

## Installation of *Pa*REx custom-built databases

There is a repository including all custom-built databases [parex-databases](https://github.com/matiasbonet/References4alignment  CAMBIAR POR: https://github.com/ARPBIG/parex-databases.git). 

 ```bash
   cd ..
   git clone https://github.com/yourusername/project.git  CAMBIAR POR: https://github.com/ARPBIG/parex-databases.git
   cd parex
   ```

## Installation of third-party bioinformatic tools, software and databases 
We have included an `install.sh` script to facilitate the installation of the required third-party tools and databases. Whith the versions tested for *Pa*REx. You can use the install.sh:  

   ```bash
   chmod 755 install.sh
   sh install.sh 
   ```


OR alternatively you can install them by your own or use your installed versions

1. **Trimmomatic**
   - Download and extract from the official source:
     ```bash
     wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
     unzip Trimmomatic-0.39.zip
     ```

2. **SPAdes**
   - Download and extract the Linux version from the official source:
     ```bash
     CAMBIAR URL! :  wget http://cab.spbu.ru/files/release3.15.3/SPAdes-3.15.3-Linux.tar.gz
     tar -xzf SPAdes-3.15.3-Linux.tar.gz
     ```

3. **Snippy**
   - Install Snippy by cloning its GitHub repository:
     ```bash
     git clone https://github.com/tseemann/snippy.git
     ```

4. **MLST**
   - MLST requires a specific download setup based on system. Instructions are available at its [official source](https://github.com/tseemann/mlst).
   - *Pa*REx validation was performed with the *P. aeruginosa* database available on November 2024.

5. **Resfinder**
   - Resfinder requires a specific download setup based on system. Instructions are available at its [official source](https://bitbucket.org/genomicepidemiology/resfinder.git)
   - *Pa*REx uses the version 2.1.0 of the resfinder database.
     

## Configuration 

We have to configure the paths for the databases and tools used in the pipeline.
Configure the `.json` files are located in the `programs_scripts/configs/` folder. 
You can see samples files in the folder `/programs_scripts/configs/samples/*.json.sample`.
You need to adjust the paths for the following configuration files: 

### Configuration Files

- **general.json**: Contains global settings, including the path to the projects and to PaREx databases.
- **trimmomatic.json**, **SPAdes.json**, **snippy.json**, **mlst.json**, **resfinder.json**: Specific configurations for each third-party tool used in the parex software.


## *Pa*REx general usage 

 --> INTRODUCIR create_project, create_sample_list, novaseq, trimmomatic, resistome, single file! 

### How to run the pipeline

The pipeline is designed to be run from the command line. You can execute individual steps or the entire pipeline as needed.

Once you have defined the configuration files, you can start using the pipeline.
The important folders to consider are:
- **PROJECTS_PATH**: Defined in `general.json`, this is where all project data
   will be stored.
- **PROJECT_NAME**: Each project will have its own folder within `PROJECTS_PATH`.
- **SAMPLES_PATH**: Each project will have a `samples/` folder where the input FASTQ files should be placed. See section input files for more details.
### Running the Complete Resistome Analysis Pipeline

The resistome pipeline can be executed from the command line with the following parameters:

python parex.py PROJECT_NAME "operation1,operation2,..."

where project_name is the name of the project (which will be created in the `PROJECTS_PATH` defined in `general.json`), and operation1, operation2, etc. are the specific operations to be performed. You can specify multiple operations separated by commas. 

First command will be to create the project structure:

```bash
python execute_pipeline.py PROJECT_NAME create_project
```

Then you add you sample FASTQ files to the `samples/` folder within the created project directory with the PROJECT_NAME as specified above.

Then, you can create the sample list automatically with the following command:

```bash
python execute_pipeline.py PROJECT_NAME create_sample_list
```

Finally, you can run the desired operations, in this case resistome analysis:

```bash
python execute_pipeline.py PROJECT_NAME resistome
```

### Optional Arguments
You can also specify optional arguments:
- `--log-level`: Set the logging level (e.g., `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). Default is `INFO`.
- `--force`: Force the execution of the program even if previous steps have not been completed.
- `--keep_output`: Keep the output of the program.
- `--clean_output`: Clean the output of the program.
- `--file`: Direct Path to the file not use sample list

### Special operations:

you can add **trimmomatic** to the process of resistome analysis:

```bash
python execute_pipeline.py PROJECT_NAME trimomatic,resistome
```

You can also unzip the files if you have them in `.gz` to speed up the process:

```bash
python execute_pipeline.py PROJECT_NAME unzip,resistome
```

You can also run the analysis from NovaSeq files:

```bash
python execute_pipeline.py PROJECT_NAME novaseq,resistome
```

Or you can you can run all together:

```bash
python execute_pipeline.py PROJECT_NAME novaseq,trimmomatic,resistome

```


## Logs

A folder with all the logs will be generated in the `logs/` folder. Each project will have its own log file named `Pparex_execution.log`.  
based on your defined `--log-level` argument.

## Contact

For questions or issues, open an *issue*

## Code Structure

- **programs_scripts/**: Contains specific scripts for each step of the pipeline.
- **programs_scripts/configs/**: Configuration files in `.json` format necessary for running the scripts. Also includes `.json.sample` example files.
- **execute_pipelines**: Run with `python execute_pipelines.py PROJECT_NAME "commands"`

## Generated Folders
When you run a new project the following folder structure will be generated within the specified `PROJECTS_PATH`:
- **PROJECT_NAME/**: Main project folder named after the project.
   - **FAST_Q_$PROJECT_NAME/**: Folder where input FASTQ files should be placed.
   - **ANALYSIS_$PROJECT_NAME/**: Folder where all results from the analysis will be stored.
      - **denovo_assemblies_SPAdes/**: 
      - **resfinder_results/**: Subfolder for ResFinder partial results.
      - **mlst_results/**: Subfolder for MLST partial results.
      - **oprd_results/**: Subfolder for OPRD partial results.
      - **snippy_results/**: Subfolder for Snippy partial results.
      - **trimmomatic_results/**: Subfolder for Trimmomatic partial results.
      - **PDF_results/**: Subfolder for PDF reports.
   - **logs/**: Subfolder for logs related to the analysis.
   - **$PROJECT_NAME_summary.xlsx**: Summary Excel file containing results from all analyses.

Each process generates a excel output individual and also a full summary excel file in the `ANALYSIS_$PROJECT_NAME/` folder.

The global summary file is named `$PROJECT_NAME_summary.xlsx`. In the PDF folder you will have a result for every inidiviaul sample.


## Input Files

The project is designed to handle a variety of input file formats, ensuring flexibility and compatibility with different sequencing workflows:

- **FASTQ Files**: Raw sequencing reads in FASTQ format are the primary input for processing workflows. These files typically contain paired-end reads (`R1` and `R2`) generated from Illumina or similar platforms.

- **NovaSeq Files**: High-throughput NovaSeq files are supported, allowing seamless integration with modern sequencing technologies. It converts in the output into FASTQ files R1/R2

- **De Novo Assembly Files**: For workflows that bypass raw read processing, preassembled de novo files can be directly utilized. This is particularly useful for tools and analyses that focus on assembled genomes, reducing the computational overhead of preprocessing.

This multi-format compatibility enables the project to adapt to different experimental setups, whether working with raw sequencing data or preassembled genomes, enhancing its usability and flexibility across genomic research pipelines.


Log files are stored in the `logs/` folder. They are automatically generated with each run and saved in the `Logs/` subdirectory within each project.

 in the repository or contact the development team at [matiasbonet@oceandrivers.com](mailto:matiasbonet@oceandrivers.com). [carla.lopez@ssib.es](mailto:carla.lopez@ssib.es)
