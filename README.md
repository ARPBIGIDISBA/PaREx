 
# *Pa*REx: an open‑source pipeline for the automated analysis of *Pseudomonas aeruginosa* resistomes from whole‑genome sequences

The *Pseudomonas aeruginosa* Resistome Explorer (*Pa*REx) is an open-source, Python-based pipeline designed for the automated analysis of *P. aeruginosa* resistomes from Illumina® paired-end whole-genome sequencing reads. 

*Pa*REx integrates multiple open-source bioinformatics tools, software and publicly available databases along with custom-built databases and scripts. It consists of two main components: the *Pa*REx pipeline and the *Pa*REx databases.

## Requirements

Parex requires:
- Git
- Python 3.8 or higher
- Other packages listed in `requirements.txt`

## Installation guide of *Pa*REx 

1. Clone the repository:

   ```bash
   git clone https://github.com/ARPBIGIDISBA/parex.git
   cd parex
   ```

2. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

## Installation guide of *Pa*REx databases

The *Pa*REx databases are hosted in a separate repository:[parex-databases](https://github.com/ARPBIGIDISBA/parex-databases.git). 

 ```bash
   cd ..
   git clone https://github.com/ARPBIGIDISBA/parex-databases.git
   cd parex
   ```

## Installation guide of third-party bioinformatic tools

*Pa*REx includes an `install.sh` script to facilitate installation of tested third-party tools and databases.

   ```bash
   chmod 755 install.sh
   sh install.sh 
   ```

OR alternatively you can install them manually as described below.

1. **Trimmomatic**
   - Download and extract:
     ```bash
     wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
     unzip Trimmomatic-0.39.zip
     ```

2. **SPAdes (v3.15.3 – tested version)**
   - Download and extract:
     ```bash
     wget http://cab.spbu.ru/files/release3.15.3/SPAdes-3.15.3-Linux.tar.gz
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

Before running PaREx, you must configure the paths to all required databases and third-party tools.

Configuration `.json` files are located in the `programs_scripts/configs/` folder. 
Sample configuration templates are available at: `/programs_scripts/configs/samples/*.json.sample`.
Copy the sample files, remove the .sample extension, and adjust the paths according to your system.

### Configuration Files

You must update the following files:

- **general.json**: Contains global settings, including the PROJECTS_PATH (Directory where all projects will be created) and the path to PaREx custom databases
- **trimmomatic.json**, **SPAdes.json**, **snippy.json**, **mlst.json**, **resfinder.json**: Each of these files contains the executable path and relevant parameters for the corresponding tool.

## *Pa*REx general usage 

The pipeline is fully command-line based and modular. You can execute individual steps or chain multiple operations together.

Important Project Variables:

- **PROJECTS_PATH**: Defined in general.json. This is the root directory where all project data will be stored.
- **PROJECT_NAME**: Each new project will have its own folder inside `PROJECTS_PATH`.
- **SAMPLES_PATH**: Each project will have a `samples/` folder where the input FASTQ files should be placed. 

### How to run the *Pa*REx Resistome Explorer Pipeline: standard workflow

General syntax: python parex.py PROJECT_NAME "operation1,operation2,..."

where project_name is the name of your project (which will be created in the `PROJECTS_PATH` defined in `general.json`), and operation1, operation2, etc. are the list of  operations. 

First operation: Create the Project Structure

Before running any analysis, you must create a new project. This step creates the complete folder structure required by *Pa*REx inside the directory defined as PROJECTS_PATH in general.json.

```bash
python parex.py PROJECT_NAME create_project
```

Second operation: Create the Sample List

Once the project structure has been created, place your input FASTQ files inside the corresponding FAST_Q folder and create the sample list automatically with the following command:

```bash
python parex.py PROJECT_NAME create_sample_list
```

Third operation: Run the Resistome Analysis 

```bash
python parex.py PROJECT_NAME resistome
```

### Optional Arguments
You may append optional flags:
- `--log-level`: Set the logging level (e.g., `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`). Default is `INFO`.
- `--force`: Force the execution of the program even if previous steps have not been completed.
- `--keep_output`: Keep the output of the program.
- `--clean_output`: Clean the output of the program.
- `--file`: Direct Path to the file not use sample list

### Special operations:

Merge L1 and L2 NovaSeq files:

```bash
python parex.py PROJECT_NAME novaseq,resistome
```

Unzip `.gz` before analysis to speed up the process:

```bash
python parex.py PROJECT_NAME unzip,resistome
```

Include **trimmomatic** preprocessing:

```bash
python parex.py PROJECT_NAME trimomatic,resistome
```

Or you can you can run all together:

```bash
python parex.py PROJECT_NAME novaseq,unzip,trimmomatic,resistome

```

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

## How to cite

If you use this software in your work, please cite the following:

López-Causapé C, et al. *Pa*REx: an open-source pipeline for the automated analysis of Pseudomonas aeruginosa resistomes from whole-genome sequences. Antimicrob Agents Chemother 0:e01326-25.https://doi.org/10.1128/aac.01326-25

## License

This pipeline is distributed under the [CC BY-NC-ND 4.0 License](http://creativecommons.org/licenses/by-nc-nd/4.0/).

## Contact

For questions or issues, open an *issue* in the repository or contact the development team at 
[carla.lopez@ssib.es](mailto:carla.lopez@ssib.es) [matiasbonet@oceandrivers.com](mailto:matiasbonet@oceandrivers.com). 
