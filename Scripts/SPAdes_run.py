'''
    Este script aplica el ensamblaje de novo con SPAdes a los ficheros fastq.gz
    Ejecuta el programa en python spades sobre los ficheros fastq
    Tiene como entrada los ficheros fastq.gz de las muestras
    Da como resultado un fichero fasta con los SPAdes.denovoassembly.fasta
    Here are the command options for spades https://github.com/ablab/spades#sec3.2

'''
import os
from modules.general_functions import read_args, execute_command


# Read command line arguments, sample list and config file
PROJECT_NAME, samples, config, logging = read_args("SPAdes_config.json")

logger = logging.getLogger(__name__)

# Parametros de configuración de este script
PROJECTS_PATH = config['PROJECTS_PATH']

# list of coma separated options https://github.com/ablab/spades#sec3.2
SPADES_PROGRAM_PATH = config['SPADES_PATH']
SPADES_OPTIONS = config['SPADES_OPTIONS']


# Crear el directorio para el lineage
PROJECT_PATH = os.path.join(PROJECTS_PATH, PROJECT_NAME)
os.makedirs(PROJECT_PATH, exist_ok=True)

TRIMMOMATIC_FILES_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{PROJECT_NAME}", "FASTQ_Trimmomatic")

OUTPUT_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{PROJECT_NAME}", "denovo_assemblies_SPAdes")
os.makedirs(OUTPUT_PATH, exist_ok=True)

for sample_name in samples:
    # Limpiar por si hay espacios en blanco
    sample_name = sample_name.strip()
    logger.info(f"Processing {sample_name}")

    # Definir los ficheros de entrada 1 y 2 

    input_r1_path = os.path.join(TRIMMOMATIC_FILES_PATH, f"{sample_name}_trim_R1_001.fastq")
    input_r2_path = os.path.join(TRIMMOMATIC_FILES_PATH, f"{sample_name}_trim_R2_001.fastq")

    command = ["python", SPADES_PROGRAM_PATH, "-o", OUTPUT_PATH, "-1", input_r1_path, "-2", input_r2_path] + SPADES_OPTIONS
    
    result = execute_command(command, logger)

    if result:
        logger.info("SPAdes assembly finished")
        # Renombrar los archivos de salida
        old_file_path = os.path.join(OUTPUT_PATH, "contigs.fasta")
        new_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}.SPAdes.denovoassembly.fasta")

        os.rename(old_file_path, new_file_path)
        logger.info(f"Rename files for other analysis {new_file_path}")
    else:
        logger.error("SPAdes assembly failed")