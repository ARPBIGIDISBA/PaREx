'''
    Este script se utiliza para aplicar el programa Trimmomatic a los ficheros fastq.gz
    Tiene como entrada los ficheros fastq.gz de las muestras
    Da como resultado dos ficheros R1_001.fastq.gz y R2_001.fastq.gz escritos en OUT
'''

import os
from modules.general_functions import read_args, execute_command

# Leer los argumentos de la línea de comandos y el fichero de configuración
PROJECT_NAME, samples, config, logging = read_args("trimmomatic_config.json")

logger = logging.getLogger(__name__)

# Parametros de configuración de este script
PROJECTS_PATH = config['PROJECTS_PATH']
TRIMMOMATIC_JAR_PATH = config['TRIMMOMATIC_JAR_PATH']

# Aqui puedes añadir opciones a trimomatic http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
TRIMMOMATIC_OPTIONS = config['TRIMMOMATIC_OPTIONS']

# Crear el directorio para el lineage
PROJECT_PATH = os.path.join(PROJECTS_PATH, PROJECT_NAME)
os.makedirs(PROJECT_PATH, exist_ok=True)

OUTPUT_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{PROJECT_NAME}", "FASTQ_Trimmomatic")
os.makedirs(OUTPUT_PATH, exist_ok=True)

for sample_name in samples:
    # Limpiar por si hay espacios en blanco
    sample_name = sample_name.strip()
    logger.info(f"Processing {sample_name}")
    # Crear el directorio para el sample
    SAMPLE_PATH = os.path.join(PROJECT_PATH, sample_name)
    os.makedirs(SAMPLE_PATH, exist_ok=True)

    # Crear los paths de entrada y salida
    input_r1_path = os.path.join(PROJECT_PATH, f"FASTQ_{PROJECT_NAME}", f"{sample_name}_R1_001.fastq.gz")
    input_r2_path = os.path.join(PROJECT_PATH, f"FASTQ_{PROJECT_NAME}", f"{sample_name}_R2_001.fastq.gz")

    # Example of output_files = /home/micro/Analysis/Trimmomatic/lineage/sample/{line}.trimmed.1P.fastq.gz
    output_files = [os.path.join(OUTPUT_PATH, f"{sample_name}.trimmed.{file}.fastq.gz") for file in ["1P", "1U", "2P", "2U"]]
    
    # Ejecutar Trimmomatic
    command = ["java", "-jar", TRIMMOMATIC_JAR_PATH, "PE", "-phred33",
               input_r1_path, input_r2_path] + output_files + TRIMMOMATIC_OPTIONS
    
    result = execute_command(command, logger)

    if result:
        # Renombrar los ficheros de salida de 1P y 2P a R1_001 y R2_001
        for suffix, new_suffix in [("1P", "R1"), ("2P", "R2")]:
            old_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}.trimmed.{suffix}.fastq.gz")
            new_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}_trim_{new_suffix}.fastq.gz")
            os.rename(old_file_path, new_file_path)
            logger.info(f"Renaming file {new_file_path}")
            os.system(f"gunzip {new_file_path}")
            logger.info(f"Unzip file {new_file_path}")

        # Mover los unpairs para futura calidad un directorio
        UNPAIRED_PATH = os.path.join(OUTPUT_PATH, "UNPAIRED")
        os.makedirs(UNPAIRED_PATH, exist_ok=True)
        for suffix in ["1U", "2U"]:
            old_file_path = os.path.join(OUTPUT_PATH, f"{sample_name}.trimmed.{suffix}.fastq.gz")
            new_file_path = os.path.join(UNPAIRED_PATH, f"{sample_name}.trimmed.{suffix}.fastq.gz")
            os.rename(old_file_path, new_file_path)
            logger.info(f"Storing unpaired file {new_file_path}")

    else:
        logger.error("There is an error in Trimmomatic execution, check the log files")