'''
    Este script se usa para hacer una nalisis bowtie
    Tiene como entrada los ficheros fastq.gz de las muestras
    
'''
import os
from modules.general_functions import read_args_bowtie, execute_command
import json

# Leer el archivo de configuración
with open('bowtie_config.json', 'r') as file:
    config = json.load(file)

# Parametros de configuración de este script
BASE_PATH = config['BASE_PATH']
BOWTIE_PATH = config['BOWTIE_PATH']
LIBS_PATH = config['LIBS_PATH']
INPUT_PATH = config['INPUT_PATH']
BOWTIE_OPTIONS = config['BOWTIE_OPTIONS']

# Leer los argumentos de la línea de comandos y el fichero de configuración
SAMPLES, LINEAGE, GENE, lines, logging = read_args_bowtie()

# Constrantes para el nombre de los ficheros de salida
files = ["1P", "1U", "2P", "2U"]

# Crear el directorio para el lineage
LINEAGE_PATH = os.path.join(BASE_PATH, LINEAGE)
os.makedirs(LINEAGE_PATH, exist_ok=True)

for line in lines:
    # Limpiar por si hay espacios en blanco
    line = line.strip()
    logging.info(f"Processing sample {line}")

    gene_index_path = os.path.join(LIBS_PATH, GENE, GENE)
    gene_index_path = os.path.join(LIBS_PATH, GENE, GENE)
    input_r1_path = os.path.join(INPUT_PATH, LINEAGE, 'new_fastq', f"new_{line}_L001_R1_001.fastq")
    input_r2_path = os.path.join(INPUT_PATH, LINEAGE, 'new_fastq', f"new_{line}_L001_R2_001.fastq")
    output_path = os.path.join(INPUT_PATH, LINEAGE, 'sam_files', f"{line}_map{GENE}.sam")

    command = [BOWTIE_PATH, "--phred33", "-x", gene_index_path, "-q", "-1", input_r1_path, "-2", input_r2_path] + BOWTIE_OPTIONS+ [output_path]


    execute_command(command, logging)

    logging.info(f"Sample {line} processed")º