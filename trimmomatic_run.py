'''
    Autor: Matias Bonet Fullana 
    Este script se utiliza para aplicar el programa Trimmomatic a los ficheros fastq.gz
    Tiene como entrada los ficheros fastq.gz de las muestras
    Da como resultado dos ficheros R1_001.fastq.gz y R2_001.fastq.gz escritos en OUT
'''
import os
from modules.general_functions import read_args, execute_command
import json


# Leer el archivo de configuración
with open('trimmomatic_config.json', 'r') as file:
    config = json.load(file)

# Parametros de configuración de este script
BASE_PATH = config['BASE_PATH']
TRIMMOMATIC_JAR_PATH = config['TRIMMOMATIC_JAR_PATH']
INPUT_PATH = config['INPUT_PATH']

# Aqui puedes añadir opciones a trimomatic http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
TRIMMOMATIC_OPTIONS = config['TRIMMOMATIC_OPTIONS']

# Leer los argumentos de la línea de comandos y el fichero de configuración
SAMPLES, LINEAGE, lines, logging = read_args()

# Constrantes para el nombre de los ficheros de salida
files = ["1P", "1U", "2P", "2U"]

# Crear el directorio para el lineage
LINEAGE_PATH = os.path.join(BASE_PATH, LINEAGE)
os.makedirs(LINEAGE_PATH, exist_ok=True)

for line in lines:
    # Limpiar por si hay espacios en blanco
    line = line.strip()

    # Crear el directorio para el sample
    SAMPLE_PATH = os.path.join(LINEAGE_PATH, line)
    os.makedirs(SAMPLE_PATH, exist_ok=True)

    # Crear los paths de entrada y salida
    input_r1_path = os.path.join(INPUT_PATH, LINEAGE, f"{line}_R1_001.fastq.gz")
    input_r2_path = os.path.join(INPUT_PATH, LINEAGE, f"{line}_R2_001.fastq.gz")

    # Example of output_files = /home/micro/Analysis/Trimmomatic/lineage/sample/{line}.trimmed.1P.fastq.gz
    output_files = [os.path.join(SAMPLE_PATH, f"{line}.trimmed.{file}.fastq.gz") for file in files]
    
    # Ejecutar Trimmomatic
    command = ["java", "-jar", TRIMMOMATIC_JAR_PATH, "PE", "-phred33",
               input_r1_path, input_r2_path] + output_files + TRIMMOMATIC_OPTIONS
    
    execute_command(command, logging)

    # Renombrar los ficheros de salida de 1P y 2P a R1_001 y R2_001
    for suffix, new_suffix in [("1P", "R1_001"), ("2P", "R2_001")]:
        old_file_path = os.path.join(SAMPLE_PATH, f"{line}.trimmed.{suffix}.fastq.gz")
        new_file_path = os.path.join(LINEAGE_PATH, f"{line}_trim_{new_suffix}.fastq.gz")
        os.rename(old_file_path, new_file_path)