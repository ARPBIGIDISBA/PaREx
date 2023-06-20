'''
    Este script aplica el ensamblaje de novo con SPAdes a los ficheros fastq.gz
    Ejecuta el programa en python spades sobre los ficheros fastq
    Tiene como entrada los ficheros fastq.gz de las muestras
    Da como resultado un fichero fasta con los SPAdes.denovoassembly.fasta
'''
import os
from modules.general_functions import get_args, execute_command

# Parametros de configuración de este script el resto no tocar
BASE_PATH = ' /home/micro/Analysis/assemblies/denovo_assemblies_SPAdes/'
SPADES_PATH = "/home/micro/Programs/SPAdes-3.15.0/bin/spades.py"
ASSEMBLIES_PATH = '/home/micro/Analysis/assemblies/fastq-MiSeq'
SPADES_OPTIONS = ["--careful"] # list of coma separated options https://github.com/ablab/spades#sec3.2
OUTPUT_PATH = '/home/micro/Analysis/assemblies/denovo_assemblies_SPAdes'


# Leer los argumentos de la línea de comandos y el fichero de configuración
SAMPLES, LINEAGE, lines, logging = get_args()


# Crear el directorio para el lineage
dir_path = os.path.join(BASE_PATH, LINEAGE)
os.makedirs(dir_path, exist_ok=True)

for line in lines:
    # Limpiar por si hay espacios en blanco
    line = line.strip()

    # Crear el directorio para el sample
    subdir_path = os.path.join(dir_path, line)
    os.makedirs(subdir_path, exist_ok=True)

    # Definir los ficheros de entrada 1 y 2 
    input_r1_path = os.path.join(ASSEMBLIES_PATH, LINEAGE, f"{line}_R1_001.fastq")
    input_r2_path = os.path.join(ASSEMBLIES_PATH, LINEAGE, f"{line}_R2_001.fastq")

    command = ["python", SPADES_PATH, "-o", subdir_path, "-1", input_r1_path, "-2", input_r2_path] + SPADES_OPTIONS 
    
    execute_command(command, logging)

    # Renombrar los archivos de salida
    old_file_path = os.path.join(subdir_path, "contigs.fasta")
    new_file_path = os.path.join(OUTPUT_PATH, LINEAGE, f"{line}.SPAdes.denovoassembly.fasta")

    os.rename(old_file_path, new_file_path)
