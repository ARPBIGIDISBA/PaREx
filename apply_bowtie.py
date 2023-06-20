'''
    Este script se usa para hacer una nalisis bowtie
    Tiene como entrada los ficheros fastq.gz de las muestras
    
'''
import os
from modules.general_functions import read_args_bowtie, execute_command

# Parametros de configuración de este script
BASE_PATH = '/home/micro/Analysis/Trimmomatic'
BOWTIE_PATH = "/home/biel/Programs/bowtie2-2.2.6/bowtie2"
INPUT_PATH = '/home/micro/Analysis/assemblies/fastq-MiSeq'
# Aqui puedes añadir opciones a trimomatic http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
BOWTIE_OPTIONS = ["-X", "1000", "-S"]
LIBS_PATH = '/home/biel/indexed_libs'

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