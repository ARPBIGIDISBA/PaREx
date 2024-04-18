'''
    Este script se utiliza para aplicar el programa Trimmomatic a los ficheros fastq.gz
    Tiene como entrada los ficheros fastq.gz de las muestras
    Da como resultado dos ficheros R1_001.fastq.gz y R2_001.fastq.gz escritos en OUT
'''

import os
import argparse
import logging
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs
import pysam
import csv
import re

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)

config = init_configs(script_directory, "snippy.json")

amino_acids = {
    'Ala': 'A', 'Gly': 'G', 'Met': 'M', 'Ser': 'S',
    'Cys': 'C', 'His': 'H', 'Asn': 'N', 'Thr': 'T',
    'Asp': 'D', 'Ile': 'I', 'Pro': 'P', 'Val': 'V',
    'Glu': 'E', 'Lys': 'K', 'Gln': 'Q', 'Trp': 'W',
    'Phe': 'F', 'Leu': 'L', 'Arg': 'R', 'Tyr': 'Y'
}

def translate_amino_acid(value):
    # Remove p. select three letters before number and after number translate "p.Asp104Glu"
    match = re.findall(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})", value)
    if match:
        before_number, number, after_number = match[0]
        before_number = amino_acids.get(before_number.capitalize(), None)
        number = number
        after_number = amino_acids.get(after_number.capitalize(), None)
        return f"p.{before_number}{number}{after_number}"
    else:
        print("No match found", value)
        return value
    
    
def process_output(vcf_path, sample_name, output_path):
    output_dir = os.path.join(output_path, "processed")
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, f"{sample_name}_snippy.csv")
    with open(csv_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        # Abrir el archivo VCF con pysam
        vcf_file = pysam.VariantFile(vcf_path)
        csv_writer.writerow(['Titulo 1', 'Titulo 2', 'Titulo 3'])

        # Iterar sobre cada registro en el archivo VCF
        for record in vcf_file:
            if 'ANN' in record.info:
                annotations = record.info['ANN']
                # Each annotation can be split by comma, and further by pipe '|'
                for annotation in annotations:
                    fields = annotation.split('|')
                    # Check if the impact field matches 'MODERATE'
                    
                    if len(fields) > 1 and fields[1] == 'missense_variant' and fields[2]!='LOW' and fields[2]!='MODIFIER':
                        if fields[3] in config["CEPAS"]:
                            
                            field = translate_amino_acid(fields[10])
                            row = [fields[3], fields[9], field]
                            print(row)
                            csv_writer.writerow(row)
                        
def snippy_run(project_name, only_output=False,  config=config):
    ''' 
        this function is used to apply the Trimmomatic program to the fastq.gz files

        parameters:
            project_name (str): Name of the project<<3
            config_json (str): Path to the config file by default is trimmomatic_config.json

        results:
            Two files R1_001.fastq.gz and R2_001.fastq.gz written in 
            PROJECTS_PATH/project_name/ANALYSIS_{project_name}/FASTQ_Trimmomatic
            filename format will be  sample_name_trim_R1.fastq.gz and sample_name_trim_R2.fastq.gz
    '''

    # Leer las muestras y ficheros de configuracion
    samples = read_args(project_name, config)

    # Parametros de configuración de este script
    PROJECTS_PATH = config['PROJECTS_PATH']
    SNIPPY_PATH = config['SNIPPY_PATH']

    # Aqui puedes añadir opciones a trimomatic http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
    SNIPPY_OPTIONS = config['SNIPPY_OPTIONS']
    SNIPPY_REFERENCE = config['REFERENCE']

    # Create project directory in case it is not created
    PROJECT_PATH = os.path.join(PROJECTS_PATH, project_name)
    os.makedirs(PROJECT_PATH, exist_ok=True)

    OUTPUT_PATH = os.path.join(PROJECT_PATH, f"ANALYSIS_{project_name}", "snippy_results")
    os.makedirs(OUTPUT_PATH, exist_ok=True)

    for sample_name in samples:
        # Limpiar por si hay espacios en blanco
        sample_name = sample_name.strip()
        logger.info("Processing %s", sample_name)
        
        # Crear los paths de entrada y salida
        input_r1_path = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}", f"{sample_name}_R1_001.fastq.gz")
        input_r2_path = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}", f"{sample_name}_R2_001.fastq.gz")
        if not os.path.exists(input_r1_path) or not os.path.exists(input_r2_path):
            logger.error(f"The fastq.gz file for {sample_name} don't exist")
            logger.error(f"One of this files does not exist:\n {input_r1_path}\n {input_r2_path}")
        
        OUTPUT_RESULTS = os.path.join(OUTPUT_PATH, "output", sample_name)
        os.makedirs(OUTPUT_RESULTS, exist_ok=True)
        VCF_FILE = os.path.join(OUTPUT_RESULTS, f"snps.vcf")
        if only_output:
            
            if os.path.exists(VCF_FILE):
                result = True
            else:
                logger.error("You have to run first the resfinder process")
                logger.error("File not found: %s", VCF_FILE)
                return False
        else:
            # Ejecutar Trimmomatic
            command = [SNIPPY_PATH, "--outdir", OUTPUT_RESULTS, "--ref", SNIPPY_REFERENCE, "--R1", input_r1_path, "--R2", input_r2_path] + SNIPPY_OPTIONS
            result = execute_command(command)
        
        if result:
            logger.info("Snippy process for %s finished", sample_name)
            process_output(VCF_FILE, sample_name, OUTPUT_PATH)


if __name__ == "__main__":
    
    # Define the arguments that the program expects
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--parse-output', action='store_true', help='Set the flag to not execute but only process output files')
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    args = parser.parse_args()
    PROJECT_NAME = args.PROJECT_NAME
    if args.json_config:
        config = init_configs(script_directory, args.json_config)   

    configure_logs(PROJECT_NAME, "snippy", config)
    logger = logging.getLogger(__name__)
    snippy_run(PROJECT_NAME, args.parse_output, config)
