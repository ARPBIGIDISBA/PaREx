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
import pandas as pd

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)

config = init_configs(script_directory, "snippy.json")

amino_acids = {
    'Ala': 'A', 'Gly': 'G', 'Met': 'M', 'Ser': 'S',
    'Cys': 'C', 'His': 'H', 'Asn': 'N', 'Thr': 'T',
    'Asp': 'D', 'Ile': 'I', 'Pro': 'P', 'Val': 'V',
    'Glu': 'E', 'Lys': 'K', 'Gln': 'Q', 'Trp': 'W',
    'Phe': 'F', 'Leu': 'L', 'Arg': 'R', 'Tyr': 'Y',
    'fs':'X', 'del':'del', 'ins':'ins', 'dup':'dup'
}

def translate_amino_acid(value, value_c=""):
    # Remove p. select three letters before number and after number translate "p.Asp104Glu"
    
    try:
        if value.find("del") > 0 or value.find("ins") > 0:
            value = value.replace("p.", "")
            # Check if it is a deletion
            if value.find("del") > 0:
                replace = "del"
            if value.find("ins") > 0:
                replace = value[value.find("ins"):]
            
            value = value.replace(replace, "")
            values = value.split("_")
            result = []
            for value in values:
                deletion = re.findall(r'([A-Za-z]{3})(\d+)', value)
                deletion = deletion[0]
                value = f"{amino_acids.get(deletion[0].capitalize(), None)}{deletion[1]}"
                result.append(value)
            return ["_".join(result)+replace]

        elif value.find("fs") > 0:
            value = value_c.replace("c.","")
            if value.find("del") > 0:
                # 240_247delGCCGGCCA add nt at the begining and remove letters after del nt240_247del
                value = [f"nt{value[:value.find('del')]}del"]
            elif value.find("ins") > 0:
                value = [f"nt{value}"]
            return value
        elif value.find("*") > 0:
            value = value.replace("p.","").replace("*","Stop")
            letter = amino_acids.get(value[0:3].capitalize(), None)
            value = f"{letter}{value[3:]}"
            return [value]
        else:
            parts = re.findall(r'(\d+)|([A-Za-z]{3})', value)
            if parts:
                previous = []
                after = []
                number = None
                for index, part in enumerate(parts):
                    if part[0]:
                        number = int(part[0])
                    if part[1]:
                        if number is None:
                            previous.append(amino_acids.get(part[1].capitalize(), None))
                        else:
                            after.append(amino_acids.get(part[1].capitalize(), None))

                result = []
                for index,part in enumerate(previous):
                    result.append(f"{previous[index]}{number}{after[index]}")
                    number+=1
                if len(result)>2:
                    result = [result[0], result[-1]]
                    
                return result
            else:
                print("No match found", value)
                return value
    except Exception as e:
        print(e)
        print(value)
        return value
    
def read_data_from_file(filename):
    hoja1_df = pd.read_excel(filename, sheet_name='All').fillna('')
    data = {}
    for index, row in hoja1_df.iterrows():
        locus_gene, polymorphisms = row
        if locus_gene.find("_")>0:
            locus, gene = locus_gene.split('_')
        else:
            locus = locus_gene
            gene = ""
            
        if gene == "-":
            gene = ""

        polymorphisms = polymorphisms.rstrip('.').split(',')
        cleaned_polymorphisms = [p.split('(')[0].strip() for p in polymorphisms]
        data[locus] = {'gene': gene, 'polymorphisms': cleaned_polymorphisms}

    return data
    
def process_output(vcf_path, sample_name, output_path):
    output_dir = os.path.join(output_path, "processed")
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, f"{sample_name}_snippy.csv")
    if os.path.exists(csv_path):
        os.remove(csv_path)
    path = config.get("MUTATION_REFERENCE_FILE")
    filter = read_data_from_file(path)

    with open(csv_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=";")
        # Abrir el archivo VCF con pysam
        vcf_file = pysam.VariantFile(vcf_path)
        csv_writer.writerow(['sample_name','gene','snippy_mutations', 'mutations', 'mutations filtered', 'C.'])
    
        # Iterar sobre cada registro en el archivo VCF
        for record in vcf_file:
            if 'ANN' in record.info:
                annotations = record.info['ANN']
                # Each annotation can be split by comma, and further by pipe '|'
                for annotation in annotations:
                    fields = annotation.split('|')
                    # Check if the impact field matches 'MODERATE'
                    if len(fields) > 1  and fields[2]!='LOW' and fields[2]!='MODIFIER':
                        locus = fields[3]
                        if locus in filter.keys(): #config["LOCUS"]:
                            mutations = translate_amino_acid(fields[10], fields[9])
                            
                            if locus in filter.keys():  
                                gene = filter[locus]['gene']
                                result_mutation = [m for m in mutations if m not in filter[locus]['polymorphisms']]
                            else:
                                result_mutation = mutations
                                gene = ""
                            row = [fields[3], gene,fields[10],",".join(mutations), ",".join(result_mutation),fields[9].replace("c.","")]
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
        if only_output or os.path.exists(VCF_FILE):
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
