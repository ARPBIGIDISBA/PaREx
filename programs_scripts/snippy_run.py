"""
This software is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).
More details: https://creativecommons.org/licenses/by-nc/4.0/"

This script is used to run the snippy program to get the polymorphisms of the samples
"""

import os
import argparse
import logging
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs
import pysam
import csv
import re
import pandas as pd
import traceback
import openpyxl

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)

config = init_configs(script_directory, "snippy.json", required_keys=["SNIPPY_PATH", "SNIPPY_OPTIONS", "REFERENCE", "POLYMORPHISMS"])

amino_acids = {
    'Ala': 'A', 'Gly': 'G', 'Met': 'M', 'Ser': 'S',
    'Cys': 'C', 'His': 'H', 'Asn': 'N', 'Thr': 'T',
    'Asp': 'D', 'Ile': 'I', 'Pro': 'P', 'Val': 'V',
    'Glu': 'E', 'Lys': 'K', 'Gln': 'Q', 'Trp': 'W',
    'Phe': 'F', 'Leu': 'L', 'Arg': 'R', 'Tyr': 'Y',
    'fs':'X', 'del':'del', 'ins':'ins', 'dup':'dup',
    "?": "?"
}

def update_dataframe(df, sample_name, name, value):
    # Ensure sample_name exists in the DataFrame index
    if sample_name not in df.index:
        df.loc[sample_name] = pd.Series(dtype=object)  # Create a new row for sample_name
    
    # Ensure name exists in the DataFrame columns
    if name not in df.columns:
        df[name] = pd.Series(dtype=object)  # Create a new column for name

    # Append or assign value
    existing_value = df.loc[sample_name, name]
    if pd.notna(existing_value):  # Append to existing value if it's not NaN
        df.loc[sample_name, name] = f"{existing_value}, {value}"
    else:  # Assign new value if cell is empty
        df.loc[sample_name, name] = value

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
                value = f"nt{value[:value.find('del')]}del"
            elif value.find("ins") > 0:
                value = f"nt{value}"
            elif value.find("dup") > 0:
                value = f"nt{value}"
            return [value]
        elif value.find("*") > 0:
            value = value.replace("p.","").replace("*","Stop")
            letter = amino_acids.get(value[0:3].capitalize(), None)
            value = f"{letter}{value[3:]}"
            return [value]
        elif value.find("?")>0:
            # Review with carla
            value = value.replace("p.","")
            return [value]
        else:
            parts = re.findall(r'(\d+)|([A-Za-z]{3}|\?)', value)
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
        traceback.print_exc()
        import sys
        sys.exit(1)
        return value


def read_data_from_file(filename):
    # Read All tab no matter if upper or lower case
    all_df = pd.read_excel(filename, sheet_name='All').fillna('')
    basic_df = pd.read_excel(filename, sheet_name='Basic').fillna('')

    files = [all_df, basic_df]
    keys = ['All', 'Basic']
    output = {}
    

    for key, df in enumerate(files):
        data = {}
        for index, row in df.iterrows():
            locus_gene, polymorphisms = row
            if locus_gene.find("_")>0:
                locus, gene = locus_gene.split('_')
            else:
                locus = locus_gene
                gene = ""

            if gene == "-":
                gene = ""

            polymorphisms = polymorphisms.rstrip('.').split(', ')
            cleaned_polymorphisms = [p.split(' (')[0].strip() for p in polymorphisms]
            data[locus] = {'gene': gene, 'polymorphisms': cleaned_polymorphisms}

        output[keys[key]] = data
    return output


def process_output(vcf_path, sample_name, output_path):
    output_dir = os.path.join(output_path, "processed")
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, f"{sample_name}_snippy.csv")
    if os.path.exists(csv_path):
        os.remove(csv_path)

    with open(csv_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=";")
        # Abrir el archivo VCF con pysam
        vcf_file = pysam.VariantFile(vcf_path)
        csv_writer.writerow(['locus', 'genes', 'P.', 'changes', 'filtered_mutations', 'C.'])

        # Iterar sobre cada registro en el archivo VCF
        path = config.get("POLYMORPHISMS")
        filter = read_data_from_file(path)
        filter_all = filter['All']


        results = []
        for record in vcf_file:
            if 'ANN' in record.info:
                annotations = record.info['ANN']
                # Each annotation can be split by comma, and further by pipe '|'
                for annotation in annotations:
                    fields = annotation.split('|')
                    # Check if the impact field matches 'MODERATE'
                    if len(fields) > 1  and fields[2]!='LOW' and fields[2]!='MODIFIER':
                        locus = fields[3]
                        if locus in filter_all.keys() :
                            p = fields[10]
                            c = fields[9]
                            changes = translate_amino_acid(p, c)
                            filter_mutations = [m for m in changes if m not in filter_all[locus]['polymorphisms']]
                            row = [locus, filter_all[locus]['gene'],p.replace("p.",""),  ",".join(changes), ",".join(filter_mutations), c.replace("c.","")]
                            results.append(row)
                            csv_writer.writerow(row)

def combined_excel_files(samples, output_path):
    output_dir = os.path.join(output_path, "processed")

    path = config.get("POLYMORPHISMS")
    filter = read_data_from_file(path)
    filter_all = filter['All']
    filter_basic = filter['Basic']
    csv_output = os.path.join(output_path, "combined_snippy.xlsx")

    columns_all = ["sample_name"]
    for locus in filter_all:
        name = f"{locus}_{filter_all[locus]['gene']}"
        columns_all.append(name)
    columns_basic = ["sample_name"]
    for locus in filter_basic:
        name = f"{locus}_{filter_basic[locus]['gene']}"
        columns_basic.append(name)
    df_all = pd.DataFrame(columns=columns_all)
    df_all_clean = pd.DataFrame(columns=columns_all)
    df_basic = pd.DataFrame(columns=columns_basic)
    df_basic_clean = pd.DataFrame(columns=columns_basic)

    # sample_name column is an index
    df_all.set_index('sample_name', inplace=True)
    df_all_clean.set_index('sample_name', inplace=True)
    df_basic.set_index('sample_name', inplace=True)
    df_basic_clean.set_index('sample_name', inplace=True)

    for sample_name in samples:
        sample_name = sample_name.strip()
        csv_path = os.path.join(output_dir, f"{sample_name.strip()}_snippy.csv")
        if not os.path.exists(csv_path):
            logger.warning("The file %s does not exist", csv_path)
            continue
        #read the csv file
        df = pd.read_csv(csv_path, delimiter=";")
        # select the columns to be added to the final file
        df = df[['locus', 'genes', 'changes', 'filtered_mutations']]
        # create a pandas with rows with different locus as columns
        # find column to insert value in df_all
        for index, row in df.iterrows():
            locus = row['locus']
            gene = row['genes']
            changes = row['changes']
            filtered_mutations = row['filtered_mutations']
            if locus in filter_all.keys():
                if name in df_all.columns:
                    if sample_name in df_all.index and pd.notna(df_all.loc[sample_name, name]):
                        changes = f"{changes},{df_all.loc[sample_name, name]}"
                    df_all.loc[sample_name, name] = changes
                    if pd.notna(filtered_mutations):
                        if sample_name in df_all_clean.index and pd.notna(df_all_clean.loc[sample_name, name]):
                            filtered_mutations = f"{filtered_mutations},{df_all_clean.loc[sample_name, name]}"
                        df_all_clean.loc[sample_name, name] = filtered_mutations
            if locus in filter_basic.keys():
                if name in df_basic.columns:
                    if sample_name in df_basic.index and pd.notna(df_basic.loc[sample_name, name]):
                        changes = f"{changes},{df_basic.loc[sample_name, name]}"
                    df_basic.loc[sample_name, name] = changes
                    if pd.notna(filtered_mutations):
                        if sample_name in df_basic_clean.index and pd.notna(df_basic_clean.loc[sample_name, name]):
                            filtered_mutations = f"{filtered_mutations},{df_basic_clean.loc[sample_name, name]}"
                        df_basic_clean.loc[sample_name, name] = filtered_mutations

    if os.path.exists(csv_output):
        os.remove(csv_output)

    # workbook = openpyxl.Workbook()
    # workbook.save(csv_output)
    workbook = openpyxl.Workbook()
    workbook.save(csv_output)
            
    # Cargar el archivo existente sin sobrescribir
    with pd.ExcelWriter(csv_output, engine='openpyxl') as writer:
        df_all.to_excel(writer, sheet_name='All', index=True)
        df_all_clean.to_excel(writer, sheet_name='All_clean', index=True)
        df_basic.to_excel(writer, sheet_name='Basic', index=True)
        df_basic_clean.to_excel(writer, sheet_name='Basic_clean', index=True)
        

def snippy_run(project_name, only_output=False,  config=config, extra_config={"force": False, "keep_output": False}):
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

    direct_file = None
    if extra_config["file"] is not None:
        direct_file = extra_config["file"]
    
    # Read command line arguments, sample list and config file  or direct file
    if not direct_file:
        samples = read_args(project_name, config)
    else:
        samples = [direct_file]

    if direct_file is None and extra_config["file"] is not None:
        direct_file = extra_config["file"]
    # Leer las muestras y ficheros de configuracion
    samples = read_args(project_name, config)

    # Parametros de configuración de este script
    PROJECTS_PATH = config['PROJECTS_PATH']
    SNIPPY_PATH = config['SNIPPY_PATH']

    # Aqui puedes añadir opciones a trimomatic http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
    SNIPPY_OPTIONS = config['SNIPPY_OPTIONS']
    SNIPPY_REFERENCE = config['REFERENCE']
    if not os.path.exists(SNIPPY_REFERENCE):
        logger.error("The reference file does not exist check snippy.json file")
        logger.error(f"This file does not exist:{SNIPPY_REFERENCE}")
        return False

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

        PROCESSED_RESULTS = os.path.join(OUTPUT_PATH, "processed")
        PROCESSED_FILE = os.path.join(PROCESSED_RESULTS, f"{sample_name}_snippy.csv")
        if only_output or (os.path.exists(PROCESSED_FILE) and not extra_config["force"]):
            if os.path.exists(VCF_FILE):
                result = True
            else:
                result = False
                logger.warning(f"Output already generated using {sample_name}_snippy.csv")
        else:
            # Ejecutar Trimmomatic
            command = [SNIPPY_PATH, "--outdir", OUTPUT_RESULTS, "--ref", SNIPPY_REFERENCE, "--R1", input_r1_path, "--R2", input_r2_path] + SNIPPY_OPTIONS
            result = execute_command(command)
        if result:
            logger.info("Snippy process for %s finished", sample_name)
            process_output(VCF_FILE, sample_name, OUTPUT_PATH)
    
    
    combined_excel_files(samples, OUTPUT_PATH)

    if not extra_config["keep_output"]:
        os.system(f"rm -r {OUTPUT_PATH}/output")

if __name__ == "__main__":

    # Define the arguments that the program expects
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--parse-output', action='store_true', help='Set the flag to not execute but only process output files')
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    parser.add_argument('--force', action='store_true', help='Force the execution of the program')
    parser.add_argument('--keep_output', action='store_true', help='Force the execution of the program')
    args = parser.parse_args()
    PROJECT_NAME = args.PROJECT_NAME
    if args.json_config:
        config = init_configs(script_directory, args.json_config, required_keys=["SNIPPY_PATH", "SNIPPY_OPTIONS", "REFERENCE", "POLYMORPHISMS"])

    configure_logs(PROJECT_NAME, "snippy", config)
    logger = logging.getLogger(__name__)
    snippy_run(PROJECT_NAME, args.parse_output, config, extra_config={"force": args.force, "keep_output": args.keep_output})
