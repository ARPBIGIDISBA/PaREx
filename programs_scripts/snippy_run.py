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
import re
import pandas as pd
import traceback

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
    sheet_names = [
        'Extended_resistome',
        'Basic_resistome',
        'Cefidorocol_resistome',
        'Hypermutation'
    ]
    output = {}

    for sheet in sheet_names:
        df = pd.read_excel(filename, sheet_name=sheet).fillna('')
        df = df.rename(columns=lambda x: x.strip())
        data = {}
        for _, row in df.iterrows():
            locus = row['GENE_LOCUS']
            gene = row['GENE_NAME']
            polymorphisms = row['POLYMORPHISMS']
            if gene == "-":
                gene = ""
            # Split and clean polymorphisms
            polymorphisms = polymorphisms.rstrip('.').split(',')
            polymorphisms = [p.strip() for p in polymorphisms if p.strip()]
            cleaned_polymorphisms = [p.split(' (')[0].strip() for p in polymorphisms]
            data[locus] = {'gene': gene, 'polymorphisms': cleaned_polymorphisms}
        output[sheet] = data

    return output


def process_output(vcf_path, sample_name, output_path):
    import pysam
    import os
    import csv

    output_dir = os.path.join(output_path, "processed")
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, f"{sample_name}_snippy.csv")
    if os.path.exists(csv_path):
        os.remove(csv_path)

    with open(csv_path, mode='w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file, delimiter=";")
        # Puedes añadir aquí el filtro si quieres: ['locus', 'genes', ..., 'filter_name']
        csv_writer.writerow(['locus', 'genes', 'P.', 'changes', 'filtered_mutations', 'C.'])

        vcf_file = pysam.VariantFile(vcf_path)
        path = os.path.join(config.get("DATABASE_PATH"), "snippy", config.get("POLYMORPHISMS"))
        filter_dict = read_data_from_file(path)

        # Construir mapping: locus -> filtro
        locus_to_filter = {}
        for filter_name, loci in filter_dict.items():
            for locus in loci.keys():
                locus_to_filter[locus] = filter_name

        results = []
        for record in vcf_file:
            if 'ANN' in record.info:
                annotations = record.info['ANN']
                for annotation in annotations:
                    fields = annotation.split('|')
                    if len(fields) > 1 and fields[2] not in ('LOW', 'MODIFIER'):
                        locus = fields[3]
                        if locus in locus_to_filter:
                            filtro_actual = locus_to_filter[locus]
                            gene_name = filter_dict[filtro_actual][locus]['gene']
                            known_polymorphisms = filter_dict[filtro_actual][locus]['polymorphisms']
                            p = fields[10]
                            c = fields[9]
                            changes = translate_amino_acid(p, c)
                            filter_mutations = [m for m in changes if m not in known_polymorphisms]
                            row = [
                                locus,
                                gene_name,
                                p.replace("p.", ""),
                                ",".join(changes),
                                ",".join(filter_mutations),
                                c.replace("c.", "")
                                # Si quieres añadir el filtro, pon: , filtro_actual
                            ]
                            results.append(row)
                            csv_writer.writerow(row)

def combined_excel_files(samples, output_path, generate_full_hyperresistome=False):
    # Set paths and files
    output_dir = os.path.join(output_path, "processed")
    csv_output = os.path.join(output_path, "combined_snippy.xlsx")
    csv_output_full = os.path.join(output_path, "combined_snippy_full.xlsx")
    path = os.path.join(config.get("DATABASE_PATH"), "snippy", config.get("POLYMORPHISMS"))

    # Define tabs to filter by 
    filter = read_data_from_file(path)
    filter_names = ['Extended_resistome', 'Basic_resistome', 'Cefidorocol_resistome', 'Hypermutation']
    filters = {name: filter[name] for name in filter_names}

    columns_per_filter = {
        "Extended_resistome": ["locus", "genes", "changes", "filtered_mutations"],
        "Basic_resistome": ["locus", "genes", "changes", "filtered_mutations"],
        "Cefidorocol_resistome": ["locus", "genes", "filtered_mutations"],
        "Hypermutation": ["locus", "genes", "changes", "filtered_mutations"]
    }
    
    dataframes = {}
    for filter_name, loci in filters.items():
        columns = ["sample_name"] + [
            f"{locus}_{info['gene']}" for locus, info in loci.items()
        ]
        dataframes[filter_name] = pd.DataFrame(columns=columns).set_index("sample_name")
        dataframes[f"{filter_name}_clean"] = pd.DataFrame(columns=columns).set_index("sample_name")

    # Initialize dataframes with empty values
    for sample_name in samples:
        sample_name = sample_name.strip()
        csv_path = os.path.join(output_dir, f"{sample_name}_snippy.csv")
        if not os.path.exists(csv_path):
            logger.warning("The file %s does not exist", csv_path)
            continue
        df = pd.read_csv(csv_path, delimiter=";", na_values=["nan", "NaN"])
        df.replace("nan", pd.NA, inplace=True)
        df.replace("NaN", pd.NA, inplace=True)
        df = df.fillna('')
        

        for filter_name, loci in filters.items():
            columns = columns_per_filter[filter_name]
            for _, row in df.iterrows():
                locus = row['locus']
                gene = row['genes']
                changes = row.get('changes', '')
                filtered_mutations = row.get('filtered_mutations', '')
                colname = f"{locus}_{gene}"
                if locus in loci and colname in dataframes[filter_name].columns:
                    # cambios
                    prev = dataframes[filter_name].at[sample_name, colname] if sample_name in dataframes[filter_name].index else ''
                    parts = [str(x) for x in [changes, prev] if not (pd.isna(x) or x in ['', 'nan', 'NaN'])]
                    val = ",".join(parts)
                    dataframes[filter_name].at[sample_name, colname] = val

                    # clean
                    prev_clean = dataframes[f"{filter_name}_clean"].at[sample_name, colname] if sample_name in dataframes[f"{filter_name}_clean"].index else ''
                    parts_clean = [str(x) for x in [filtered_mutations, prev_clean] if not (pd.isna(x) or x in ['', 'nan', 'NaN'])]
                    val_clean = ",".join(parts_clean)
                    dataframes[f"{filter_name}_clean"].at[sample_name, colname] = val_clean
                
    if os.path.exists(csv_output):
        os.remove(csv_output)
    
    with pd.ExcelWriter(csv_output, engine='openpyxl') as writer:
        dataframes['Extended_resistome'].to_excel(writer, sheet_name='Extended_resistome', index=True)
        dataframes['Extended_resistome_clean'].to_excel(writer, sheet_name='Extended_resistome_clean', index=True)
        dataframes['Basic_resistome'].to_excel(writer, sheet_name='Basic_resistome', index=True)
        dataframes['Basic_resistome_clean'].to_excel(writer, sheet_name='Basic_resistome_clean', index=True)
    logger.info("Saving the file %s", csv_output)  

    if generate_full_hyperresistome:
        if os.path.exists(csv_output_full):
            os.remove(csv_output_full)

        with pd.ExcelWriter(csv_output_full, engine='openpyxl') as writer:
            for name, df in dataframes.items():
                df.to_excel(writer, sheet_name=name, index=True)
        logger.info("Saving the file full hyperresistome %s", csv_output_full)  
        

def snippy_run(project_name, only_output=False,  config=config, extra_config={"force": False, "keep_output": True, "file": None, "add_full_hyperresistome": False}):
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

    logger.info("Samples to process %s", samples)


    # Parametros de configuración de este script
    PROJECTS_PATH = config['PROJECTS_PATH']
    gi_PATH = config['SNIPPY_PATH']

    # Aqui puedes añadir opciones a trimomatic http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
    SNIPPY_OPTIONS = config['SNIPPY_OPTIONS']
    SNIPPY_REFERENCE = os.path.join(config['DATABASE_PATH'], 'snippy', config['REFERENCE'])
    
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
            input_r1_path = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}", f"{sample_name}_R1_001.fastq")
            input_r2_path = os.path.join(PROJECT_PATH, f"FASTQ_{project_name}", f"{sample_name}_R2_001.fastq")
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

    combined_excel_files(samples, OUTPUT_PATH, generate_full_hyperresistome=extra_config["add_full_hyperresistome"])

    if not extra_config["keep_output"]:
        os.system(f"rm -r {OUTPUT_PATH}/output")

if __name__ == "__main__":

    # Define the arguments that the program expects
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--parse-output', action='store_true', help='Set the flag to not execute but only process output files')
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)
    parser.add_argument('--force', action='store_true', help='Force the execution of the program')
    parser.add_argument('--keep-output', action='store_true', help='Force the execution of the program')
    parser.add_argument('--file', type=str, help='Direct file to process', default=None)
    parser.add_argument('--add-full-hyperresistome', action='store_true', help='Add the full hyperresistome to the output')
    args = parser.parse_args()
    PROJECT_NAME = args.PROJECT_NAME
    if args.json_config:
        config = init_configs(script_directory, args.json_config, required_keys=["SNIPPY_PATH", "SNIPPY_OPTIONS", "REFERENCE", "POLYMORPHISMS"])

    configure_logs(PROJECT_NAME, "snippy", config)
    logger = logging.getLogger(__name__)
    snippy_run(PROJECT_NAME, args.parse_output, config, extra_config={"force": args.force, "keep_output": args.keep_output, "file": args.file, "add_full_hyperresistome": args.add_full_hyperresistome})
