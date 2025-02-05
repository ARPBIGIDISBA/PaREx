'''
    Este script aplica el ensamblaje de novo con SPAdes a los ficheros fastq.gz
    Ejecuta el programa en python spades sobre los ficheros fastq
    Tiene como entrada los ficheros fastq.gz de las muestras
    Da como resultado un fichero fasta con los SPAdes.denovoassembly.fasta
    Here are the command options for spades https://github.com/ablab/spades#sec3.2

'''
import os
import argparse
import logging
import csv
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs
import glob
import pandas as pd

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory)


def read_csv_results(csv_path, columns,  sample_id_col="sample_name", delimiter=";"):
    if not os.path.exists(csv_path):
        logger.error(f"File not found {csv_path}")
        logger.warning("Execute first the analysis")
        return None

    try:
        df = pd.read_csv(csv_path, delimiter=delimiter)
        if sample_id_col not in df.columns:
            logger.error(f"Column '{sample_id_col}' not found in the file.")
            return None
        df.rename(columns={sample_id_col: "STRAIN ID"}, inplace=True)
        df =  df.set_index("STRAIN ID")
        
        if columns == ["all"]:
            return df
        else:
            return df[columns]
    except Exception as e:
        logger.error(f"Error reading CSV file: {e}")
        logger.warning(f"Check the file format {csv_path}")
        return None


def process_resfinder_samples(resfinder_path, sample_id_col="name"):
    # Crear una lista vacía para almacenar las filas del DataFrame
    rows = []

    for file in glob.glob(os.path.join(resfinder_path, "*fullcoverage.csv")):
        sample_name = os.path.basename(file).replace(".fullcoverage.csv", "")
        
        # Lee el archivo como DataFrame
        df = pd.read_csv(file, delimiter=";")
        # Asegúrate de que la columna `sample_id_col` y `phenotypes` existen
        if sample_id_col in df.columns and 'phenotypes' in df.columns:
            df["phenotypes"] = df["phenotypes"].fillna("").str.split(",")  # Divide en listas
            
            # Inicializa la estructura de datos para este sample
            sample_data = {
                "STRAIN ID": sample_name,
                "beta": [],
                "aminoglycoside": [],
                "fluoroquinolones": [],
                "other": []
            }
            
            # Clasifica los genes por categorías
            for _, row in df.iterrows():
                gene = row[sample_id_col]
                phenotypes = [p.strip() for p in row["phenotypes"]]  # Elimina espacios en blanco
                # Categorización de acuerdo a los fenotipos
                if any(phenotype in ["tobramycin", "gentamycin", "amikacin", "aph", "aad"] for phenotype in phenotypes):
                    sample_data["aminoglycoside"].append(gene)
                elif any(phenotype in ["fluoroquinolones", "ciprofloxacin"] for phenotype in phenotypes):
                    sample_data["fluoroquinolones"].append(gene)
                elif gene.startswith("bla"):
                    identity = row["identity"]
                    sample_data["beta"].append(f"{gene} ({identity}%)")
                else:
                    sample_data["other"].append(gene)

            # Añadir la fila a la lista de filas
            rows.append(sample_data)

    # Si no se encontraron resultados, devolver None
    if not rows:
        return None

    # Crear el DataFrame directamente desde la lista de filas
    resfinder_df = pd.DataFrame(rows)
    resfinder_df.set_index("STRAIN ID", inplace=True)
    
    # Combina las listas en cada celda en una sola cadena separada por comas
    for col in resfinder_df.columns:
        resfinder_df[col] = resfinder_df[col].apply(lambda x: ", ".join(x) if isinstance(x, list) else x)

    return resfinder_df

def generate_excel_run(project_name, config=config, extra_config=None):
    PROJECTS_PATH = config["PROJECTS_PATH"]
    OUTPUT_PATH = os.path.join(PROJECTS_PATH, project_name, f"ANALYSIS_{project_name}", "")

    # Load data from CSVs
    csv_files = {
                    "oprD": {
                        "path": os.path.join(OUTPUT_PATH, "oprD_results", f"{project_name}_oprD_results.csv"),
                        "columns": ["oprD", "oprD_REFERENCE"]
                    },
                    "PDC": {
                        "path": os.path.join(OUTPUT_PATH, "PDC_results", f"{project_name}_PDC_results.csv"),
                        "columns": ["PDC", "PDC_REFERENCE"]
                    },
                    "mlst": {
                        "path": os.path.join(OUTPUT_PATH, "mlst_results", f"{project_name}_mlst_results.csv"),
                        "columns": ["sequence_type", "alleles" ]
                    }
                }

    sample_results = {name: read_csv_results(value["path"], value["columns"]) for name, value in csv_files.items()}

    resfinder_path = os.path.join(OUTPUT_PATH, "resfinder_results", "csv_samples")
    resfinder_samples = process_resfinder_samples(resfinder_path)
    sample_results['resfinder'] = resfinder_samples

    # Remove None values from sample_results    
    
    concat = []
    # Concatenate all the DataFrames by the index (STRAIN ID)
    for key, value in sample_results.items():
        if value is not None:
            concat.append(value)
    
    combined_df = pd.concat(concat, axis=1)
    
    # SNIPPY Section
    output_file = os.path.join(OUTPUT_PATH, f"{project_name}_summary.xlsx")

    snippy_csv = os.path.join(OUTPUT_PATH, "snippy_results", f"combined_snippy.xlsx")
    if os.path.exists(snippy_csv):
        sheets = ["All", "All_clean", "Basic", "Basic_clean"]
        for sheet in sheets:
            df_snippy = pd.read_excel(snippy_csv, sheet_name=sheet)
            ## rename sample_name to STRAIN ID16835951_S9_L001 16835951_S9_L001
            df_snippy.rename(columns={"sample_name": "STRAIN ID"}, inplace=True)
            df_snippy['STRAIN ID'] = df_snippy['STRAIN ID'].str.strip().str.replace('"', '')
            df_snippy = df_snippy.set_index("STRAIN ID")
            
            df_snippy = pd.concat([combined_df, df_snippy], axis=1)
            if sheet == "All":
                df_all = df_snippy
            elif sheet == "All_clean":
                df_all_clean = df_snippy
            elif sheet == "Basic":
                df_basic = df_snippy
            elif sheet == "Basic_clean":
                df_basic_clean = df_snippy

        # Cargar el archivo existente sin sobrescribir
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            df_all.to_excel(writer, sheet_name='All', index=True)
            df_all_clean.to_excel(writer, sheet_name='All_clean', index=True)
            df_basic.to_excel(writer, sheet_name='Basic', index=True)
            df_basic_clean.to_excel(writer, sheet_name='Basic_clean', index=True)

    else:
        combined_df.to_excel(output_file, sheet_name="Summary", index=True)
        ## Add to combined_df the index of the samples

        combined_df.to_excel(output_file, sheet_name="Summary", index=False)

    logger.info(f"Results written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('--log-level', type=str, help='Log levels DEBUG, INFO, WARNING, ERROR', default=None)
    parser.add_argument('--json-config', type=str, help='Json file in the config directory', default=None)

    args = parser.parse_args()
    project_name = args.PROJECT_NAME

    if args.json_config:
        config = init_configs(script_directory, f"{args.json_config}.json")

    # Start the python logging variable to generate a file
    configure_logs(project_name, "generate_excel", config, log_level=args.log_level)

    logger = logging.getLogger(__name__)

    generate_excel_run(project_name)

