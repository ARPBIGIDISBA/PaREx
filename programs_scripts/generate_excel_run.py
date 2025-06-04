'''
This software is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).
More details: https://creativecommons.org/licenses/by-nc/4.0/

This script is used to generate an excel file with the results of the pipeline
'''

import os
import argparse
import logging
import csv
from modules.general_functions import read_args, execute_command
from modules.general_functions import configure_logs, init_configs
import glob
import pandas as pd
import re

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory)


def read_csv_results(csv_path, columns,  sample_id_col="sample_name", delimiter=";"):
    if not os.path.exists(csv_path):
        logger.debug(f"File not found {csv_path}")
        logger.debug("Execute first the analysis")
        return None

    try:
        df = pd.read_csv(csv_path, delimiter=delimiter)
        if sample_id_col not in df.columns:
            logger.error(f"Column '{sample_id_col}' not found in the file.")
            return None
        df.rename(columns={sample_id_col: "STRAIN ID"}, inplace=True)
        df = df.set_index("STRAIN ID")
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
            
            added_aminoglycoside = []
            added_fluoroquinolones = []
            added_beta = []
            added_other = []

            def check_if_exist(row, list):
                for sample in list:

                    if row["name"]==sample["name"] and row["query_start_pos"]==sample["query_start_pos"] and row["query_end_pos"]==sample["query_end_pos"]:
                        return False
                    
                    if sample["name"].startswith("blaOXA") and row["name"].startswith("blaOXA"):
                        return False
                return True

            # Clasifica los genes por categorías
            for _, row in df.iterrows():
                gene = row[sample_id_col]
                phenotypes = [p.strip() for p in row["phenotypes"]]  # Elimina espacios en blanco
                # Identity only 2 decimal places
                identity = f"{float(row['identity'].replace(',', '.')):.2f}" if row['identity'] is not None else ""
                
                # Categorización de acuerdo a los fenotipos
                if any(phenotype in ["tobramycin", "gentamycin", "amikacin", "aph", "aad"] for phenotype in phenotypes):
                    if check_if_exist(row, added_aminoglycoside):
                        added_aminoglycoside.append(row)
                        sample_data["aminoglycoside"].append(f"{gene} ({identity}%)")

                elif any(phenotype in ["fluoroquinolones", "ciprofloxacin"] for phenotype in phenotypes):
                    if check_if_exist(row, added_fluoroquinolones):
                        added_fluoroquinolones.append(row)
                        sample_data["fluoroquinolones"].append(f"{gene} ({identity}%)")
                elif gene.startswith("bla"):
                    if check_if_exist(row, added_beta):
                        added_beta.append(row)
                        sample_data["beta"].append(f"{gene} ({identity}%)")
                else:
                    if check_if_exist(row, added_other):
                        added_other.append(row)
                        sample_data["other"].append(f"{gene} ({identity}%)")

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

def reorder_columns(df, moving_cols, ref_col):
        """
        Reorders specific columns in the DataFrame to ensure they appear in a desired order.
        """
        cols = list(df.columns)

        # Quita las columnas PDC si ya están
        for col in moving_cols:
            if col in cols:
                cols.remove(col)

        # Busca la posición de la columna de referencia
        if ref_col in cols:
            idx = cols.index(ref_col)
            # Inserta las columnas después, manteniendo orden
            for col in reversed(moving_cols):  # reversed para que queden en el orden correcto
                cols.insert(idx + 1, col)

        # Reorder the columns
        df = df.reindex(columns=cols)
        return df

def rename_columns(df):
    """
    Renames columns in a DataFrame based on a mapping dictionary.
    
    Parameters:
    df (pd.DataFrame): The DataFrame whose columns are to be renamed.
    columns_mapping (dict): A dictionary where keys are current column names and values are new column names.
    
    Returns:
    pd.DataFrame: The DataFrame with renamed columns.
    """
    columns_mapping = {
        "sequence_type": "ST",
        "alleles": "MLST allelic profile",
        "beta": "Acquired beta-lactamases",
        "aminoglycoside": "Acquired aminoglycoside modifying enzymes",
        "fluoroquinolones": "Acquired quinolones resistance genes",
        "other": "Other acquired resistance genes",
        "oprD_REFERENCE": "oprD_reference-strain",
        "PDC": "aminoacid substitutions (vs PDC-1)",
        "PDC_REFERENCE": "PDC variant (RefSeq protein ID)"
    }
    # Strain ID is the index of the DataFrame, so we need to rename it separately
    if "STRAIN ID" in df.index.names:
        df.index.rename("Isolate ID", inplace=True)
    
    df = df.rename(columns=columns_mapping)

    primeras = ["ST", "MLST allelic profile", "Acquired beta-lactamases", "Acquired aminoglycoside modifying enzymes", "Acquired quinolones resistance genes", "Other acquired resistance genes"]  # las que quieras primero

    resto = [col for col in df.columns if col not in primeras]
    df = df.reindex(columns=primeras + resto)      

    # Ordenar PDC after P4A110_ampc 
    pdc_cols = ["aminoacid substitutions (vs PDC-1)", "PDC variant (RefSeq protein ID)"]
    ref_col = "PA4110_ampC"
    df = reorder_columns(df, pdc_cols, ref_col)

    # Reordenamos oprD	oprD_reference-strain lo mas cercano a PA0958
    reference_number = 958
    candidatas = []
    for col in df.columns:
        m = re.match(r"PA(\d+)", col)
        if m:
            num = int(m.group(1))
            candidatas.append((abs(num - reference_number), col))

    # Si hay candidatas, escoge la más cercana
    if candidatas:
        _, mejor_col = min(candidatas, key=lambda x: x[0])
        logger.debug(f"La columna más cercana a PA{reference_number} es: {mejor_col}")
    else:
        logger.debug("No se encontró ninguna columna PAxxxx")

    oprD_cols = ["oprD", "oprD_reference-strain"]
    ref_col = mejor_col if mejor_col else "PA0958"
    df = reorder_columns(df, oprD_cols, ref_col)

    return df

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
    output_file = os.path.join(OUTPUT_PATH, f"{project_name}_summary.xlsx")

    # SNIPPY Section
    snippy_csv = os.path.join(OUTPUT_PATH, "snippy_results", f"combined_snippy.xlsx")
    if os.path.exists(snippy_csv):
        sheets = ["Extended_resistome", "Basic_resistome"]
        # Add also the clean sheets
        sheets += [f"{sheet}_clean" for sheet in sheets]
        dfs = {}
        for sheet in sheets:
            df_snippy = pd.read_excel(snippy_csv, sheet_name=sheet)
            df_snippy.rename(columns={"sample_name": "STRAIN ID"}, inplace=True)
            df_snippy['STRAIN ID'] = df_snippy['STRAIN ID'].astype(str).str.strip().str.replace('"', '')
            df_snippy = df_snippy.set_index("STRAIN ID")
            df_snippy = pd.concat([combined_df, df_snippy], axis=1)
            dfs[sheet] = df_snippy
            
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            for sheet in dfs:
                df = rename_columns(dfs[sheet])
                df.to_excel(writer, sheet_name=sheet, index=True)

    else:
        ## Add to combined_df the index of the samples
        combined_df = rename_columns(combined_df)
        combined_df.to_excel(output_file, sheet_name="Summary", index=True)

    logger.info(f"Results written to {output_file}")


def generate_pdf_from_excel(project_name, config=config, extra_config=None):
    import os
    import pandas as pd
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, PageBreak, Paragraph, Spacer
    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.lib.units import cm

    PROJECTS_PATH = config["PROJECTS_PATH"]
    OUTPUT_PATH = os.path.join(PROJECTS_PATH, project_name, f"ANALYSIS_{project_name}", "")
    excel_file = os.path.join(OUTPUT_PATH, f"{project_name}_summary.xlsx")
    pdf_file = os.path.join(OUTPUT_PATH, f"{project_name}_summary.pdf")

    df = pd.read_excel(excel_file, sheet_name="All")

    styles = getSampleStyleSheet()
    style_title = styles['Heading2']

    class PDFDocTemplate(SimpleDocTemplate):
        def afterPage(self):
            self.canv.saveState()
            self.canv.setFont('Helvetica-Bold', 9)
            self.canv.drawString(cm, A4[1] - 1.5 * cm, "ARPBIG IDISBA")
            self.canv.setFont('Helvetica', 8)
            self.canv.drawRightString(A4[0] - cm, 1.2 * cm, f"Page {self.page}")
            self.canv.restoreState()

    doc = PDFDocTemplate(pdf_file, pagesize=A4)
    elements = []

    for _, row in df.iterrows():
        strain_id = row.get('STRAIN ID', 'Unknown')
        elements.append(Paragraph(f"Strain ID: {strain_id}", style_title))
        elements.append(Spacer(1, 6))

        data = [[k, str(v) if pd.notna(v) else ''] for k, v in row.items() if k != 'STRAIN ID']
        table = Table([['Field', 'Value']] + data, colWidths=[150, 380])
        table.setStyle(TableStyle([
            ('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey),  # Encabezado gris
            ('GRID', (0, 0), (-1, -1), 0.25, colors.black),
            ('FONTSIZE', (0, 0), (-1, -1), 8),
            ('VALIGN', (0, 0), (-1, -1), 'TOP'),
        ]))
        elements.append(table)
        elements.append(PageBreak())

    doc.build(elements)
    logger.info(f"PDF generated: {pdf_file}")
    return pdf_file



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

