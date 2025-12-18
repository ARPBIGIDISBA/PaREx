'''
This software is licensed under the Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0).
More details: https://creativecommons.org/licenses/by-nc/4.0/

This script is used to generate an excel file with the results of the pipeline
'''

import os
import argparse
import logging
from modules.general_functions import configure_logs, init_configs
import glob
import pandas as pd
import datetime

logger = logging.getLogger(__name__)
script_path = os.path.abspath(__file__)
script_directory = os.path.dirname(script_path)
config = init_configs(script_directory)


def read_csv_results(csv_path, columns,  sample_id_col="sample_name", delimiter=";"):
    ''' 
        Reads a CSV file and returns a DataFrame with the specified columns.
        If columns is ["all"], returns all columns.
        The sample_id_col is renamed to "STRAIN ID" and set as index.
        If the file does not exist or there is an error, returns None.
        The delimiter is used to read the CSV file.
    '''
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
    '''
        Processes all the resfinder results in the given path and returns a DataFrame with the results.
        The sample_id_col is used to identify the sample name in the resfinder files.
        The function looks for files ending with "*fullcoverage.csv" in the given path.
        The resulting DataFrame has the following columns:
            - STRAIN ID
            - beta
            - aminoglycoside
            - fluoroquinolones
            - other
        Each column contains a list of genes found in that category, separated by commas.
        If no results are found, returns None.
    '''
    # Create an empty list to store the rows of the DataFrame
    rows = []
    for file in glob.glob(os.path.join(resfinder_path, "*fullcoverage.csv")):
        sample_name = os.path.basename(file).replace(".fullcoverage.csv", "")
        # Read the file as a DataFrame
        df = pd.read_csv(file, delimiter=";")
        # Ensure that the `sample_id_col` and `phenotypes` columns exist
        if sample_id_col in df.columns and 'phenotypes' in df.columns:
            df["phenotypes"] = df["phenotypes"].fillna("").str.split(",")  # Split into lists

            # Initialize the data structure for this sample
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

            def check_if_exist(row, existing_list, existing_list_name=""):
                '''
                    Check if a gene already exists in the existing_list based on
                    query_start_pos, query_end_pos and identity.
                    If it exists, merge the names and return True.
                    If not, return False.
                '''
                # round to 2 decimals
                identity_row = float(row['identity'].replace(',', '.'))
                identity_row = f"{identity_row:.2f}"

                for idx, sample in enumerate(existing_list):    
                    identity_sample = float(sample['identity'].replace(',', '.'))
                    identity_sample = f"{identity_sample:.2f}"
                    logger.debug(f"Checking gene: {row['name']} ({identity_row}%) against existing gene: {sample['name']} ({identity_sample}%)")

                    if (
                        row['query_start_pos'] == sample['query_start_pos'] and
                        row['query_end_pos'] == sample['query_end_pos'] and
                        identity_row == identity_sample):

                        logger.debug(f"Duplicate gene found: {row['name']} with positions {row['query_start_pos']}-{row['query_end_pos']} and identity {row['identity']}")
                        removed_sample = existing_list.pop(idx) 
                        gene_existing = removed_sample["name"]
                        gene_new = row["name"]
                        
                        if gene_new not in gene_existing.split(" / "):
                            logger.debug(f"Merging gene names for sample {removed_sample['name']} and {row['name']}")
                            new_name = f"{gene_existing} / {gene_new}"
                            # Modify de existing_list_name to include both genes
                            for element in existing_list_name:
                                if element == f"{gene_existing} ({identity_row}%)":
                                    logger.debug(f"Removing {element} from existing_list_name")
                                    existing_list_name.remove(element)
                                    existing_list_name.append(f"{new_name} ({identity_row}%)")
                                    break

                            removed_sample["name"] = new_name
                            
                        existing_list.append(removed_sample)
                        logger.debug("Merged sample:", existing_list_name)
                        
                        # Devolvemos False e inmediatamente salimos (comportamiento original)
                        return existing_list, existing_list_name
                        
                # Si el bucle termina sin coincidencia, retornamos True
                existing_list_name.append(f"{row['name']} ({identity_row}%)")
                existing_list.append(row)
                return existing_list, existing_list_name

            # Classify genes by categories
            for _, row in df.iterrows():
                gene = row[sample_id_col]
                phenotypes = [p.strip() for p in row["phenotypes"]]  # Remove whitespace
                
                # Categorization based on phenotypes and gene names
                # Avoid duplicates based on gene name and positions
                if any(phenotype in ["tobramycin", "gentamycin", "amikacin", "aph", "aad"] for phenotype in phenotypes):
                    added_aminoglycoside, sample_data["aminoglycoside"] = check_if_exist(row, added_aminoglycoside, sample_data["aminoglycoside"])
                elif any(phenotype in ["fluoroquinolones", "ciprofloxacin"] for phenotype in phenotypes):
                    added_fluoroquinolones, sample_data["fluoroquinolones"] = check_if_exist(row, added_fluoroquinolones, sample_data["fluoroquinolones"])
                elif gene.startswith("bla"):
                    added_beta, sample_data["beta"] = check_if_exist(row, added_beta, sample_data["beta"])
                    
                else:
                    added_other, sample_data["other"] = check_if_exist(row, added_other, sample_data["other"])

            # Add the row to the list of rows
            rows.append(sample_data)

    # If no results were found, return None
    if not rows:
        return None

    # Create the DataFrame directly from the list of rows
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

        # Remove PDC columns if they already exist
        for col in moving_cols:
            if col in cols:
                cols.remove(col)

        # Find the position of the reference column
        if ref_col in cols:
            idx = cols.index(ref_col)
            # Insert the columns after, maintaining order
            for col in reversed(moving_cols):  # reversed to keep the correct order
                cols.insert(idx + 1, col)

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
        "aminoglycoside": "Acquired aminoglycoside modifying genes",
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

    primeras = ["ST", "MLST allelic profile", "Acquired beta-lactamases", "Acquired aminoglycoside modifying genes", "Acquired quinolones resistance genes", "Other acquired resistance genes"]  # las que quieras primero
    
    primeras_ok = [c for c in primeras if c in df.columns]
    resto = [col for col in df.columns if col not in primeras_ok]
    
    df = df.reindex(columns=primeras_ok + resto)      

    # Sort PDC after P4A110_ampc
    pdc_cols = ["aminoacid substitutions (vs PDC-1)", "PDC variant (RefSeq protein ID)"]
    ref_col = "PA4110_ampC"
    df = reorder_columns(df, pdc_cols, ref_col)

    ref_col = "_oprD"
    oprD_cols = ["oprD", "oprD_reference-strain"]
    df = reorder_columns(df, oprD_cols, ref_col)

    # Reorder oprD oprD_reference-strain to be closest to PA0958
    # First remove column from _oprD file
    if "_oprD" in df.columns:
        df = df.drop(columns=["_oprD"])

    return df

def add_gene_absence_results(GENE_ABSENCE_PATH, combined_df, snippy_run=False):
    '''
        Adds gene absence results to the combined DataFrame.
        If the gene absence results already exist in the combined DataFrame, they are replaced.
    '''
    
    genes_absence_samples = read_csv_results(GENE_ABSENCE_PATH, ["all"])
    if genes_absence_samples is None:
        return combined_df
    
    columns_mapping = {
        "PA2020": "PA2020_mexZ",
        "PA2019": "PA2019_mexX",
        "PA2018": "PA2018_mexY",
        "PA2023": "PA2023_galU",
        "PA4522": "PA4522_ampD",
    }
    
    genes_absence_samples = genes_absence_samples.rename(columns=columns_mapping)
    
    if genes_absence_samples is not None:
        # For each STRAIN ID in combined_df, check if it exists in genes_absence_samples and add the columns or replace them if exist
        for strain_id in combined_df.index:
            if strain_id in genes_absence_samples.index:
                for col in genes_absence_samples.columns:
                    if col not in combined_df.columns:
                        combined_df[col] = ""  # crea la columna si falta
                        
                    if combined_df.at[strain_id, col] and pd.notna(combined_df.at[strain_id, col]) and combined_df.at[strain_id, col] != "" \
                        and genes_absence_samples.at[strain_id, col] and pd.notna(genes_absence_samples.at[strain_id, col]) and genes_absence_samples.at[strain_id, col] != "":
                        combined_df.at[strain_id, col] = f"{genes_absence_samples.at[strain_id, col]} ({combined_df.at[strain_id, col]})"  
                    else:
                        combined_df.at[strain_id, col] = genes_absence_samples.at[strain_id, col]  

    return combined_df

def add_piuAD_results(PIUAD_PATH, combined_df, snippy_run=False):
    ''' Add column piuA/D and piuA/D_REFERENCE to the combined DataFrame from the piuAD results file.
        If the piuAD results already exist in the combined DataFrame, they are replaced.
        Position in the DataFrame: after PA4514. search in combineted_df columns for PA4514 and insert 
        before it the columns piuA/D and piuA/D_REFERENCE
    '''
    
    piuAD_columns = ["piuA/D", "piuA/D_REFERENCE"]
    piuAD_samples = read_csv_results(PIUAD_PATH, piuAD_columns)
    
    if piuAD_samples is not None:
        
        #Insert columns in pandas DataFrame if exist after PA4514_piuA
        df_cols = list(combined_df.columns)
        if "PA4514_piuA" in df_cols:
            pa4514_index = df_cols.index("PA4514_piuA")
            df_cols.remove("PA4514_piuA")
        else:
            if snippy_run:
                return combined_df  # If snippy run and PA4514_piuA not found, skip
            pa4514_index = len(df_cols)  # If not found, append at the end

        for col in reversed(piuAD_columns):
            if col not in df_cols:
                df_cols.insert(pa4514_index, col)
        
        combined_df = combined_df.reindex(columns=df_cols)
    
        # For each STRAIN ID in combined_df, check if it exists in piuAD_samples and add the columns or replace them if exist
        for strain_id in combined_df.index:
            if strain_id in piuAD_samples.index:
                for col in piuAD_samples.columns:
                    if combined_df.at[strain_id, col] and pd.notna(combined_df.at[strain_id, col]) and combined_df.at[strain_id, col] != "" \
                        and piuAD_samples.at[strain_id, col] and pd.notna(piuAD_samples.at[strain_id, col]) and piuAD_samples.at[strain_id, col] != "":
                        combined_df.at[strain_id, col] = f"{piuAD_samples.at[strain_id, col]} ({combined_df.at[strain_id, col]})"
                    else:
                        combined_df.at[strain_id, col] = piuAD_samples.at[strain_id, col]


    return combined_df


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
    
    GENE_ABSENCE_PATH = os.path.join(OUTPUT_PATH, "gene_absence_results", f"{project_name}_gene_absence_results.csv")
    PIUAD_PATH = os.path.join(OUTPUT_PATH, "piuAD_results", f"{project_name}_piuAD_results.csv")

    # SNIPPY Section
    snippy_csv = os.path.join(OUTPUT_PATH, "snippy_results", f"combined_snippy.xlsx")
    if os.path.exists(snippy_csv):
        sheets = ["Extended_resistome", "Extended_resistome_clean", "Basic_resistome", "Basic_resistome_clean"]
        dfs = {}
        for sheet in sheets:
            df_snippy = pd.read_excel(snippy_csv, sheet_name=sheet)
            df_snippy.rename(columns={"sample_name": "STRAIN ID"}, inplace=True)
            df_snippy['STRAIN ID'] = df_snippy['STRAIN ID'].astype(str).str.strip().str.replace('"', '')
            df_snippy = df_snippy.set_index("STRAIN ID")
            df_snippy = pd.concat([combined_df, df_snippy], axis=1)
            df_snippy = add_gene_absence_results(GENE_ABSENCE_PATH, df_snippy, snippy_run=True)
            df_snippy = add_piuAD_results(PIUAD_PATH, df_snippy, snippy_run=True)
            dfs[sheet] = df_snippy

        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            for sheet in dfs:
                df = rename_columns(dfs[sheet])
                df.to_excel(writer, sheet_name=sheet, index=True)

    else:
        ## Add to combined_df the index of the samples
        combined_df = rename_columns(combined_df)
        combined_df = add_gene_absence_results(GENE_ABSENCE_PATH, combined_df, snippy_run=False)
        combined_df = add_piuAD_results(PIUAD_PATH, combined_df, snippy_run=False)
        combined_df.to_excel(output_file, sheet_name="Summary", index=True)

    logger.info(f"Results written to {output_file}")


def generate_pdf_from_excel(project_name, config=config, extra_config=None):
    '''
        Generates a PDF file for each isolate in the excel file.
        The PDF file contains a table with the results of the isolate.
        The PDF file is saved in the OUTPUT_PATH/PDF_results/ folder.
    '''

    import os
    import pandas as pd
    from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, PageBreak, Paragraph, Spacer, Image
    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import cm
    from reportlab.lib import colors

    PROJECTS_PATH = config["PROJECTS_PATH"]
    LOGO_ARPBIG = "./programs_scripts/assets/logoarpbig.png"
    LOGO_PAREX = "./programs_scripts/assets/logoparex.png"
    OUTPUT_PATH = os.path.join(PROJECTS_PATH, project_name, f"ANALYSIS_{project_name}", "")
    excel_file = os.path.join(OUTPUT_PATH, f"{project_name}_summary.xlsx")
    PDF_OUTPUT_PATH = os.path.join(OUTPUT_PATH, "PDF_results", "")
    if not os.path.exists(PDF_OUTPUT_PATH):
        os.makedirs(PDF_OUTPUT_PATH)
    else:
        # Si existe, eliminamos los PDFs anteriores
        for file in os.listdir(PDF_OUTPUT_PATH):
            if file.endswith(".pdf"):
                os.remove(os.path.join(PDF_OUTPUT_PATH, file))

    # Check if sheet_name Basic_resistome_clean exists
    xls = pd.ExcelFile(excel_file)
    tab_name = "Basic_resistome_clean"
    if tab_name not in xls.sheet_names:
        logger.warning(f"Sheet 'Basic_resistome_clean' not found in the excel file {excel_file}. Please generate the excel file with snippy results first.")
        tab_name="Summary"

    df = pd.read_excel(excel_file, sheet_name=tab_name, index_col="Isolate ID")

    styles = getSampleStyleSheet()
    title_style = styles['Title']
    subtitle_style = styles['Heading2']
    section_style = styles['Heading3']
    normal_style = styles['BodyText']

    # Color #234356 #26878b
    small_style = ParagraphStyle(
        name="SmallText",
        parent=normal_style,
        fontSize=10,
    )
        
    class PDFDocTemplate(SimpleDocTemplate):

        def afterPage(self):
            # Logo header and footer for each page
            self.canv.saveState()
            # LOGO ARPBIG (left)
            self.canv.drawImage(LOGO_PAREX, cm-1.5*cm, A4[1] - 3.5*cm, width=9*cm, height=4*cm, preserveAspectRatio=True, mask='auto')
            
            # Título centrado
            x = A4[0]/2 - 80  # Ajusta para centrar todo el texto
            y = A4[1] - 1.3*cm
            font_size = 16
            self.canv.setFont("Times-BoldItalic", font_size)
            self.canv.drawString(x, y, "Pseudomonas aeruginosa")
            width_bold = self.canv.stringWidth("Pseudomonas aeruginosa", "Times-Bold", font_size)
            self.canv.setFont("Times-Bold", font_size)
            self.canv.drawString(x + width_bold, y, " Resistome Explorer")

            # LOGO ARPBIG (bottom left)
            self.canv.drawImage(
                    LOGO_ARPBIG,
                    x=1*cm,              # Left margin
                    y=2*cm,              # Bottom margin
                    width=5.5*cm,
                    height=2.6*cm,
                    preserveAspectRatio=True,
                    mask='auto'
                )
            # Footer
            self.canv.setFont('Helvetica', 8)
            self.canv.drawRightString(A4[0] - cm, 1.2 * cm, f"Page {self.page}")
            self.canv.setFont('Helvetica', 8)
            self.canv.drawRightString(6*cm, 1.2 * cm, f"Generated at {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
            self.canv.restoreState()

    def get_isolate_table(row):
        # Replace all the NaN values with empty strings
        row = row.fillna("")
        pdc_variant = row.get("PDC variant (RefSeq protein ID)", "")
        aminoacid_substitutions = row.get("aminoacid substitutions (vs PDC-1)", "")
                
        if aminoacid_substitutions and isinstance(aminoacid_substitutions, str):
            aminoacid_substitutions = aminoacid_substitutions.replace(",", ",")
            
        # Iterate columns starting with PA example PA4225_pchF  get the part after _ type it in cursiva y entre parentesis el valor, solo quiero una columna
        mutational_results = {}
        for col in row.index:
            if col.startswith("PA") and col not in ["PA4110_ampC", "PA0958"]:
                value = row.get(col, "")
                if value != "":
                    # Get the part after _ and put it in cursive
                    gene_name = col.split("_")[-1]
                    if gene_name == "":
                        gene_name = col.split("_")[-2]
                    # Key is lower case of gene_name
                    mutational_results[gene_name.lower()] = f"<i>{gene_name}</i> ({str(value)})"
                    
            if col == "oprD_reference-strain":
                # If it is oprD_reference-strain, it is added
                continue
            if col == "piuA/D_REFERENCE":
                continue
            if col == "piuA/D":
                mutational_results["piuA/D"] = f"<i>{col}</i> ({row.get('piuA/D_REFERENCE', '')})"
            if col=="oprD":
                value = row.get(col, "")
                if value != "WT" and value != "":
                    # If it is oprD and not WT, we add it at the end
                    oprD_reference = row.get("oprD_reference-strain", "")
                    if oprD_reference:
                        value = f"{value} [{oprD_reference}]"
                    mutational_results["oprD"] = f"<i>oprD</i> ({value} )"
        
        # Order by key in mutational_results
        mutational_results = {k: v for k, v in sorted(mutational_results.items(), key=lambda item: item[0])}
        mutational_list = []
        for key, value in mutational_results.items():
           mutational_list.append(value)
  
        data = [
            # Section MLST
            [Paragraph("<b>Sequence Type</b>", normal_style), ""],
            [Paragraph("<b>ST</b>"), row.get("ST", "")],
            [Paragraph("<b>MLST allelic profile</b>"), Paragraph(row.get("MLST allelic profile", ""))],

            # Section PDC
            [Paragraph("<b>PDC</b>", normal_style), ""],
            [Paragraph("<b>Aminoacid substitutions (vs PDC-1)</b>"), Paragraph(aminoacid_substitutions.replace(",", ", "))],
            [Paragraph("<b>PDC variant (RefSeq protein ID)</b>"), Paragraph(pdc_variant)],

            # Section adquired resistome
            [Paragraph("<b>Horizontally acquired resistome</b>", normal_style), ""],
            [Paragraph("<b>Beta-lactamases</b>"), Paragraph(row.get("Acquired beta-lactamases", ""))],
            [Paragraph("<b>Aminoglycosides resistance genes</b>"), Paragraph(row.get("Acquired aminoglycoside modifying enzymes", ""))],
            [Paragraph("<b>Quinolones resistance genes</b>"), Paragraph(row.get("Acquired quinolones resistance genes", ""))],
            [Paragraph("<b>Other resistance genes</b>"), Paragraph(row.get("Other acquired resistance genes", ""))],

            # Section mutational resistome
            [Paragraph("<b>Mutational resistome</b>", normal_style), ""],
        ]


        if mutational_list:
            data.append([Paragraph(", ".join(mutational_list), small_style), ""])
        else:
            data += ["", ""]
                    
        table = Table(data, colWidths=[5.2*cm, 12*cm])

        COLORS_LIGHT_BLUE = colors.HexColor("#9bd1d3")
        # Merged cells y styles
        style = TableStyle([
            ('SPAN', (0,0), (1,0)),  # Sequence Type
            ('SPAN', (0,3), (1,3)),  # PDC
            ('SPAN', (0,6), (1,6)),  # Horizontally acquired resistome
            ('SPAN', (0,11), (1,11)), # Mutational resistome
            ('SPAN', (0,12), (1,12)), # Mutational value
            
            ('BACKGROUND', (0,0), (1,0), COLORS_LIGHT_BLUE),
            ('BACKGROUND', (0,3), (1,3), COLORS_LIGHT_BLUE),
            ('BACKGROUND', (0,6), (1,6), COLORS_LIGHT_BLUE),
            ('BACKGROUND', (0,11), (1,11), COLORS_LIGHT_BLUE),

            ('VALIGN', (0,0), (-1,-1), 'TOP'),
            ('FONTSIZE', (0,0), (-1,-1), 9),
            ('GRID', (0,0), (-1,-1), 0.25, colors.grey),
            ('LEFTPADDING', (0,0), (-1,-1), 8),
            ('RIGHTPADDING', (0,0), (-1,-1), 8),
            ('BOTTOMPADDING', (0,0), (-1,-1), 6),
            ('TOPPADDING', (0,0), (-1,-1), 6),
        ])
        table.setStyle(style)
        return table

    for _, row in df.iterrows():
        pdf_file = os.path.join(PDF_OUTPUT_PATH, f"{row.name}_summary.pdf")
        doc = PDFDocTemplate(pdf_file, pagesize=A4)
        elements = []
        elements.append(Spacer(1, 1.8*cm))
        elements.append(Paragraph(f"Isolate ID: {row.name}", styles['Heading2']))
        elements.append(Spacer(1, 0.2*cm))
        elements.append(get_isolate_table(row))
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

