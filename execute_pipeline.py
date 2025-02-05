import argparse
import logging
import sys
import os
import glob

# Add path Scripts to sys.path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), 'programs_scripts')))

from programs_scripts.modules.general_functions import read_config, check_project
from programs_scripts.trimmomatic_run import trimmomatic_run
from programs_scripts.SPAdes_run import SPAdes_run
from programs_scripts.bowtie_run import bowtie_run
from programs_scripts.resfinder_run import resfinder_run
from programs_scripts.oprD_run import oprD_run
from programs_scripts.mlst_run import mlst_run
from programs_scripts.snippy_run import snippy_run
from programs_scripts.PDC_run import PDC_run
from programs_scripts.generate_excel_run import generate_excel_run
from programs_scripts.novasec_run import novasec_run


logger = logging.getLogger(__name__)

if __name__ == "__main__":

    OPERATIONS_DEVELOPED = ["create_project", "create_sample_list", "generate_excel", "trimmomatic",
                             "SPAdes", "bowtie", "resfinder", "oprD", "mlst", 
                             "all_sequence", "snippy", "PDC", "novasec", "projects"]
    
    parser = argparse.ArgumentParser(description='Execute pipeline scripts.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')
    parser.add_argument('operation', type=str, help=f"Existing operations: {' | '.join(OPERATIONS_DEVELOPED)}")
    parser.add_argument('--log-level', type=str, help='Log levels DEBUG, INFO, WARNING, ERROR', default="INFO")
    parser.add_argument('--force', action='store_true', help='Force the execution of the program')
    parser.add_argument('--keep_output', action='store_true', help='Force the execution of the program')
    parser.add_argument('--file', type=str, help='Direct Path to the file not use sample list', default=None)
    args = parser.parse_args()


    extra_config = {
        "force": args.force,
        "keep_output": args.keep_output,
        "log_level": args.log_level,
        "file": args.file
    }

    PROJECT_NAME = args.PROJECT_NAME
    OPERATIONS = args.operation.split(",")
    
    general_config = os.path.join("programs_scripts", "configs", "general.json")
    config_general = read_config(general_config, ["PROJECTS_PATH"])
    PROJECTS_PATH = config_general["PROJECTS_PATH"]
    project_path = os.path.join(PROJECTS_PATH, PROJECT_NAME)
    
    if PROJECT_NAME != "list":
        check_project(project_path)
        logging.basicConfig(level=args.log_level,
                        format='%(asctime)s - %(levelname)s - %(message)s - %(filename)s:%(lineno)d',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        handlers=[
                            logging.FileHandler(
                                os.path.join(project_path, "Logs", f'execute_pipeline.log'),
                                mode="w"),
                            logging.StreamHandler()
                        ])

    logger.info("executing operations: '%s'", " ".join(OPERATIONS))
    for operation in OPERATIONS:
        logger.info("Executing operation %s", operation)
        if operation not in OPERATIONS_DEVELOPED:
            logger.warning(f"Operation not found {operation}")
            logger.info(f"Operations available: {OPERATIONS_DEVELOPED}")
            sys.exit()
        if operation == "create_project":
            # comand line question if you are sure to continue
            logger.info(f"Creating project {PROJECT_NAME} estructure")
            project_path = os.path.join(PROJECTS_PATH, PROJECT_NAME)
            os.makedirs(project_path, exist_ok=True)
            os.makedirs(os.path.join(project_path, f"FASTQ_{PROJECT_NAME}"), exist_ok=True)
            os.makedirs(os.path.join(project_path, f"ANALYSIS_{PROJECT_NAME}"), exist_ok=True)
            with open(os.path.join(project_path, f"SAMPLES_LIST_{PROJECT_NAME}"), 'w') as file:
                pass 
        elif operation == "projects":
            logger.info("Listing projects")
            # List folder name in PROJECTS_PATH
            projects = os.listdir(PROJECTS_PATH)
            for project in projects:
                print("project %s" % project)

        elif operation == "create_sample_list":
            file_name = f"SAMPLES_LIST_{PROJECT_NAME}"
            path = os.path.join(project_path, file_name)
            logger.info(f"Creating {file_name} file")
            logger.debug(f"File path: {path}")
            with open(path, 'w') as file:
                path_fasta = os.path.join(project_path, f"FASTQ_{PROJECT_NAME}")
                fastaq_files = glob.glob(f'{path_fasta}/*.fastq')
                fastaq_gz_files = glob.glob(f'{path_fasta}/*.fastq.gz')
                fastaq_files = fastaq_files + fastaq_gz_files
                if len(fastaq_files) == 0:
                    logger.warning(f"No files found in {path_fasta} Using denovo SPAdes path now")
                    path_fasta = os.path.join(project_path, f"ANALYSIS_{PROJECT_NAME}", "denovo_assemblies_SPAdes")
                    logger.debug(f'{path_fasta}/*.fasta')
                    fastaq_files = glob.glob(f'{path_fasta}/*.fasta')
                    fastaq_files = [os.path.basename(f).replace(".SPAdes.denovoassembly.fasta", "") for f in fastaq_files]

                else:
                    fastaq_files = [os.path.basename(f) for f in fastaq_files]
                    fastaq_files = [f.split("_R1")[0] for f in fastaq_files]
                    fastaq_files = [f.split("_R2")[0] for f in fastaq_files]
                    fastaq_files = list(set(fastaq_files))
                # Print out all filenames
                logger.debug("Files found:")
                for filename in fastaq_files:
                    logger.debug(filename)
                    file.write(filename + '\n')
        
        elif operation == "generate_excel":
            logger.info(f"Running generate_excell for project {PROJECT_NAME}")
            generate_excel_run(PROJECT_NAME, extra_config=extra_config)
        elif operation == "trimmomatic":
            logger.info(f"Running trimmomatic for project {PROJECT_NAME}")
            trimmomatic_run(PROJECT_NAME, extra_config=extra_config)
        elif operation == "SPAdes":
            logger.info(f"Running SPAdes for project {PROJECT_NAME}")
            SPAdes_run(PROJECT_NAME, extra_config=extra_config)
        elif operation == "bowtie":
            logger.info(f"Running bowtie for project {PROJECT_NAME}")
            reference = args.reference
            bowtie_run(PROJECT_NAME, reference, extra_config=extra_config)
        elif operation == "resfinder":
            logger.info(f"Running resfinder for project {PROJECT_NAME}")
            resfinder_run(PROJECT_NAME, extra_config=extra_config)
        elif operation == "oprD":
            logger.info(f"Running oprD for project {PROJECT_NAME}")
            oprD_run(PROJECT_NAME, extra_config=extra_config)
        elif operation == "mlst":
            logger.info(f"Running mlst for project {PROJECT_NAME}")
            mlst_run(PROJECT_NAME, extra_config=extra_config)
        elif operation == "snippy":
            logger.info(f"Running snippy for project {PROJECT_NAME}")
            snippy_run(PROJECT_NAME, extra_config=extra_config)
        elif operation == "PDC":
            logger.info(f"Running PDC for project {PROJECT_NAME}")
            PDC_run(PROJECT_NAME, extra_config=extra_config)
        elif operation == "all_sequence":
            logger.info(f"Running all for project {PROJECT_NAME}")
            logger.info(f"short for SPAdes, resfinder, oprD, mlst, generate_excel")
            SPAdes_run(PROJECT_NAME, extra_config=extra_config)
            reference = args.reference
            resfinder_run(PROJECT_NAME, extra_config=extra_config)
            oprD_run(PROJECT_NAME, extra_config=extra_config)
            mlst_run(PROJECT_NAME, extra_config=extra_config)
            generate_excel_run(PROJECT_NAME, extra_config=extra_config)
        elif operation == "novasec":
            logger.info(f"Running novasec for project {PROJECT_NAME}")
            novasec_run(PROJECT_NAME, extra_config=extra_config)
        else:
            logger.warning("Operation not found %s", operation)
            logger.info("Operations available: create_project")
    
    
    
    # reference = args.REFERENCE
    
    # makedirs(path.join(PROJECTS_PATH, "Logs"), exist_ok=True)


    # # Execute trimmomatic process
    # trimmomatic_run(PROJECT_NAME)

    # # # Execute SPAdes analisis
    # SPAdes_run(PROJECT_NAME)

    # # # Execute bowtie analisis
    # bowtie_run(PROJECT_NAME, reference)

    # # # Execute resfinder analisis
    # resfinder_run(PROJECT_NAME)

    # # # Execute oprD analisis
    # oprD_run(PROJECT_NAME)
