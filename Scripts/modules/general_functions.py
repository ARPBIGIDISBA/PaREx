# Funciones para todos los scripts
import argparse
import subprocess
import logging
import json
import os

LOG_MODE = "w" # Si pones W sobre escribe en cada ejecución y con "a" añade al log uno

def read_args(config_json):
    ''' Función para leer los argumentos de la línea de comandos, el fichero de samples y la configuración del json'''
    
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('PROJECT_NAME', type=str, help='Nombre del projecto')

    args = parser.parse_args()

    PROJECT_NAME = args.PROJECT_NAME
    script_name = config_json.split("_")[0]
    # Leer el archivo catde configuración
    with open(config_json, 'r') as file:
        config = json.load(file)

    # Crear el directorio para el lineage
    PROJECT_PATH = os.path.join(config["PROJECTS_PATH"], PROJECT_NAME)
    os.makedirs(PROJECT_PATH, exist_ok=True)

    SAMPLE_FILE = os.path.join(PROJECT_PATH, f"SAMPLES_LIST_{PROJECT_NAME}")
    
    logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=[
                        logging.FileHandler(f'{PROJECT_NAME}_{script_name}.log', mode=LOG_MODE),
                        logging.StreamHandler()
                    ])
    
    try:
        with open(SAMPLE_FILE, 'r') as file:
            samples = file.readlines()

        return PROJECT_NAME, samples, config, logging
    except Exception:
    
        print("Error al leer el fichero de samples que has pasado como argumento")
        print(SAMPLE_FILE)
        exit(1)


def read_args_bowtie():
    # Definir los argumentos que el programa espera
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('SAMPLES', type=str, help='Path para el fichero con la lista de samples')
    parser.add_argument('LINEAGE', type=str, help='Lineage a procesar')
    parser.add_argument('GENE', type=str, help='GENE como librería de alineamiento')

    # Leer los argumentos de la línea de comandos
    args = parser.parse_args()

    SAMPLES = args.SAMPLES
    LINEAGE = args.LINEAGE
    GENE = args.GENE
    
    logging.basicConfig(level=logging.INFO, filename=f'{LINEAGE}.log',
                        filemode='w',
                        format='%(name)s - %(levelname)s - %(message)s')

    try:
        with open(SAMPLES, 'r') as file:
            lines = file.readlines()

        return SAMPLES, LINEAGE, GENE, lines, logging
    except Exception:
        print("Error al leer el fichero de samples que has pasado como argumento")
        exit(1)


# Execute a command and log the output
def execute_command(command, logging):
    # Enseñar el comando que se va a ejecutar
    logging.info(f"Ejecutando el comando: {' '.join(command)}")

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Capturar las salidas del comando y añadirlas al fichero de salida
    stdout, stderr = process.communicate()

    # Log the output
    if stdout:
        logging.info(stdout.decode())
    if stderr:
        logging.info(stderr.decode())

    # Return True if the process was successful
    return (process.returncode == 0)