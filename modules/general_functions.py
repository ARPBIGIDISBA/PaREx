import argparse
import subprocess
import logging
# Se escribe el log en un fichero con el nombre del lineage


# Funcion para leer dos parametros de la linea de comandos par luego leer un fichero con la lista de samples y lineage
def read_args():
    # Definir los argumentos que el programa espera
    parser = argparse.ArgumentParser(description='Procesa algunos argumentos.')
    parser.add_argument('SAMPLES', type=str, help='Path para el fichero con la lista de samples')
    parser.add_argument('LINEAGE', type=str, help='Lineage a procesar')

    # Leer los argumentos de la línea de comandos
    args = parser.parse_args()

    SAMPLES = args.SAMPLES
    LINEAGE = args.LINEAGE
    logging.basicConfig(level=logging.INFO, filename=f'{LINEAGE}.log', filemode='w', format='%(name)s - %(levelname)s - %(message)s')

    try:
        with open(SAMPLES, 'r') as file:
            lines = file.readlines()

        return SAMPLES, LINEAGE, lines, logging
    except Exception:
        print("Error al leer el fichero de samples que has pasado como argumento")
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


# Ejecutar un comando y loggear la salida
def execute_comman(command, logging):
    # Enseñar el comando que se va a ejecutar
    logging.info(f"Ejecutando el comando: {' '.join(command)}")

    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Capturar las salidas del comando y añadirlas al fichero de salida
    stdout, stderr = process.communicate()

    # Log the output
    if stdout:
        logging.info(stdout.decode())
    if stderr:
        logging.error(stderr.decode())