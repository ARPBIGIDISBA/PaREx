import os

base_path = '/home/mbonet/microbiologia/Reference4alignment/PDC/protein/'
file = '/home/mbonet/microbiologia/Reference4alignment/PDC/PDCXX.fasta'
def split_fasta_file(file_path):
    try:
        with open(file_path, 'r') as file:
            current_file = None
            for line in file:
                if line.startswith('>'):  # Nueva secuencia encontrada
                    if current_file:
                        current_file.close()  # Cerrar el archivo actual antes de abrir uno nuevo
                    
                    name = line.split("[")[0]
                    identifier = name[1:].split(" ") # Extraer el identificador (eliminando '>')        
                    file_name = "_".join(identifier[:-1])  # Asignar el nombre del archivo basado en el identificador
                    file_name = file_name.replace(",", "")  # Reemplazar '/' por '_' en el nombre del archivo
                    current_file = open(os.path.join(base_path,f"{file_name}.fasta"), 'w')  # Abrir nuevo archivo para escritura
                    
                if current_file:
                    current_file.write(line)  # Escribir línea en el archivo actual

            if current_file:
                current_file.close()  # Asegurarse de cerrar el último archivo abierto
    except IOError as e:
        print(f"Error al abrir o leer el archivo: {e}")

# Llamar a la función con la ruta del archivo
split_fasta_file(file)