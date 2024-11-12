# Proyecto de Análisis de Resistencias Bacterianas

Este proyecto es un pipeline en Python diseñado para realizar análisis genómicos y de resistencia antimicrobiana en muestras de *Pseudomonas aeruginosa*. Incluye el procesamiento de archivos FASTQ, ensamblado de novo, alineación, identificación de fenotipos de resistencia y generación de reportes.

## Estructura del Proyecto

- **programs_scripts/**: Contiene scripts específicos para cada paso del pipeline.
- **configs/**: Archivos de configuración en formato `.json` necesarios para la ejecución de los scripts. También incluye archivos de ejemplo `.json.sample`.
- **logs/**: Carpeta donde se almacenan los logs generados durante la ejecución.

## Requisitos

- Python 3.8 o superior
- Paquetes de Python:
  - `pandas`
  - `openpyxl`
  - `pysam`
  - Otros paquetes listados en `requirements.txt`

## Instalación

1. Clona el repositorio:

   ```bash
   git clone https://github.com/tuusuario/proyecto.git
   cd proyecto```

2. Instala las dependencias:

   ```bash
   pip install -r requirements.txt```

3. Configura los archivos `.json` en la carpeta `configs/`. Si los archivos `.json` no existen, el programa generará automáticamente estos archivos a partir de los archivos `.json.sample`.

## Configuración

Cada archivo de configuración `.json` contiene parámetros necesarios para diferentes pasos del pipeline. Asegúrate de personalizar los valores en función de tus necesidades antes de ejecutar el programa.

### Archivos de Configuración

- **general.json**: Contiene configuraciones globales, como la ruta a los proyectos y archivos de referencia.
- **trimmomatic.json**, **snippy.json**, etc.: Configuraciones específicas para cada herramienta usada en el pipeline.

Si no existen los archivos `.json`, puedes crearlos automáticamente desde los `.json.sample` proporcionados. El programa emitirá una advertencia y los generará.

## Uso

### Ejecución del Pipeline Completo

El pipeline se puede ejecutar desde la línea de comandos con los siguientes parámetros:

```bash
python execute_pipeline.py PROJECT_NAME operation [--reference REFERENCE_PATH] [--log-level LOG_LEVEL] [--force]
```

- **PROJECT_NAME**: Nombre del proyecto. Se usará como nombre de carpeta en el `PROJECTS_PATH` especificado en `general.json`.
- **operation**: Operación a ejecutar. Puede ser una sola operación o varias, separadas por comas.
- **--reference**: Ruta al archivo de referencia (opcional).
- **--log-level**: Nivel de logging (`DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`).
- **--force**: Fuerza la ejecución de procesos.

### Ejemplo

Para ejecutar la operación `trimmomatic` en el proyecto `Muestra_A` con nivel de log en `DEBUG`:

```bash
python execute_pipeline.py Muestra_A trimmomatic --log-level DEBUG
```

### Operaciones Disponibles

Las operaciones incluyen:

- `create_project`: Crea la estructura inicial de un proyecto.
- `create_sample_list`: Genera una lista de muestras.
- `generate_excel`: Genera un reporte consolidado en Excel.
- `trimmomatic`: Realiza trimming de las lecturas.
- `SPAdes`: Realiza ensamblado de novo.
- `bowtie`: Realiza alineación de lecturas.
- `resfinder`, `oprD`, `mlst`, `snippy`, `PDC`: Otros análisis de resistencia.

Para más información sobre cada operación, consulta la documentación en cada script en `programs_scripts/`.

## Logs

Los archivos de log se almacenan en `logs/`. Se generan de forma automática en cada ejecución y se almacenan en la subcarpeta `Logs/` dentro de cada proyecto.

## Contacto

Para preguntas o problemas, abre un *issue* en el repositorio o contacta con el equipo de desarrollo en [matiasbonet@oceandrivers.com](matiasbonet@oceandrivers.com).

