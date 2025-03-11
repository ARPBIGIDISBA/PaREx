<<<<<<< HEAD
#!/bin/bash

# Directorio de instalación por defecto
INSTALL_DIR=/root/code/Programs
=======
#!/bin/bash -x
VERBOSE=0
OS_NAME=
OS_VERSION=
OS_CODE_NAME=
ARCH=
export IP=
export DEBIAN_FRONTEND=noninteractive
[[ -z "$CREATE_AMI" ]] && export CREATE_AMI
[[ -z "$DB_ENGINE" ]] && export DB_ENGINE="MYSQL_8.0"
>>>>>>> add63062d4fd9cc195a972d4b8c0c5611714c98c

# Si se proporciona un argumento, usarlo como directorio de instalación
if [ ! -z "$1" ]; then
    INSTALL_DIR=$1
fi

echo "📂 Creando directorio de instalación en: $INSTALL_DIR"
mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

# Preguntar qué programas instalar
echo "Seleccione los programas que desea instalar (s/n):"
read -p "¿Instalar FastQC? (s/n): " INSTALL_FASTQC
read -p "¿Instalar Trimmomatic? (s/n): " INSTALL_TRIMMOMATIC
read -p "¿Instalar SPAdes? (s/n): " INSTALL_SPADES
read -p "¿Instalar Resfinder? (s/n): " INSTALL_RESFINDER
read -p "¿Instalar SeqMonk? (s/n): " INSTALL_SEQMONK
read -p "¿Instalar Snippy? (s/n): " INSTALL_SNIPPY
read -p "¿Instalar MLST? (s/n): " INSTALL_MLST

# Actualizar el sistema e instalar dependencias
echo "🔄 Actualizando el sistema y instalando dependencias..."
sudo apt update -y && sudo apt install -y default-jre default-jdk unzip build-essential cmake \
    zlib1g-dev libbz2-dev python3 python3-pip libgsl-dev libsqlite3-dev libxml2-dev libxslt1-dev git && \
    echo "✅ Dependencias instaladas correctamente." || echo "❌ Error instalando dependencias."

if [[ "$INSTALL_FASTQC" == "s" ]]; then
    # Instalar FastQC
    echo "📥 Instalando FastQC..."
    wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip -q fastqc_v0.12.1.zip && \
    chmod +x FastQC/fastqc && \
    sudo ln -sf "$INSTALL_DIR/FastQC/fastqc" /usr/local/bin/fastqc && \
    echo "✅ FastQC instalado." || echo "❌ Error instalando FastQC."
fi

if [[ "$INSTALL_TRIMMOMATIC" == "s" ]]; then
    # Instalar Trimmomatic
    echo "📥 Instalando Trimmomatic..."
    wget -q http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip -q Trimmomatic-0.39.zip && \
    sudo ln -sf "$INSTALL_DIR/Trimmomatic-0.39/trimmomatic-0.39.jar" /usr/local/bin/trimmomatic && \
    echo "✅ Trimmomatic instalado." || echo "❌ Error instalando Trimmomatic."
fi

if [[ "$INSTALL_SPADES" == "s" ]]; then
    # Instalar SPAdes
    echo "📥 Instalando SPAdes..."
    wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz && \
    tar -xzf SPAdes-4.0.0-Linux.tar.gz && \
    mv SPAdes-4.0.0-Linux SPAdes-4.0.0 && \
    sudo ln -sf "$INSTALL_DIR/SPAdes-4.0.0/bin/spades.py" /usr/local/bin/spades.py && \
    cd .. && echo "✅ SPAdes instalado." || echo "❌ Error instalando SPAdes."
fi

if [[ "$INSTALL_RESFINDER" == "s" ]]; then

   # Instalar Resfinder
    echo "📥 Instalando Resfinder..."

    # Instalar dependencias necesarias
    pip3 install --upgrade pandas cgelib resfinder tabulate biopython cgecore gitpython python-dateutil && \

    # Clonar el repositorio principal de Resfinder
    git clone https://bitbucket.org/genomicepidemiology/resfinder.git "$INSTALL_DIR/resfinder" && \

    # Moverse al directorio de Resfinder
    cd "$INSTALL_DIR/resfinder" && \

    # Clonar bases de datos y KMA dentro del directorio correcto
    git clone https://bitbucket.org/genomicepidemiology/kma.git kma && \
    cd kma && make && cd .. && \
    git clone https://bitbucket.org/genomicepidemiology/resfinder_db.git db_resfinder && \
    git clone https://bitbucket.org/genomicepidemiology/pointfinder_db.git db_pointfinder && \

    # Actualizar base de datos Resfinder
    cd db_resfinder && python3 update_db.py && cd "$INSTALL_DIR/resfinder"

    # Configurar variables de entorno en ~/.bashrc (solo si no existen)
    echo "📌 Configurando variables de entorno..."
    if ! grep -q "export CGE_BLASTN=" ~/.bashrc; then
        echo "export CGE_BLASTN=\"/usr/bin/blastn\"" >> ~/.bashrc
    fi
    if ! grep -q "export CGE_KMA=" ~/.bashrc; then
        echo "export CGE_KMA=\"$INSTALL_DIR/resfinder/kma/\"" >> ~/.bashrc
    fi
    if ! grep -q "export CGE_RESFINDER_RESGENE_DB=" ~/.bashrc; then
        echo "export CGE_RESFINDER_RESGENE_DB=\"$INSTALL_DIR/resfinder/db_resfinder/\"" >> ~/.bashrc
    fi
    if ! grep -q "export CGE_RESFINDER_RESPOINT_DB=" ~/.bashrc; then
        echo "export CGE_RESFINDER_RESPOINT_DB=\"$INSTALL_DIR/resfinder/db_pointfinder/\"" >> ~/.bashrc
    fi

    source ~/.bashrc
    
    # Aplicar cambios de variables de entorno en la sesión actual
    export CGE_KMA="$INSTALL_DIR/resfinder/kma/"
    export CGE_RESFINDER_RESGENE_DB="$INSTALL_DIR/resfinder/db_resfinder/"
    export CGE_RESFINDER_RESPOINT_DB="$INSTALL_DIR/resfinder/db_pointfinder/"

    echo "✅ Resfinder instalado correctamente." || echo "❌ Error instalando Resfinder."
fi

if [[ "$INSTALL_SEQMONK" == "s" ]]; then
  
    ##############
    # Instalar SeqMonk

    # Instalar R
    echo "📥 Instalando R..."
    sudo sudo apt install -y r-base && \
    echo "✅ R instalado." || echo "❌ Error instalando R."

    # Instalar dependencias de desarrollo para R
    echo "📥 Instalando dependencias de desarrollo para R..."
    sudo apt install -y build-essential libxml2-dev && \
    echo "✅ Dependencias de desarrollo instaladas." || echo "❌ Error instalando dependencias de desarrollo."

    # Instalar SeqMonk
    echo "📥 Instalando SeqMonk..."
    wget -q https://www.bioinformatics.babraham.ac.uk/projects/seqmonk/seqmonk_v1.48.1_linux64.tar.gz && \
    tar -xzf seqmonk_v1.48.1_linux64.tar.gz && \
    sudo ln -sf "$INSTALL_DIR/SeqMonk/seqmonk" /usr/local/bin/seqmonk && \
    echo "✅ SeqMonk instalado." || echo "❌ Error instalando SeqMonk."
fi

if [[ "$INSTALL_SNIPPY" == "s" ]]; then

    # Instalar Snippy
    echo "📥 Instalando Snippy..."

    # Instalar dependencias
    sudo sudo apt install -y samtools && \
    sudo cpan App::cpanminus && \
    sudo cpanm Bio::SeqIO && \
    echo "✅ Dependencias de Snippy instaladas." || echo "❌ Error instalando dependencias de Snippy."

    # Clonar repositorio de Snippy
    cd "$INSTALL_DIR"
    git clone https://github.com/tseemann/snippy.git && \
    echo "✅ Snippy clonado correctamente." || echo "❌ Error clonando Snippy."

    # Crear un enlace simbólico para ejecutar Snippy fácilmente
    sudo ln -sf "$INSTALL_DIR/snippy/bin/snippy" /usr/local/bin/snippy && \
    echo "✅ Enlace simbólico de Snippy creado." || echo "❌ Error creando enlace de Snippy."

    echo "🎉 Snippy instalado correctamente."
fi

if [[ "$INSTALL_MLST" == "s" ]]; then
    # Instalar MLST
    echo "📥 Instalando MLST..."

    # Instalar dependencias necesarias
    sudo sudo apt install -y ncbi-blast+ && \
    sudo cpan install Moo && \
    sudo cpan install JSON && \
    echo "✅ Dependencias de MLST instaladas." || echo "❌ Error instalando dependencias de MLST."

    # Clonar repositorio de MLST y Any2Fasta
    cd "$INSTALL_DIR"
    git clone https://github.com/tseemann/mlst.git && \
    git clone https://github.com/tseemann/any2fasta.git && \
    echo "✅ MLST y Any2Fasta clonados correctamente." || echo "❌ Error clonando MLST o Any2Fasta."

    # Crear enlaces simbólicos para ejecución rápida
    sudo ln -sf "$INSTALL_DIR/mlst/bin/mlst" /usr/local/bin/mlst && \
    sudo ln -sf "$INSTALL_DIR/any2fasta/any2fasta" /usr/local/bin/any2fasta && \
    echo "✅ Enlaces simbólicos de MLST y Any2Fasta creados." || echo "❌ Error creando enlaces simbólicos."

    echo "🎉 MLST instalado correctamente."
fi
