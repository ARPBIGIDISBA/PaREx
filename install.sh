#!/bin/bash

# Directorio de instalación por defecto
# INSTALL_DIR=/root/code/Programs

# Si no se proporciona un argumento, salir y pedirlo
if [ -z "$1" ]; then
    echo "❌  You must specify an installation directory as an argument."
    echo "Usage: $0 /path/to/directory"
    exit 1
fi

# Si se proporciona un argumento, usarlo como directorio de instalación
if [ ! -z "$1" ]; then
    INSTALL_DIR=$1
fi
# translate to english
echo "📂 Creating installation directory at: $INSTALL_DIR"

mkdir -p "$INSTALL_DIR"
cd "$INSTALL_DIR"

# Ask which programs to install
echo "Select the programs you want to install (y/n):"
read -p "Install FastQC? (y/n): " INSTALL_FASTQC
read -p "Install Trimmomatic? (y/n): " INSTALL_TRIMMOMATIC
read -p "Install SPAdes? (y/n): " INSTALL_SPADES
read -p "Install Resfinder? (y/n): " INSTALL_RESFINDER
# read -p "Install SeqMonk? (y/n): " INSTALL_SEQMONK
read -p "Install Snippy? (y/n): " INSTALL_SNIPPY
read -p "Install MLST? (y/n): " INSTALL_MLST

# Update the system and install dependencies
echo "🔄 Updating the system and installing dependencies..."
sudo apt update -y && sudo apt install -y default-jre default-jdk unzip build-essential cmake \
    zlib1g-dev libbz2-dev python3 python3-pip libgsl-dev libsqlite3-dev libxml2-dev libxslt1-dev git && \
    echo "✅ Dependencies installed successfully." || echo "❌ Error installing dependencies."

if [[ "$INSTALL_FASTQC" == "y" ]]; then
    # Install FastQC
    echo "📥 Installing FastQC..."
    wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip -q fastqc_v0.12.1.zip && \
    chmod +x FastQC/fastqc && \
    sudo ln -sf "$INSTALL_DIR/FastQC/fastqc" /usr/local/bin/fastqc && \
    echo "✅ FastQC installed." || echo "❌ Error installing FastQC."
fi

if [[ "$INSTALL_TRIMMOMATIC" == "y" ]]; then
    # Install Trimmomatic
    echo "📥 Installing Trimmomatic..."
    wget -q http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip -q Trimmomatic-0.39.zip && \
    sudo ln -sf "$INSTALL_DIR/Trimmomatic-0.39/trimmomatic-0.39.jar" /usr/local/bin/trimmomatic && \
    echo "✅ Trimmomatic installed." || echo "❌ Error installing Trimmomatic."
fi

if [[ "$INSTALL_SPADES" == "y" ]]; then
    # Install SPAdes 3.15.5
    echo "📥 Installing SPAdes..."
    wget https://github.com/ablab/spades/releases/download/v3.15.5/SPAdes-3.15.5-Linux.tar.gz && \
    tar -xzf SPAdes-3.15.5-Linux.tar.gz && \
    mv SPAdes-3.15.5-Linux SPAdes-3.15.5 && \
    sudo ln -sf "$INSTALL_DIR/SPAdes-3.15.5/bin/spades.py" /usr/local/bin/spades.py && \
    cd .. && echo "✅ SPAdes installed." || echo "❌ Error installing SPAdes."
fi

if [[ "$INSTALL_RESFINDER" == "y" ]]; then
    # Install Resfinder
    echo "📥 Installing Resfinder..."

    # Install required dependencies
    pip3 install --upgrade pandas cgelib resfinder tabulate biopython cgecore gitpython python-dateutil && \

    # Clone the main Resfinder repository
    git clone https://bitbucket.org/genomicepidemiology/resfinder.git "$INSTALL_DIR/resfinder" && \

    # Move to the Resfinder directory
    cd "$INSTALL_DIR/resfinder" && \

    # Clone databases and KMA into the correct directory
    git clone https://bitbucket.org/genomicepidemiology/kma.git kma && \
    cd kma && make && cd .. && \
    git clone --branch 2.1.0 https://bitbucket.org/genomicepidemiology/resfinder_db.git db_resfinder && \
    git clone https://bitbucket.org/genomicepidemiology/pointfinder_db.git db_pointfinder && \

    # Update Resfinder database
    cd db_resfinder && python3 update_db.py && cd "$INSTALL_DIR/resfinder"

    # Set environment variables in ~/.bashrc (only if not present)
    echo "📌 Setting environment variables..."
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
    
    # Apply environment variable changes to current session
    export CGE_KMA="$INSTALL_DIR/resfinder/kma/"
    export CGE_RESFINDER_RESGENE_DB="$INSTALL_DIR/resfinder/db_resfinder/"
    export CGE_RESFINDER_RESPOINT_DB="$INSTALL_DIR/resfinder/db_pointfinder/"

    echo "✅ Resfinder installed successfully." || echo "❌ Error installing Resfinder."
fi

# if [[ "$INSTALL_SEQMONK" == "y" ]]; then
#     # Install SeqMonk

#     # Install R
#     echo "📥 Installing R..."
#     sudo apt install -y r-base && \
#     echo "✅ R installed." || echo "❌ Error installing R."

#     # Install development dependencies for R
#     echo "📥 Installing development dependencies for R..."
#     sudo apt install -y build-essential libxml2-dev && \
#     echo "✅ Development dependencies installed." || echo "❌ Error installing development dependencies."

#     # Install SeqMonk
#     echo "📥 Installing SeqMonk..."
#     wget -q https://www.bioinformatics.babraham.ac.uk/projects/seqmonk/seqmonk_v1.48.1_linux64.tar.gz && \
#     tar -xzf seqmonk_v1.48.1_linux64.tar.gz && \
#     sudo ln -sf "$INSTALL_DIR/SeqMonk/seqmonk" /usr/local/bin/seqmonk && \
#     echo "✅ SeqMonk installed." || echo "❌ Error installing SeqMonk."
# fi

# Install Snippy
if [[ "$INSTALL_SNIPPY" == "y" ]]; then
    # Install Snippy
    echo "📥 Installing Snippy..."

    # Install dependencies
    sudo apt install -y samtools && \
    sudo cpan App::cpanminus && \
    sudo cpanm Bio::SeqIO && \
    echo "✅ Snippy dependencies installed." || echo "❌ Error installing Snippy dependencies."

    # Clone Snippy repository
    cd "$INSTALL_DIR"
    git clone https://github.com/tseemann/snippy.git && \
    echo "✅ Snippy cloned successfully." || echo "❌ Error cloning Snippy."

    # Create a symlink for easy execution
    sudo ln -sf "$INSTALL_DIR/snippy/bin/snippy" /usr/local/bin/snippy && \
    echo "✅ Snippy symlink created." || echo "❌ Error creating Snippy symlink."

    echo "🎉 Snippy installed successfully."
fi

if [[ "$INSTALL_MLST" == "y" ]]; then
    # Install MLST
    echo "📥 Installing MLST..."

    # Install required dependencies
    sudo apt install -y ncbi-blast+ && \
    sudo cpan install Moo && \
    sudo cpan install JSON && \
    echo "✅ MLST dependencies installed." || echo "❌ Error installing MLST dependencies."

    # Clone MLST and Any2Fasta repositories
    cd "$INSTALL_DIR"
    git clone https://github.com/tseemann/mlst.git && \
    git clone https://github.com/tseemann/any2fasta.git && \
    echo "✅ MLST and Any2Fasta cloned successfully." || echo "❌ Error cloning MLST or Any2Fasta."

    # Create symlinks for quick execution
    sudo ln -sf "$INSTALL_DIR/mlst/bin/mlst" /usr/local/bin/mlst && \
    sudo ln -sf "$INSTALL_DIR/any2fasta/any2fasta" /usr/local/bin/any2fasta && \
    echo "✅ MLST and Any2Fasta symlinks created." || echo "❌ Error creating symlinks."

    echo "🎉 MLST installed successfully."
fi
