#!/bin/sh

# Download all barrel information from IsItABarrelDB
wget -O iiab.tsv https://isitabarrel.ku.edu/download

# Extract AlphaFoldDB IDs from each IsItABarrel entry
tail +2 iiab.tsv | cut -f 16 | sort -u > afdb_ids.txt

# Create folder to store PDB structures
mkdir structures


# Download AFDB structures (remove to run local version below)
while read afdb_id; do
    wget -O structures/$afdb_id.pdb https://alphafold.ebi.ac.uk/files/AF-$afdb_id-F1-model_v4.pdb
done < afdb_ids.txt

# Alternatively, if you already have a local version of AlphaFoldDB
# you can copy the structures needed from AFDB into the subfolder.

# Copy local AFDB structures (uncomment and add AFDB source folder to run)
# while read afdb_id; do
#    cp <local AFDB folder>/$afdb_id.pdb structures/
# done < afdb_ids.txt