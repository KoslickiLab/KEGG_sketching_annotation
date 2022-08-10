#!/usr/bin/env bash
set -e
set -u
dataDir="../data"
mkdir -p $dataDir
# File extensions must be fa/fasta/faa otherwise bbmap doesn't know how to parse them
genomeDatabase="$dataDir/genome_ref.fa"
proteinDatabase="$dataDir/protein_ref.faa"
simulatedMetagenome="$dataDir/simulatedMetagenome.fastq"

# set variables
numGenomes=10
numReads=10000
readLen=150
numGenes=100
kSize=11
refScale=1
queryScale=1
# download genomes
#./get_reference_genomes.py -n $numGenomes -s $dataDir  -u

# create the genome reference database
#./create_genome_ref_db.py $dataDir $genomeDatabase

# create the mapping files required for the protein database
#./make_mapping_file.py "$dataDir/reference_genomes"

# create the protein reference database
#./create_gene_ref_db.py "$dataDir/reference_genomes" $proteinDatabase

# simulate a metagenome
#./simulate_metagenome.py -r $genomeDatabase -o $simulatedMetagenome -n $numReads -l $readLen --num_orgs $numGenes

# get the abundance estimates for the simulated metagenome
#./find_genes_in_sim.py --database_dir "$dataDir/reference_genomes" --simulation $simulatedMetagenome --output_file "$dataDir/ground_truth.csv"

# Run sourmash
/usr/bin/time ./classify_sourmash.py -r $proteinDatabase -m $simulatedMetagenome -o $dataDir -k $kSize --ref_scale_size $refScale --query_scale_size $queryScale --query_translate

