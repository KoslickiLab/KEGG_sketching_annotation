#!/usr/bin/env bash
set -e
set -u
scriptDir="../../scripts"
dataDir="data"
mkdir -p $dataDir
# File extensions must be fa/fasta/faa otherwise bbmap doesn't know how to parse them
genomeDatabase="$dataDir/genome_ref.fa"
proteinDatabase="$dataDir/protein_ref.faa"
simulatedMetagenome="$dataDir/simulatedMetagenome.fastq"

# set variables
numGenomes=10
numReads=1000
readLen=150
numGenes=10
kSize=7
thresholdBP=100

# download genomes
#$scriptDir/./get_reference_genomes.py -n $numGenomes -s $dataDir  -u
# create the genome reference database
$scriptDir/./create_genome_ref_db.py $dataDir $genomeDatabase
# create the mapping files required for the protein database
$scriptDir/./make_mapping_file.py "$dataDir/reference_genomes"
# create the protein reference database
$scriptDir/./create_gene_ref_db.py "$dataDir/reference_genomes" $proteinDatabase
# simulate a metagenome
$scriptDir/./simulate_metagenome.py -r $genomeDatabase -o $simulatedMetagenome -n $numReads -l $readLen --num_orgs $numGenes
# get the abundance estimates for the simulated metagenome
$scriptDir/./find_genes_in_sim.py --database_dir "$dataDir/reference_genomes" --simulation $simulatedMetagenome --output_file "$dataDir/ground_truth.csv"

# Then iterate over the scale values and run sourmash
for refScale in 1 10 100 1000; do
    echo "Running with refScale = $refScale"
    for queryScale in 1 10 100 1000; do
        echo "Running with queryScale = $queryScale"
        /usr/bin/time $scriptDir/./classify_sourmash.py -r $proteinDatabase -m $simulatedMetagenome -o $dataDir -k $kSize --ref_scale_size $refScale --query_scale_size $queryScale --query_translate -t $thresholdBP
        gatherFile="$dataDir/$(basename $simulatedMetagenome)_k_${kSize}_scale_${queryScale}.sig_$(basename $proteinDatabase)_k_${kSize}_scale_${refScale}.sig_gather.csv"
#echo "gatherFile: $gatherFile"
        $scriptDir/./calculate_sourmash_performance.py -g $dataDir/ground_truth.csv -s $gatherFile -o $dataDir/sourmash_performance_metrics_refScale_${refScale}_queryScale_${queryScale}.csv
    done
done
