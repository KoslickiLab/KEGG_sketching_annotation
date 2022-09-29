#!/usr/bin/env bash
set -e
set -u
scriptDir="../../scripts"
dataDir="data"
mkdir -p $dataDir

# File extensions must be fa/fasta/faa otherwise bbmap doesn't know how to parse them
genomeDatabaseFull="$dataDir/genome_ref_full.fa"
genomeDatabaseTruncated="$dataDir/genome_ref_truncated.fa"
proteinDatabaseFull="$dataDir/protein_ref_full.faa"
proteinDatabaseTruncated="$dataDir/protein_ref_truncated.faa"
simulatedMetagenome="$dataDir/simulatedMetagenome.fastq"
genomePath="../reference_genomes"

# set variables
numGenomes=5
numGenomesFullDB=$numGenomes
numGenomesTruncatedDB=3
numReads=1000000
readLen=150
numGenes=10000
kSize=7  # decreasing this increases sensitivity at cost of FP's
refScale=10  # the higher this number, the faster things run, the smaller the database, at the cost of less sensitivity
queryScale=1  # likely will want to keep this at one (no down-sampling of the query)
thresholdBP=50  # this has the largest impact on FNs and FPs: setting it higher filters out more false positives, at the cost of more false negatives

# create the full genome reference database
echo "$scriptDir/create_genome_ref_db.py $genomePath $genomeDatabaseFull $numGenomesFullDB"
$scriptDir/create_genome_ref_db.py $genomePath $genomeDatabaseFull $numGenomesFullDB

# create the truncated genome reference database
echo "$scriptDir/create_genome_ref_db.py $genomePath $genomeDatabaseTruncated $numGenomesTruncatedDB"
$scriptDir/create_genome_ref_db.py $genomePath $genomeDatabaseTruncated $numGenomesTruncatedDB

# create the mapping files required for the protein database
echo "$scriptDir/make_mapping_file.py "$genomePath""
#$scriptDir/make_mapping_file.py "$genomePath"
# no need because all genomes have mapping files ready

# create the protein reference database
echo "$scriptDir/create_gene_ref_db.py "$genomePath" $proteinDatabaseFull $numGenomesFullDB"
$scriptDir/create_gene_ref_db.py "$genomePath" $proteinDatabaseFull $numGenomesFullDB
echo "$scriptDir/create_gene_ref_db.py "$genomePath" $proteinDatabaseTruncated $numGenomesTruncatedDB"
$scriptDir/create_gene_ref_db.py "$genomePath" $proteinDatabaseTruncated $numGenomesTruncatedDB

for numGenes in 500 1000; do
    # simulate a metagenome
    echo "$scriptDir/simulate_metagenome.py -r $genomeDatabaseFull -o $simulatedMetagenome -n $numReads -l $readLen --num_orgs $numGenes"
    $scriptDir/simulate_metagenome.py -r $genomeDatabaseFull -o $simulatedMetagenome -n $numReads -l $readLen --num_orgs $numGenes

    # get the abundance estimates for the simulated metagenome
    echo "$scriptDir/find_genes_in_sim.py --database_dir $genomePath --simulation $simulatedMetagenome --output_file $dataDir/ground_truth.csv --num_genomes $numGenomesFullDB"
    $scriptDir/find_genes_in_sim.py --database_dir $genomePath --simulation $simulatedMetagenome --output_file $dataDir/ground_truth.csv --num_genomes $numGenomesFullDB

    # Run sourmash
    echo "/usr/bin/time $scriptDir/classify_sourmash.py -r $proteinDatabaseTruncated -m $simulatedMetagenome -o $dataDir -k $kSize --ref_scale_size $refScale --query_scale_size $queryScale --query_translate -t $thresholdBP"
    /usr/bin/time $scriptDir/classify_sourmash.py -r $proteinDatabaseTruncated -m $simulatedMetagenome -o $dataDir -k $kSize --ref_scale_size $refScale --query_scale_size $queryScale --query_translate -t $thresholdBP

    # Calculate sourmash performance metrics
    gatherFile="$dataDir/$(basename $simulatedMetagenome)_k_${kSize}_scale_${queryScale}.sig_$(basename $proteinDatabaseTruncated)_k_${kSize}_scale_${refScale}.sig_gather.csv"
    echo "gatherFile: $gatherFile"
    echo "$scriptDir/calculate_sourmash_performance.py -g $dataDir/ground_truth.csv -s $gatherFile -o $dataDir/sourmash_performance_metrics.csv"
    $scriptDir/calculate_sourmash_performance.py -g $dataDir/ground_truth.csv -s $gatherFile -o $dataDir/sourmash_performance_metrics_numGenes_"$numGenes".csv

    # Run Diamond
    /usr/bin/time $scriptDir/classify_diamond.py -r $proteinDatabaseTruncated -m $simulatedMetagenome -o $dataDir
    diamondFile="$dataDir/$(basename $simulatedMetagenome)_$(basename $proteinDatabaseTruncated).dmnd_matches.csv"
    echo "diamondFile: $diamondFile"

    # Calculate diamond performance metrics
    $scriptDir/calculate_diamond_performance.py -g $dataDir/ground_truth.csv -s $diamondFile -o $dataDir/diamond_performance_metrics_numGenes_"$numGenes".csv
done
