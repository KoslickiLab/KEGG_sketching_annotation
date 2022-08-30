#!/usr/bin/env bash
set -e
set -u
dataDir="../data"
mkdir -p $dataDir
# File extensions must be fa/fasta/faa otherwise bbmap doesn't know how to parse them
genomeDatabaseFull="$dataDir/genome_ref_full.fa"
genomeDatabaseTruncated="$dataDir/genome_ref_truncated.fa"
proteinDatabaseFull="$dataDir/protein_ref_full.faa"
proteinDatabaseTruncated="$dataDir/protein_ref_truncated.faa"
simulatedMetagenome="$dataDir/simulatedMetagenome.fastq"
genomePath="$dataDir/reference_genomes"

# set variables
numGenomes=10
numGenomesFullDB=$numGenomes
numGenomesTruncatedDB=8
numReads=1000
readLen=150
numGenes=10
kSize=7  # decreasing this increases sensitivity at cost of FP's
refScale=10  # the higher this number, the faster things run, the smaller the database, at the cost of less sensitivity
queryScale=1  # likely will want to keep this at one (no down-sampling of the query)
thresholdBP=100  # this has the largest impact on FNs and FPs: setting it higher filters out more false positives, at the cost of more false negatives
# download genomes
echo "./get_reference_genomes.py -n $numGenomes -s $dataDir  -u"
#./get_reference_genomes.py -n $numGenomes -s $dataDir  -u

# create the full genome reference database
echo "./create_genome_ref_db.py $genomePath $genomeDatabaseFull $numGenomesFullDB"
#./create_genome_ref_db.py $dataDir $genomeDatabaseFull $numGenomesFullDB

# create the truncated genome reference database
echo "./create_genome_ref_db.py $genomePath $genomeDatabaseTruncated $numGenomesTruncatedDB"
#./create_genome_ref_db.py $dataDir $genomeDatabaseTruncated $numGenomesTruncatedDB

# create the mapping files required for the protein database
echo "./make_mapping_file.py "$genomePath""
#./make_mapping_file.py "$dataDir/reference_genomes

# create the protein reference database
echo "./create_gene_ref_db.py "$dataDir/reference_genomes" $proteinDatabaseFull $numGenomesFullDB"
#./create_gene_ref_db.py "$dataDir/reference_genomes" $proteinDatabaseFull $numGenomesFullDB
echo "./create_gene_ref_db.py "$dataDir/reference_genomes" $proteinDatabaseTruncated $numGenomesTruncatedDB"
#$dataDir/reference_genomes" $proteinDatabaseTruncated $numGenomesTruncatedDB

# simulate a metagenome
echo "./simulate_metagenome.py -r $genomeDatabaseFull -o $simulatedMetagenome -n $numReads -l $readLen --num_orgs $numGenes"
#./simulate_metagenome.py -r $genomeDatabaseFull -o $simulatedMetagenome -n $numReads -l $readLen --num_orgs $numGenes

# get the abundance estimates for the simulated metagenome
echo "./find_genes_in_sim.py --database_dir $genomePath --simulation $simulatedMetagenome --output_file $dataDir/ground_truth.csv"
#./find_genes_in_sim.py --database_dir $genomePath --simulation $simulatedMetagenome --output_file $dataDir/ground_truth.csv

# Run sourmash
echo "/usr/bin/time ./classify_sourmash.py -r $proteinDatabaseTruncated -m $simulatedMetagenome -o $dataDir -k $kSize --ref_scale_size $refScale --query_scale_size $queryScale --query_translate -t $thresholdBP"
#/usr/bin/time ./classify_sourmash.py -r $proteinDatabaseTruncated -m $simulatedMetagenome -o $dataDir -k $kSize --ref_scale_size $refScale --query_scale_size $queryScale --query_translate -t $thresholdBP

# Calculate sourmash performance metrics
gatherFile="$dataDir/$(basename $simulatedMetagenome)_k_${kSize}_scale_${queryScale}.sig_$(basename $proteinDatabase)_k_${kSize}_scale_${refScale}.sig_gather.csv"
#echo "gatherFile: $gatherFile"
echo "./calculate_sourmash_performance.py -g $dataDir/ground_truth.csv -s $gatherFile -o $dataDir/sourmash_performance_metrics.csv"
#./calculate_sourmash_performance.py -g $dataDir/ground_truth.csv -s $gatherFile -o $dataDir/sourmash_performance_metrics.csv

# Run Diamond
/usr/bin/time ./classify_diamond.py -r $proteinDatabaseTruncated -m $simulatedMetagenome -o $dataDir
diamondFile="$dataDir/$(basename $simulatedMetagenome)_$(basename $proteinDatabase).dmnd_matches.csv"
echo "diamondFile: $diamondFile"

# Calculate diamond performance metrics
./calculate_diamond_performance.py -g $dataDir/ground_truth.csv -s $diamondFile -o $dataDir/diamond_performance_metrics.csv
