#!/usr/bin/env bash
set -e
set -u
scriptDir=$1
dataDir="data"
mkdir -p $dataDir

# File extensions must be fa/fasta/faa otherwise bbmap doesn't know how to parse them
genomeDatabaseFull="$dataDir/genome_ref_full.fa"
genomeDatabaseTruncated="$dataDir/genome_ref_truncated.fa"
proteinDatabaseFull="$dataDir/protein_ref_full.faa"
proteinDatabaseTruncated="$dataDir/protein_ref_truncated.faa"
simulatedMetagenome="$dataDir/simulatedMetagenome.fastq"
genomePath="../extracted_genomes_from_kegg"

# set variables
numGenomes=30
numGenomesFullDB=$numGenomes
numGenomesTruncatedDB=20
#numReads=100000
readLen=150
numGenes=1000
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
