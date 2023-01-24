script_dir = '../../scripts'
data_dir = 'data'
read_length = 150
k = 7
ref_scale = 10
query_scale = 1
threshold_bp = 50
num_organisms = 30
num_genomes_full_db = 30
num_genomes_truncated_db = 20
genome_path = '../extracted_genomes_from_kegg'
num_reads_list = ['50000', '100000', '150000']
seeds_list = [str(i) for i in range(10)]
kmer_sizes = ['5', '7', '9']

data_dir+'/diamond_performance_metrics_num_reads_{nr}_seed_{seed}'
data_dir+'/sourmash_performance_metrics_num_reads_{nr}_seed_{seed}_k_{k}'
data_dir+'/diamond_benchmark_num_reads_{num_reads}_seed_{seed}'
data_dir + "/sourmash_gather_benchmark_{num_reads}_seed_{seed}_k_{k}"

def get_diamond_running_time(num_reads, seed):
    filename = data_dir+f'/diamond_benchmark_num_reads_{num_reads}_seed_{seed}'
    f = open(filename, 'r')
    time_str = f.readlines()[1].split('\t')[0]
    f.close()
    return float(time_str)

def get_all_diamond_running_times(num_reads):
    for seed in seeds_list:
        print( get_diamond_running_time(num_reads, seed) )

if __name__ == "__main__":
    for num_reads in num_reads_list:
        print(f"Num reads: {num_reads}")
        get_all_diamond_running_times(num_reads)
    
