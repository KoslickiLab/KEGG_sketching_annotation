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
    return [ get_diamond_running_time(num_reads, seed) for seed in seeds_list ]

def get_sourmash_running_time(num_reads, ksize, seed):
    filename = data_dir + f"/sourmash_gather_benchmark_{num_reads}_seed_{seed}_k_{ksize}"
    f = open(filename, 'r')
    time_str = f.readlines()[1].split('\t')[0]
    f.close()
    return float(time_str)

def get_all_sourmash_running_times(num_reads, ksize):
    return [get_sourmash_running_time(num_reads, ksize, seed) for seed in seeds_list]

def get_diamond_precision(num_reads, seed):
    filename = data_dir+f'/diamond_performance_metrics_num_reads_{num_reads}_seed_{seed}'
    f = open(filename, 'r')
    prec_str = f.readlines()[2].split(',')[3]
    f.close()
    return float(prec_str)

def get_diamond_recall(num_reads, seed):
    filename = data_dir+f'/diamond_performance_metrics_num_reads_{num_reads}_seed_{seed}'
    f = open(filename, 'r')
    prec_str = f.readlines()[2].split(',')[4]
    f.close()
    return float(prec_str)

def get_diamond_F1(num_reads, seed):
    filename = data_dir+f'/diamond_performance_metrics_num_reads_{num_reads}_seed_{seed}'
    f = open(filename, 'r')
    prec_str = f.readlines()[2].split(',')[5]
    f.close()
    return float(prec_str)


if __name__ == "__main__":
    num_reads = 100000
    for seed in seeds_list:
        print(get_diamond_recall(num_reads, seed))
