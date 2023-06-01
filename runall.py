import os

def execute(cmd):
    print(cmd)
    os.system(cmd)

source_dir = 'clean_reads_assembly/'
genes_dir = source_dir[:-1] + '_genes/'
log_dir = source_dir[:-1] + '_log/'
os.makedirs(genes_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)

files = []
do_gzip = False
for file in os.listdir(source_dir):
    if '.gz' in file: execute("gunzip " + source_dir + '*')
    files.append(file)
files.sort()
print("Found", files)


for file in files:
    if os.path.exists(log_dir + file.split('.')[0] + '_Results.txt'): continue
    print('\nRunning MetaCompare for', source_dir + file)
    file_genes = file[:-3] + '_genes.fa'
    file_log_prodigal = file[:-3] + '_plog.txt'
    file_log_metacompare = file[:-3] + '_mlog.txt'
    execute("prodigal -i %s -d %s -p meta > %s" % (source_dir + file, genes_dir + file_genes, log_dir + file_log_prodigal))
    execute("python metacmp.py -t 24 -c %s -g %s -l %s > %s" % 
        (source_dir + file, genes_dir + file_genes, log_dir, log_dir + file_log_metacompare))