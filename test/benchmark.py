import glob
import random
import time
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from dnaapler_wrapper import run_dnaapler


def get_all_contigs(pattern):
    assemblies = glob.glob(pattern)
    all_contigs = []
    for file in assemblies:
        with open(file) as f:
            contigs = f.read().split('>')
            contigs = [c.split('\n', 1)[1] for c in contigs if c]
            all_contigs.extend(contigs)
    return all_contigs


def create_test_dataset(contigs, n, output_file):
    test_contigs = random.sample(contigs, n)
    with open(output_file, 'w') as f:
        for i, contig in enumerate(test_contigs):
            f.write(f">{i}\n{contig}\n")


all_contigs = get_all_contigs('data-pb-new/*/assembly-curator/hybrid.fasta')
times = {}


def measure_dnaapler_time(test_fasta, n_runs=5):
    times = []
    for i in range(n_runs):
        print(f"Running dnaapler {n=} run={i}")
        create_test_dataset(all_contigs, n, test_fasta)
        start = time.time()
        fasta_dnaapler, output_table = run_dnaapler(test_fasta)
        end = time.time()
        times.append(end - start)
        print(output_table)
    return times


for n in range(10, 101, 10):
    test_fasta = f'test/test-ds/ds_{n:02d}.fasta'
    times[n] = measure_dnaapler_time(test_fasta)

# Turn times into a pandas DataFrame
times_df = pd.DataFrame(times)

times_df.to_csv('test/test-ds/times.csv')

# Plot the times using seaborn
sns.boxplot(data=times_df)
plt.xlabel('Number of contigs')
plt.ylabel('Time (s)')
plt.title('Time taken by dnaapler')
plt.show()
