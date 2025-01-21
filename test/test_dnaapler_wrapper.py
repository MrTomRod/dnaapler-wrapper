import gzip
import logging
import urllib.request
from dnaapler_wrapper import *
from unittest import TestCase

logging.basicConfig(level=logging.INFO)

test_fasta = 'data/input.fasta'
workdir = 'data/workdir'
headers = [
    '>contig_1 [topology=circular]',
    '>contig_2 [topology=linear]',
    '>contig_3 [topology=circular]',
    '>contig_4 [topology=linear]',
    '>contig_5',
    '>contig_6',
    '>contig_7'
]
header_ids = [h.strip('>').split(' ')[0] for h in headers]

test_output_fasta = f'{workdir}/out/dnaapler_reoriented.fasta'
test_output_table = f'{workdir}/out/dnaapler_all_reorientation_summary.tsv'
expected_columns = ['Gene_Reoriented', 'Start', 'Strand', 'Top_Hit', 'Top_Hit_Length', 'Covered_Length',
                    'Coverage', 'Identical_AAs', 'Identity_Percentage', 'Overlapping_Contig_End']


def get_test_fasta():
    # download and create new headers:
    url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/020/362/645/GCA_020362645.3_ASM2036264v3/GCA_020362645.3_ASM2036264v3_genomic.fna.gz'
    os.makedirs('data', exist_ok=True)
    urllib.request.urlretrieve(url, test_fasta + '.gz')
    with (gzip.open(test_fasta + '.gz', 'rb') as fi, open(test_fasta, 'w') as fo):
        contig_counter = 0
        for line in fi:
            line = line.decode('utf-8')
            if line.startswith('>'):
                fo.write(f'{headers[contig_counter]}\n')
                contig_counter += 1
            else:
                fo.write(line)
    os.remove(test_fasta + '.gz')


if not os.path.exists(test_fasta):
    get_test_fasta()

if not all(os.path.exists(f) for f in [test_output_fasta, test_output_table]):
    subprocess.run(['dnaapler', 'all', '-i', 'input.fasta', '-o', 'out', '-t', '8', '-f'], cwd=workdir)


class Test(TestCase):
    def test_read_input_fasta(self):
        fasta = read_fasta(test_fasta)
        self.assertEqual(header_ids, [c for c, d in fasta.items()])
        self.assertEqual(headers, ['>' + d['header'] for c, d in fasta.items()])

    def test_read_output_table(self):
        output_table = read_output_table(test_output_table)
        print(output_table.to_string())
        self.assertEqual(output_table.columns.tolist(), expected_columns)
        self.assertEqual([f'>{h}' for h in output_table.index], headers)

    def test_read_output_fasta(self):
        fasta = read_fasta(test_output_fasta)
        for hi, ho in zip(headers, fasta.keys()):
            self.assertTrue(hi.startswith(f'>{ho}'), msg=f'"{ho}" does not start with "{hi}"')

    def test_merge_fastas(self):
        out_fasta = read_fasta(test_output_fasta)
        out_table = read_output_table(test_output_table)
        merged_fasta = merge_fastas(read_fasta(test_fasta), out_fasta, out_table, default_topology='circular')
        self.assertEqual(header_ids, list(merged_fasta))
        expected_result_both = [
            'contig_1 [topology=circular] [dnaapler=already-reoriented]',
            'contig_2 [topology=linear] [dnaapler=reverse-complement] [dnaapler-gene=repA]',
            'contig_3 [topology=circular] [dnaapler=rotated] [dnaapler-gene=repA]',
            'contig_4 [topology=linear] [dnaapler=reverse-complement] [dnaapler-gene=repA]',
            'contig_5 [dnaapler=failed-no-hits]',
        ]
        expected_result = expected_result_both + [
            'contig_6 [dnaapler=rotated] [dnaapler-gene=repA]',
            'contig_7 [dnaapler=rotated] [dnaapler-gene=repA]'
        ]
        for i, (c, d) in enumerate(merged_fasta.items()):
            self.assertEqual(expected_result[i], d['header'])
        merged_fasta = merge_fastas(read_fasta(test_fasta), out_fasta, out_table, default_topology='linear')
        self.assertEqual(header_ids, list(merged_fasta))
        expected_result = expected_result_both + [
            'contig_6 [dnaapler=reverse-complement] [dnaapler-gene=repA]',  # different
            'contig_7 [dnaapler=reverse-complement] [dnaapler-gene=repA]'  # different
        ]
        for i, (c, d) in enumerate(merged_fasta.items()):
            self.assertEqual(expected_result[i], d['header'])

    def test_merge_fastas_fail(self):
        in_fasta = read_fasta(test_fasta)
        out_fasta = read_fasta(test_output_fasta)
        out_table = read_output_table(test_output_table)
        with self.assertRaises(ValueError):
            merge_fastas(in_fasta, out_fasta, out_table, default_topology='raise')

    def test_write_fasta(self):
        out_fasta = read_fasta(test_output_fasta)
        out_table = read_output_table(test_output_table)
        merged_fasta = merge_fastas(read_fasta(test_fasta), out_fasta, out_table, default_topology='circular')
        write_fasta(merged_fasta, 'data/output.fasta')

    # def test_run_dnaapler(self):  # SLOW
    #     fasta, output_table = run_dnaapler(test_fasta)
    #     self.assertEqual(output_table.columns.tolist(), expected_columns)
    #     self.assertEqual([f'>{h}' for h in output_table.index], headers)
    #     for hi, ho in zip(headers, fasta.keys()):
    #         self.assertTrue(ho.startswith(hi), msg=f'"{ho}" does not start with "{hi}"')
