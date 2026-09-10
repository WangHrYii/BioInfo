"""Create small synthetic teaching inputs, without third-party packages."""
import argparse
import csv
import gzip
from pathlib import Path
import random


def write_table(path, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        csv.writer(handle, delimiter='\t').writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path('demo-data'))
    out = parser.parse_args().out
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(42)
    reference = ''.join(rng.choice('ACGT') for _ in range(3000))
    (out / 'reference.fa').write_text(
        '>chrDemo synthetic teaching chromosome\n' +
        '\n'.join(reference[i:i + 60] for i in range(0, len(reference), 60)) + '\n', encoding='utf-8')
    complement = str.maketrans('ACGT', 'TGCA')
    for sample in ['control', 'treated']:
        records = [[], []]
        for i in range(20):
            start = rng.randint(0, 2700)
            fragment = reference[start:start + 200]
            sequences = [fragment[:75], fragment[-75:].translate(complement)[::-1]]
            for mate, sequence in enumerate(sequences):
                quality = 'I' * 60 + '!' * 15 if sample == 'treated' and i == 0 else 'I' * 75
                records[mate].append('@{}_{}/{}\n{}\n+\n{}\n'.format(
                    sample, i + 1, mate + 1, sequence, quality))
        for mate, lines in enumerate(records, 1):
            with gzip.open(out / '{}_R{}.fastq.gz'.format(sample, mate), 'wt', encoding='ascii') as handle:
                handle.write(''.join(lines))
    write_table(out / 'genes.bed', [['chrDemo', 100, 400, 'geneA'], ['chrDemo', 800, 1200, 'geneB']])
    write_table(out / 'peaks.bed', [
        ['chrDemo', 150, 230, 'peak1'], ['chrDemo', 1100, 1250, 'peak2'], ['chrDemo', 2000, 2100, 'peak3']])
    write_table(out / 'genes.gtf', [
        ['chrDemo', 'demo', 'exon', 101, 400, '.', '+', '.', 'gene_id "geneA"; transcript_id "txA";'],
        ['chrDemo', 'demo', 'exon', 801, 1200, '.', '+', '.', 'gene_id "geneB"; transcript_id "txB";']])
    write_table(out / 'counts.tsv', [
        ['gene', 'C1', 'C2', 'C3', 'T1', 'T2', 'T3'],
        ['geneA', 100, 95, 110, 400, 380, 420], ['geneB', 200, 210, 195, 190, 205, 200],
        ['geneC', 80, 75, 90, 20, 25, 18]])
    write_table(out / 'samples.tsv', [
        ['sample', 'condition', 'r1', 'r2'],
        ['control', 'control', 'control_R1.fastq.gz', 'control_R2.fastq.gz'],
        ['treated', 'treated', 'treated_R1.fastq.gz', 'treated_R2.fastq.gz']])
    ref_base = reference[1199]
    alt_base = next(base for base in 'ACGT' if base != ref_base)
    (out / 'variants.vcf').write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=chrDemo,length=3000>\n'
        '##FILTER=<ID=PASS,Description="All filters passed">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">\n'
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tDEMO\n'
        'chrDemo\t1200\t.\t{}\t{}\t60\tPASS\t.\tGT:DP:AD\t0/1:20:11,9\n'.format(ref_base, alt_base),
        encoding='utf-8')
    (out / 'README.txt').write_text(
        'SYNTHETIC teaching data, seed=42. Not a biological experiment.\n'
        'FASTQ: 20 pairs/sample; 75 bases/read; 200-base fragments.\n'
        'Counts and VCF are independent examples, NOT inferred from these reads.\n'
        'The two FASTQ samples are for file/workflow exercises, not differential testing.\n', encoding='utf-8')
    print('Created teaching inputs in {}'.format(out.resolve()))


if __name__ == '__main__':
    main()
