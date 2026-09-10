"""Validate paired FASTQ records and summarize reads and bases (Phred+33)."""
import argparse
import csv
import gzip
from itertools import zip_longest
from pathlib import Path


def records(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii') as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline().rstrip('\r\n')
            separator = handle.readline()
            quality_line = handle.readline()
            quality = quality_line.rstrip('\r\n')
            if (not header.startswith('@') or not separator.startswith('+') or
                    not quality_line or not sequence or len(sequence) != len(quality)):
                raise ValueError('Invalid or truncated four-line FASTQ record in {}'.format(path))
            name = header[1:].split()[0]
            if name.endswith(('/1', '/2')):
                name = name[:-2]
            scores = [ord(c) - 33 for c in quality]
            if any(q < 0 or q > 93 for q in scores):
                raise ValueError('Invalid Phred+33 quality in {}'.format(path))
            yield name, len(sequence), sum(q >= 30 for q in scores)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--r1', required=True, type=Path)
    parser.add_argument('--r2', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    pairs = bases = high_quality = 0
    for left, right in zip_longest(records(args.r1), records(args.r2)):
        if left is None or right is None or left[0] != right[0]:
            raise ValueError('Mate count/order/name mismatch for {}'.format(args.sample))
        pairs += 1
        bases += left[1] + right[1]
        high_quality += left[2] + right[2]
    if not pairs:
        raise ValueError('Empty paired FASTQ input')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(['sample', 'pairs', 'reads', 'bases', 'q30_fraction'])
        writer.writerow([args.sample, pairs, pairs * 2, bases, high_quality / bases])


if __name__ == '__main__':
    main()
