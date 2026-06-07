import urllib.request
import gzip
import random
import os

URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz"
REF_OUTPUT = "Inputs/Reference.txt"
READS_OUTPUT = "Inputs/Reads.txt"

READS_COUNT = 5000000
READ_LEN = 100

# Empirical Error Model Parameters
# Illumina sequencing error + human SNP rate combined approx 0.5% per base.
BASE_ERROR_RATE = 0.005 

def setup_dirs():
    os.makedirs("Inputs", exist_ok=True)
    os.makedirs("Outputs", exist_ok=True)

def download_and_extract_ref():
    print("Downloading chr21.fa.gz from UCSC...")
    tmp_file = "chr21.fa.gz"
    urllib.request.urlretrieve(URL, tmp_file)
    print("Extracting and processing sequence...")
    seq_lines = []
    with gzip.open(tmp_file, 'rt') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                seq_lines.append(line.upper())
    seq = "".join(seq_lines)
    seq = seq.replace('N', '')
    print(f"Total Reference Length (without N): {len(seq)} bp")
    
    print("Writing to Reference.txt...")
    with open(REF_OUTPUT, 'w') as f:
        f.write(">chr21\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")
            
    os.remove(tmp_file)
    return seq

def generate_reads_empirical(ref_seq):
    print(f"Generating {READS_COUNT} empirical reads (length={READ_LEN}, error_rate={BASE_ERROR_RATE})...")
    bases = ['A', 'C', 'G', 'T']
    
    with open(READS_OUTPUT, 'w') as f:
        for i in range(READS_COUNT):
            start = random.randint(0, len(ref_seq) - READ_LEN - 1)
            read = list(ref_seq[start:start+READ_LEN])
            
            # 독립 시행(Bernoulli Trial)을 통한 현실적 염기 단위 에러 모델링
            for j in range(READ_LEN):
                if random.random() < BASE_ERROR_RATE:
                    orig = read[j]
                    mut_base = random.choice([b for b in bases if b != orig])
                    read[j] = mut_base
                
            read_str = "".join(read)
            f.write(f">read_{i}\n{read_str}\n")
            
            if (i + 1) % 1000000 == 0:
                print(f"Generated {i + 1} reads...")

if __name__ == '__main__':
    setup_dirs()
    ref = download_and_extract_ref()
    generate_reads_empirical(ref)
    print("Empirical data generation complete!")
