import sys
import time
from fasta_parser import parse_config, read_reference_fasta, read_reads_fasta
from fmindex import FMCheckpointing

def main():
    config_path = "../Inputs/config.txt"
    if len(sys.argv) > 1:
        config_path = sys.argv[1]

    print(f"Loading config from: {config_path}")
    config = parse_config(config_path)
    
    if not config['ref_path'] or not config['reads_path'] or not config['output_path']:
        print("Config parsing failed. Please check config.txt")
        sys.exit(1)

    # 1. 레퍼런스 게놈 로드
    reference, ref_name = read_reference_fasta(config['ref_path'])
    if not reference:
        print("Reference genome is empty.")
        sys.exit(1)
    if not ref_name:
        ref_name = "unknown_reference"

    print(f"Reference Genome ({ref_name}) length: {len(reference)}")

    # 2. 인덱스 빌드 또는 로드 (Serialization)
    index_path = "../Inputs/index.pkl"
    fm = FMCheckpointing()
    
    start_time = time.time()
    
    if fm.load_from_file(index_path):
        print(f"Successfully loaded FM-Index from {index_path}")
        elapsed = (time.time() - start_time) * 1000
        print(f"FM-Index loaded in {elapsed:.0f} ms.")
    else:
        print("Index file not found. Building FM-Index...")
        fm = FMCheckpointing(reference, occ_step=64, sa_step=32)
        elapsed = (time.time() - start_time) * 1000
        print(f"FM-Index built in {elapsed:.0f} ms.")
        
        if fm.save_to_file(index_path):
            print(f"FM-Index saved to {index_path}")
        else:
            print(f"Failed to save FM-Index to {index_path}")

    # 3. Reads 로드
    print("Loading reads...")
    reads = read_reads_fasta(config['reads_path'])
    print(f"Loaded {len(reads)} reads.")

    # 4. 매핑 수행 및 SAM 파일 출력
    print(f"Mapping reads with max {config['D']} mismatches...")
    start_time = time.time()
    mapped_count = 0
    
    try:
        with open(config['output_path'], 'w', encoding='utf-8') as out_sam:
            # SAM 헤더 작성
            out_sam.write("@HD\tVN:1.0\tSO:unsorted\n")
            out_sam.write(f"@SQ\tSN:{ref_name}\tLN:{len(reference)}\n")
            
            for read_name, read_seq in reads:
                # 각 Read에 대해 미스매치 탐색 수행
                results = fm.search_mismatch(read_seq, config['D'])
                
                if not results:
                    # 매핑 실패 (Flag = 4)
                    out_sam.write(f"{read_name}\t4\t*\t0\t0\t*\t*\t0\t0\t{read_seq}\t*\n")
                else:
                    # 매핑 성공 (첫 번째로 발견된 최상단 결과 리포트)
                    best = results[0]
                    flag = 0 # 정방향 매치
                    out_sam.write(f"{read_name}\t{flag}\t{ref_name}\t{best['position']}\t255\t{best['cigar']}\t*\t0\t0\t{read_seq}\t*\tNM:i:{best['mismatches']}\n")
                    mapped_count += 1
    except Exception as e:
        print(f"Failed to write SAM output: {e}")
        sys.exit(1)

    elapsed = (time.time() - start_time) * 1000
    print(f"Mapping completed in {elapsed:.0f} ms.")
    print(f"Mapped {mapped_count} / {len(reads)} reads successfully.")
    print(f"SAM output saved to: {config['output_path']}")

if __name__ == "__main__":
    main()
