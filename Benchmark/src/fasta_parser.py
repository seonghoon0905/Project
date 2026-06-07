import os

def parse_config(config_path):
    # 설정 파일(config.txt)을 읽어 L, D 및 입출력 파일 경로를 반환
    config = {'L_min': 0, 'L_max': 0, 'D': 0, 'ref_path': '', 'reads_path': '', 'output_path': 'Outputs/result.sam'}
    
    with open(config_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(':', 1)
            if len(parts) == 2:
                key = parts[0].strip()
                val = parts[1].strip()
                if key == 'L':
                    # L: 10 또는 L:[32, 100] 등 다양할 수 있음
                    val = val.replace('[', '').replace(']', '').replace(' ', '')
                    bounds = val.split(',')
                    if len(bounds) == 2:
                        config['L_min'] = int(bounds[0])
                        config['L_max'] = int(bounds[1])
                    else:
                        config['L_min'] = int(bounds[0])
                        config['L_max'] = int(bounds[0])
                elif key == 'D':
                    config['D'] = int(val)
                elif key == 'Reference':
                    config['ref_path'] = val
                elif key == 'Reads':
                    config['reads_path'] = val
                elif key == 'Output':
                    config['output_path'] = val
    return config

def read_reference_fasta(ref_path):
    # FASTA 포맷의 레퍼런스 시퀀스를 읽어 단일 문자열로 반환
    seq = ""
    name = ""
    with open(ref_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                name = line[1:].split()[0]
            else:
                seq += line
    return seq, name

def read_reads_fasta(reads_path):
    # FASTA 포맷의 Read 집합을 읽어 이름과 시퀀스의 튜플 리스트로 반환
    reads = []
    current_name = ""
    current_seq = ""
    
    with open(reads_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_name:
                    reads.append((current_name, current_seq))
                current_name = line[1:].split()[0]
                current_seq = ""
            else:
                current_seq += line
        if current_name:
            reads.append((current_name, current_seq))
            
    return reads
