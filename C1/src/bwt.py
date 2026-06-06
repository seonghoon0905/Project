# BWT 생성 및 LF Mapping에 필요한 핵심 도우미 함수들

def make_bwt(s):
    # Suffix Array를 이용해 BWT를 O(N log N)으로 빠르게 생성합니다.
    s += '\0' # 종료 마커 추가
    
    # 각 인덱스부터 시작하는 접미사의 시작 위치를 저장한 후 사전순으로 정렬
    rotations = sorted(range(len(s)), key=lambda i: s[i:])
    
    # 정렬된 접미사들의 바로 앞 글자 추출 (LF Mapping의 L열)
    r = []
    for i in rotations:
        if i == 0:
            r.append('\0')
        else:
            r.append(s[i-1])
    return "".join(r)

def calc_first_occ(s):
    # BWT 문자열(L열)을 통해 가상의 정렬된 F열에서 각 알파벳이 시작하는 위치를 계산합니다.
    A = {}
    for c in s:
        A[c] = A.get(c, 0) + 1
        
    letters = sorted(A.keys())
    occ = {}
    idx = 0
    for c in letters:
        occ[c] = idx
        idx += A[c]
    return occ

def calc_checkpoints(s, step):
    # 메모리 절약을 위해 step 간격으로만 각 알파벳의 누적 등장 횟수를 저장합니다.
    A = {}
    C = []
    for i, c in enumerate(s):
        if i % step == 0:
            C.append(A.copy())
        A[c] = A.get(c, 0) + 1
    return C

def count_letter_with_checkpoints(C, step, data, idx, qc):
    # 체크포인트 배열 C와 BWT 원본 문자열을 사용하여 특정 위치(idx)까지 문자 qc의 누적 등장 횟수를 계산합니다.
    if idx == 0:
        return 0
        
    # 가장 가까운 과거 체크포인트 인덱스 계산
    chk_idx = idx // step
    c = C[chk_idx].get(qc, 0)
    
    # 체크포인트 위치부터 현재 위치 바로 앞까지 직접 세기
    start_pos = chk_idx * step
    for i in range(start_pos, idx):
        if data[i] == qc:
            c += 1
            
    return c
