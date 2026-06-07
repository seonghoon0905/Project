import PySAIS

def build_suffix_array(text):
    return PySAIS.sais(text)

def make_bwt(text):
    text += "$"
    sa = build_suffix_array(text)
    bwt_chars = []
    
    for pos in sa:
        if pos == 0:
            bwt_chars.append("$")
        else:
            bwt_chars.append(text[pos - 1])
            
    bwt_result = "".join(bwt_chars)
    return bwt_result, sa

def calc_first_occ(s):
    # F열 글자 시작 위치 찾기
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
    # 메모리 절약을 위해 일정 간격으로만 저장
    A = {}
    C = []
    for i, c in enumerate(s):
        if i % step == 0:
            C.append(A.copy())
        A[c] = A.get(c, 0) + 1
    return C

def count_letter_with_checkpoints(C, step, data, idx, qc):
    if idx == 0:
        return 0
        
    chk_idx = idx // step
    c = C[chk_idx].get(qc, 0)
    
    start_pos = chk_idx * step
    for i in range(start_pos, idx):
        if data[i] == qc:
            c += 1
            
    return c
