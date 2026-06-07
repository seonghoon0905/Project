def make_bwt(text):
    # 1. 문자열 끝에 종료 마커 추가
    text = text + "$"
    
    # 2. 순환 문자열 리스트 
    rotations = [text[i:] + text[:i] for i in range(len(text))]
    
    # 3. 알파벳 순서로 정렬
    table = sorted(rotations)
    
    # 4. 정렬된 테이블의 마지막 글자만 떼어냄
    last_column = [row[-1] for row in table]
    
    # 5. 리스트에 담긴 글자들을 하나의 문자열로 합쳐서 반환
    bwt_result = "".join(last_column)
    
    return bwt_result

def get_first_occ(bwt_str):
    # BWT 문자열에 각 글자가 몇 개씩 있는지 확인
    counts = {}
    for char in bwt_str:
        counts[char] = counts.get(char, 0) + 1
        
    # 알파벳 순서대로 정렬
    sorted_chars = sorted(counts.keys())
    
    first_occ = {}
    current_index = 0
    
    for char in sorted_chars:
        first_occ[char] = current_index 
        current_index += counts[char] 
        
    return first_occ

def get_occ_count(bwt_str, char, idx):
    count = 0
    for i in range(idx):
        if bwt_str[i] == char:
            count += 1
    return count
