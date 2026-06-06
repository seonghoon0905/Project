def make_bwt(text):
    # 1. 문자열 끝에 종료 마커 달기
    text = text + "$"
    
    # 2. 순환 문자열(Rotations) 리스트 만들기
    rotations = [text[i:] + text[:i] for i in range(len(text))]
    
    # 3. 알파벳 순서로 정렬하기
    table = sorted(rotations)
    
    # 4. 정렬된 테이블의 마지막 글자만 똑똑 떼어내기
    last_column = [row[-1] for row in table]
    
    # 5. 리스트에 담긴 글자들을 하나의 문자열로 합쳐서 반환
    bwt_result = "".join(last_column)
    
    return bwt_result

# C[c]
def get_first_occ(bwt_str):
    # 1. BWT 문자열에 각 글자가 몇 개씩 있는지 셉니다.
    counts = {}
    for char in bwt_str:
        # 딕셔너리에 글자가 없으면 0부터 시작하고, 1을 더합니다.
        counts[char] = counts.get(char, 0) + 1
        
    # 2. 알파벳 순서대로 정렬합니다. (예: '$', 'a', 'c', 'g')
    sorted_chars = sorted(counts.keys())
    
    # 3. 누적합을 이용해 시작 인덱스를 계산합니다.
    first_occ = {}
    current_index = 0
    
    for char in sorted_chars:
        first_occ[char] = current_index  # 현재 위치를 시작점으로 기록!
        current_index += counts[char]    # 다음 글자는 내 개수만큼 뒤에서 시작함
        
    return first_occ

# 각 글자가 몇번째 인덱스에 처음 등장하는가?

def get_occ_count(bwt_str, char, idx):
    count = 0
    # 0번 인덱스부터 내 위치(idx) 바로 앞까지 훑어봅니다.
    for i in range(idx):
        if bwt_str[i] == char:
            count += 1
    return count

# 특정 문자가 idx까지 총 몇번 나왔을까?

def get_lf_mapping(bwt_str, first_occ, char, idx):
    """ 마법의 공식: LF(i, c) = C[c] + Occ(c, i) """
    c_value = first_occ[char]                       # 1. 시작점
    occ_value = get_occ_count(bwt_str, char, idx)   # 2. 내 앞의 개수
    return c_value + occ_value                      # 3. 다음 위치 반환!

def inverse_bwt(bwt_str):
    """ BWT 문자열을 원래 문자열로 완벽히 복원하는 함수 """
    # 1. 전체 지도의 기준이 될 이정표(C[c])를 만듭니다.
    first_occ = get_first_occ(bwt_str)
    
    # 2. 복원된 글자들을 담을 빈 리스트를 준비합니다. (종료 마커 '$'는 제외)
    original_len = len(bwt_str) - 1
    result = [''] * original_len
    
    # 3. 길 찾기 시작점: 항상 0번 인덱스부터 시작!
    current_idx = 0
    
    # 4. 뒤에서부터 앞으로 거꾸로 글자를 채워 넣습니다.
    for k in range(original_len - 1, -1, -1):
        char = bwt_str[current_idx]  # 현재 위치(L열)의 글자를 꺼냅니다.
        result[k] = char             # 결과 리스트의 빈칸에 뒤에서부터 채워 넣습니다.
        # ★ 핵심: LF-Mapping을 사용해 다음 글자가 있는 위치로 텔레포트!
        current_idx = get_lf_mapping(bwt_str, first_occ, char, current_idx)
        
    return "".join(result)

def get_bounds(bwt_str, first_occ, query):
    """ 쿼리(query)가 존재하는 F열의 시작점(top)과 끝점(bot)을 찾는 함수 """
    top = 0
    bot = len(bwt_str)  # 처음엔 전체 범위를 잡고 시작합니다.
    
    # 파이썬 문법 [::-1]을 사용해 쿼리를 거꾸로 뒤집어서 한 글자씩 꺼냅니다.
    # 예: 'caa' -> 'a', 'a', 'c' 순서로 반복
    for char in query[::-1]:
        # 현재 범위(top ~ bot) 안에서, 찾고자 하는 글자를 따라 다음 범위로 점프!
        top = get_lf_mapping(bwt_str, first_occ, char, top)
        bot = get_lf_mapping(bwt_str, first_occ, char, bot)
        
        print(f"[{char}] 검색 후 범위 -> top: {top}, bot: {bot}") # 진행 상황 확인용
        
        # 만약 top과 bot이 같아져서 범위가 사라졌다면? 매칭 실패!
        if top == bot:
            return 0, 0 
            
    return top, bot

def walk(bwt_str, first_occ, idx):
    """ F열의 인덱스(idx)를 원본 텍스트의 인덱스로 변환하는 함수 """
    steps = 0  # 걸어간 발걸음 수
    current_idx = idx
    
    # 내 바로 앞글자(L열의 글자)가 마커 '$'가 아닐 때까지 계속 왼쪽으로 걸어갑니다!
    while bwt_str[current_idx] != '$':
        char = bwt_str[current_idx]  # 내 바로 앞글자를 확인
        steps += 1                   # 발걸음 수 1 증가
        
        # LF-Mapping을 이용해 앞글자의 위치로 점프!
        current_idx = get_lf_mapping(bwt_str, first_occ, char, current_idx)
        
    return steps  # 마커를 만나면 여태 걸은 발걸음 수를 반환


# 테스트
my_text = "acaacg"
my_text = make_bwt(my_text);
print(my_text)

original_text = inverse_bwt(my_text)
print("복원된 문자열:", original_text)


first_occ = get_first_occ(my_text)
query = "acg"
top, bot = get_bounds(my_text, first_occ, query)
print(f"top:{top}, bot:{bot}, query :{query}")
original_index = walk(my_text, first_occ, top)

print(f"{query}는 원본 텍스트의 {original_index}번 인덱스")