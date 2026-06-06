import bwt
import pickle
import os

class FMCheckpointing:
    def __init__(self, text=None, step=50):
        # text가 주어지면 새 인덱스를 빌드합니다.
        if text is not None:
            self.step = step
            self.ref_seq = text
            self.data = bwt.make_bwt(text)
            self.occ = bwt.calc_first_occ(self.data)
            self.C = bwt.calc_checkpoints(self.data, step)
            self.offset_cache = {}

    def save_to_file(self, filename):
        # 생성된 인덱스 구조체를 파일로 직렬화(Serialization)하여 영구 보관합니다.
        # offset_cache는 휘발성 데이터이므로 저장하지 않습니다.
        data_to_save = {
            'step': self.step,
            'ref_seq': self.ref_seq,
            'data': self.data,
            'occ': self.occ,
            'C': self.C
        }
        try:
            with open(filename, 'wb') as f:
                pickle.dump(data_to_save, f)
            return True
        except Exception:
            return False

    def load_from_file(self, filename):
        # 디스크에 저장된 인덱스 파일(.pkl)을 읽어와 즉시 검색 가능한 상태로 복원합니다.
        if not os.path.exists(filename):
            return False
        try:
            with open(filename, 'rb') as f:
                loaded_data = pickle.load(f)
            self.step = loaded_data['step']
            self.ref_seq = loaded_data['ref_seq']
            self.data = loaded_data['data']
            self.occ = loaded_data['occ']
            self.C = loaded_data['C']
            self.offset_cache = {}
            return True
        except Exception:
            return False

    def _count(self, idx, qc):
        # idx까지 qc가 몇 번 등장했는지 카운트합니다.
        return bwt.count_letter_with_checkpoints(self.C, self.step, self.data, idx, qc)

    def _lf(self, idx, qc):
        # LF Mapping: L열의 특정 글자가 F열에서 차지하는 절대 위치를 계산합니다.
        if qc not in self.occ:
            return 0
        return self.occ[qc] + self._count(idx, qc)

    def _walk(self, idx):
        # F열의 행 번호(idx)를 시작으로, 문자열 맨 앞('\0')이 나올 때까지 
        # 거꾸로 걸어가며(Backward Walk) 원본 문자열 상의 진짜 인덱스를 추적합니다.
        r = 0
        i = idx
        
        while self.data[i] != '\0':
            # 캐시에 이미 계산된 위치가 있다면 즉시 리턴 (속도 최적화)
            if i in self.offset_cache:
                r += self.offset_cache[i]
                break
            r += 1
            i = self._lf(i, self.data[i])
            
        # 추적한 발걸음(인덱스)을 캐시에 저장
        if idx not in self.offset_cache:
            self.offset_cache[idx] = r
            
        return r

    def bounds(self, q):
        # 찾고자 하는 쿼리를 뒤에서부터 검색하며 F열 상의 행 범위(top ~ bot)를 좁혀나갑니다.
        top = 0
        bot = len(self.data)
        
        for qc in reversed(q):
            top = self._lf(top, qc)
            bot = self._lf(bot, qc)
            if top == bot:
                return -1, -1 # 범위가 사라지면 탐색 실패
        return top, bot

    def search(self, q):
        # 완전 일치(Exact Match) 탐색을 수행하고, 일치한 모든 원본 인덱스 리스트를 반환합니다.
        top, bot = self.bounds(q)
        if top == -1 and bot == -1:
            return []
            
        matches = []
        for i in range(top, bot):
            matches.append(self._walk(i))
        return sorted(matches)

    def search_mismatch(self, q, max_mismatches):
        # 비둘기집 원리(Indexing Algorithm)를 사용해 최대 허용 미스매치(D) 이내로 일치하는 위치를 찾습니다.
        # D개의 미스매치가 발생하더라도 D+1개의 조각 중 최소 1개는 반드시 완전 일치해야 함을 이용합니다.
        L = len(q)
        segments = max_mismatches + 1
        seg_len = L // segments
        
        found_positions = set()
        results = []
        
        # 1. 쿼리를 조각내기
        for i in range(segments):
            start = i * seg_len
            length = (L - start) if (i == segments - 1) else seg_len
            segment = q[start:start+length]
            
            # 2. 각 조각별로 완전 일치하는 위치들을 검색
            exact_matches = self.search(segment)
            
            for pos in exact_matches:
                # 3. 조각이 일치한 위치를 바탕으로, 전체 쿼리가 시작해야 할 위치를 역산
                read_start_in_ref = pos - start
                
                # 범위를 벗어나면 무시
                if read_start_in_ref < 0 or read_start_in_ref + L > len(self.ref_seq):
                    continue
                    
                # 이미 검사한 위치면 스킵
                if read_start_in_ref in found_positions:
                    continue
                found_positions.add(read_start_in_ref)
                
                # 4. 전체 위치를 알아냈으니 원본 레퍼런스 시퀀스와 1:1 대조하여 실제 미스매치 수 검증
                mismatches = 0
                for j in range(L):
                    if self.ref_seq[read_start_in_ref + j] != q[j]:
                        mismatches += 1
                        if mismatches > max_mismatches:
                            break
                            
                # 최종 검증을 통과했다면 결과 리스트에 추가
                if mismatches <= max_mismatches:
                    results.append({
                        'position': read_start_in_ref + 1, # SAM 포맷은 1-based index 사용
                        'mismatches': mismatches,
                        'cigar': f"{L}M"
                    })
                    
        return sorted(results, key=lambda x: x['position'])
