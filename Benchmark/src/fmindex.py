import bwt

class FMCheckpointing:
    def __init__(self, text):
        self.ref_seq = text
        self.data = bwt.make_bwt(text)
        self.first_occ = bwt.get_first_occ(self.data)

    def _lf(self, char, idx):
        if char not in self.first_occ:
            return 0
        c_value = self.first_occ[char]
        occ_value = bwt.get_occ_count(self.data, char, idx)
        return c_value + occ_value

    def _walk(self, idx):
        steps = 0
        current_idx = idx
        
        while self.data[current_idx] != '$':
            char = self.data[current_idx]
            steps += 1
            current_idx = self._lf(char, current_idx)
            
        return steps

    def bounds(self, q):
        top = 0
        bot = len(self.data)
        
        for char in reversed(q):
            top = self._lf(char, top)
            bot = self._lf(char, bot)
            if top == bot:
                return -1, -1
        return top, bot

    def search(self, q):
        top, bot = self.bounds(q)
        if top == -1 and bot == -1:
            return []
            
        matches = []
        for i in range(top, bot):
            matches.append(self._walk(i))
        return sorted(matches)

    def search_mismatch(self, q, max_mismatches):
        L = len(q)
        segments = max_mismatches + 1
        seg_len = L // segments
        
        found_positions = set()
        results = []
        
        for i in range(segments):
            start = i * seg_len
            length = (L - start) if (i == segments - 1) else seg_len
            segment = q[start:start+length]
            
            exact_matches = self.search(segment)
            
            for pos in exact_matches:
                read_start_in_ref = pos - start
                
                if read_start_in_ref < 0 or read_start_in_ref + L > len(self.ref_seq):
                    continue
                    
                if read_start_in_ref in found_positions:
                    continue
                found_positions.add(read_start_in_ref)
                
                mismatches = 0
                for j in range(L):
                    if self.ref_seq[read_start_in_ref + j] != q[j]:
                        mismatches += 1
                        if mismatches > max_mismatches:
                            break
                            
                if mismatches <= max_mismatches:
                    results.append({
                        'position': read_start_in_ref + 1,
                        'mismatches': mismatches,
                        'cigar': f"{L}M"
                    })
                    
        return sorted(results, key=lambda x: x['position'])
