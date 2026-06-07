import bwt
import pickle
import os

class FMCheckpointing:
    def __init__(self, text=None, occ_step=64, sa_step=32):
        if text is not None:
            self.occ_step = occ_step
            self.sa_step = sa_step
            self.ref_seq = text
            self.data, sa = bwt.make_bwt(text)
            self.occ = bwt.calc_first_occ(self.data)
            self.C = bwt.calc_checkpoints(self.data, occ_step)
            self.sa_sample = {i: sa[i] for i in range(len(sa)) if sa[i] % sa_step == 0}

    def save_to_file(self, filename):
        data_to_save = {
            'occ_step': self.occ_step,
            'sa_step': self.sa_step,
            'ref_seq': self.ref_seq,
            'data': self.data,
            'occ': self.occ,
            'C': self.C,
            'sa_sample': self.sa_sample
        }
        try:
            with open(filename, 'wb') as f:
                pickle.dump(data_to_save, f)
            return True
        except Exception:
            return False

    def load_from_file(self, filename):
        if not os.path.exists(filename):
            return False
        try:
            with open(filename, 'rb') as f:
                loaded_data = pickle.load(f)
            self.occ_step = loaded_data['occ_step']
            self.sa_step = loaded_data['sa_step']
            self.ref_seq = loaded_data['ref_seq']
            self.data = loaded_data['data']
            self.occ = loaded_data['occ']
            self.C = loaded_data['C']
            self.sa_sample = loaded_data['sa_sample']
            return True
        except Exception:
            return False

    def _count(self, idx, qc):
        return bwt.count_letter_with_checkpoints(self.C, self.occ_step, self.data, idx, qc)

    def _lf(self, idx, qc):
        if qc not in self.occ:
            return 0
        return self.occ[qc] + self._count(idx, qc)

    def _walk(self, idx):
        steps = 0
        current_idx = idx
        
        while current_idx not in self.sa_sample:
            char = self.data[current_idx]
            steps += 1
            current_idx = self._lf(current_idx, char)
            
        return self.sa_sample[current_idx] + steps

    def bounds(self, q):
        top = 0
        bot = len(self.data)
        
        for qc in reversed(q):
            top = self._lf(top, qc)
            bot = self._lf(bot, qc)
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
