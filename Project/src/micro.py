
import time, json
import fmindex
with open("C:/Users/Seonghoon/Desktop/Study/Algorithm/Project/Inputs/micro_data.json") as f:
    data = json.load(f)
t0 = time.time()
fm = fmindex.FMCheckpointing(data['ref'], occ_step=64, sa_step=32)
build_time = (time.time() - t0) * 1000
t0 = time.time()
mapped = 0
for r in data['reads']:
    if fm.search_mismatch(r, data['D']): mapped += 1
search_time = (time.time() - t0) * 1000
with open("C:/Users/Seonghoon/Desktop/Study/Algorithm/Project/Outputs/proj_res.json", 'w') as f:
    json.dump({"build": build_time, "search": search_time, "mapped": mapped}, f)
