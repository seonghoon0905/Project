import sys, os, time, random, json, subprocess
import matplotlib.pyplot as plt

def main():
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    inputs_dir = os.path.join(root_dir, 'Inputs')
    outputs_dir = os.path.join(root_dir, 'Outputs')
    bench_dir = os.path.join(root_dir, 'Benchmark', 'src')
    proj_dir = os.path.join(root_dir, 'Project', 'src')

    N_LEN = 5000
    M_READS = 100
    L_LEN = 50
    D_MAX = 1

    print(f"Generating N={N_LEN} ref and M={M_READS} reads...")
    bases = ['A', 'C', 'G', 'T']
    ref = "".join(random.choices(bases, k=N_LEN))
    reads = []
    for _ in range(M_READS):
        start = random.randint(0, N_LEN - L_LEN - 1)
        read = list(ref[start:start+L_LEN])
        if random.random() < 0.1:
            mut_idx = random.randint(0, L_LEN - 1)
            read[mut_idx] = random.choice([b for b in bases if b != read[mut_idx]])
        reads.append("".join(read))
        
    data_path = os.path.join(inputs_dir, 'micro_data.json')
    with open(data_path, 'w') as f:
        json.dump({'ref': ref, 'reads': reads, 'D': D_MAX}, f)
        
    bench_res_path = os.path.join(outputs_dir, 'bench_res.json').replace('\\', '/')
    proj_res_path = os.path.join(outputs_dir, 'proj_res.json').replace('\\', '/')
    data_path_str = data_path.replace('\\', '/')

    bench_script = f'''
import time, json
import fmindex
with open("{data_path_str}") as f:
    data = json.load(f)
t0 = time.time()
fm = fmindex.FMCheckpointing(data['ref'])
build_time = (time.time() - t0) * 1000
t0 = time.time()
mapped = 0
for r in data['reads']:
    if fm.search_mismatch(r, data['D']): mapped += 1
search_time = (time.time() - t0) * 1000
with open("{bench_res_path}", 'w') as f:
    json.dump({{"build": build_time, "search": search_time, "mapped": mapped}}, f)
'''

    proj_script = f'''
import time, json
import fmindex
with open("{data_path_str}") as f:
    data = json.load(f)
t0 = time.time()
fm = fmindex.FMCheckpointing(data['ref'], occ_step=64, sa_step=32)
build_time = (time.time() - t0) * 1000
t0 = time.time()
mapped = 0
for r in data['reads']:
    if fm.search_mismatch(r, data['D']): mapped += 1
search_time = (time.time() - t0) * 1000
with open("{proj_res_path}", 'w') as f:
    json.dump({{"build": build_time, "search": search_time, "mapped": mapped}}, f)
'''

    with open(os.path.join(bench_dir, 'micro.py'), 'w') as f: f.write(bench_script)
    with open(os.path.join(proj_dir, 'micro.py'), 'w') as f: f.write(proj_script)

    print("Running Benchmark...")
    subprocess.run([sys.executable, 'micro.py'], cwd=bench_dir)
    print("Running Project...")
    subprocess.run([sys.executable, 'micro.py'], cwd=proj_dir)

    with open(bench_res_path) as f: b_res = json.load(f)
    with open(proj_res_path) as f: p_res = json.load(f)

    print(f"Benchmark: Build={b_res['build']:.2f}ms, Search={b_res['search']:.2f}ms, Mapped={b_res['mapped']}")
    print(f"Project: Build={p_res['build']:.2f}ms, Search={p_res['search']:.2f}ms, Mapped={p_res['mapped']}")

    labels = ['Index Build Time (ms)', 'Mapping Search Time (ms)']
    bench_times = [b_res['build'], b_res['search']]
    proj_times = [p_res['build'], p_res['search']]
    
    x = [0, 1]
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(8, 6))
    rects1 = ax.bar([i - width/2 for i in x], bench_times, width, label='Benchmark', color='#F44336')
    rects2 = ax.bar([i + width/2 for i in x], proj_times, width, label='Project', color='#4CAF50')
    
    ax.set_ylabel('Time (ms) - Log Scale', fontsize=12)
    ax.set_title(f'Micro-scale Time Complexity Comparison (N={N_LEN:,}, M={M_READS:,})', fontsize=14, pad=20)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=11)
    ax.legend()
    ax.set_yscale('log')
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    
    for rect in rects1:
        height = rect.get_height()
        ax.annotate(f'{height:,.1f}', xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=10)
    for rect in rects2:
        height = rect.get_height()
        ax.annotate(f'{height:,.1f}', xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=10)
                    
    plt.tight_layout()
    output_png = os.path.join(outputs_dir, 'micro_benchmark.png')
    plt.savefig(output_png, dpi=300)
    print(f"Graph saved to {output_png}")

if __name__ == '__main__':
    main()
