import os
import matplotlib.pyplot as plt

sam_file = "Outputs/result.sam"
exact_match = 0
one_mismatch = 0
unmapped = 0

print("Parsing SAM file...")
with open(sam_file, 'r', encoding='utf-8') as f:
    for line in f:
        if line.startswith('@'):
            continue
        parts = line.split('\t')
        flag = int(parts[1])
        if flag == 4:
            unmapped += 1
        else:
            nm = 0
            for p in parts[11:]:
                if p.startswith('NM:i:'):
                    nm = int(p.split(':')[2])
                    break
            if nm == 0:
                exact_match += 1
            elif nm == 1:
                one_mismatch += 1

total = exact_match + one_mismatch + unmapped

# Plot Bar Chart
labels = ['Exact Match\n(0 Mismatch)', '1 Mismatch', 'Unmapped\n(2+ Mismatches)']
counts = [exact_match, one_mismatch, unmapped]
colors = ['#4CAF50', '#2196F3', '#F44336']

plt.figure(figsize=(8, 6))
bars = plt.bar(labels, counts, color=colors, width=0.6)

# Add values on top of bars
for i, bar in enumerate(bars):
    count = counts[i]
    pct = (count / total) * 100
    plt.text(bar.get_x() + bar.get_width()/2., count + (total*0.01),
             f'{count:,}\n({pct:.2f}%)', 
             ha='center', va='bottom', fontsize=11)

plt.title('Read Mapping Results vs Error Distribution (N=5,000,000)', fontsize=14, pad=20)
plt.ylabel('Number of Reads', fontsize=12)
plt.ylim(0, max(counts) * 1.2) # Add some headroom for text
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()

output_png = 'Outputs/mapping_distribution.png'
plt.savefig(output_png, dpi=300)
print(f"Graph saved to {output_png}")
print(f"Results: Exact: {exact_match}, 1MM: {one_mismatch}, Unmapped: {unmapped}")
