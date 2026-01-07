#!/usr/bin/env python3
"""
Analyze characteristics of missed bridge triplets to find selection patterns
"""

def parse_triplets(filepath):
    """Parse triplets from pairs.txt"""
    triplets = []
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip() and not line.startswith('#')]

    i = 0
    # Skip dense pairs
    if i < len(lines):
        count = int(lines[i])
        i += 1 + count
    # Skip refine pairs
    if i < len(lines):
        count = int(lines[i])
        i += 1 + count
    # Read triplets
    if i < len(lines):
        count = int(lines[i])
        i += 1
        for _ in range(count):
            if i < len(lines):
                parts = lines[i].split()
                if len(parts) >= 3:
                    triplets.append(tuple(sorted([int(parts[0]), int(parts[1]), int(parts[2])])))
                i += 1
    return triplets

# Load triplets
ref = parse_triplets('Datas/1/pairs.txt')
gen = parse_triplets('Datas/1/pairs_generated.txt')

ref_set = set(ref)
gen_set = set(gen)
common = ref_set & gen_set
ref_only = list(ref_set - gen_set)
gen_only = list(gen_set - ref_set)

print('='*80)
print('  BRIDGE TRIPLET ANALYSIS')
print('='*80)
print(f'\nTotal triplets: Ref={len(ref)}, Gen={len(gen)}, Common={len(common)}')
print(f'Ref-only={len(ref_only)}, Gen-only={len(gen_only)}')

# Focus on bridges (span > 300)
ref_bridges = [t for t in ref if t[2] - t[0] > 300]
gen_bridges = [t for t in gen if t[2] - t[0] > 300]
common_bridges = list(set(ref_bridges) & set(gen_bridges))
missed_bridges = [t for t in ref_only if t[2] - t[0] > 300]

print(f'\nBridge triplets (span>300):')
print(f'  Reference: {len(ref_bridges)}')
print(f'  Generated: {len(gen_bridges)}')
print(f'  Matched: {len(common_bridges)}')
print(f'  Missed: {len(missed_bridges)}')

# Analyze span distribution of missed bridges
print(f'\nMissed bridges span distribution:')
spans = [t[2] - t[0] for t in missed_bridges]
bins = [(300, 500), (500, 700), (700, 1000), (1000, 2000)]
for low, high in bins:
    count = sum(1 for s in spans if low <= s < high)
    print(f'  [{low:4d}, {high:4d}): {count:2d} ({count*100/len(spans):5.1f}%)')

# Print all missed bridges for manual inspection
print(f'\nAll {len(missed_bridges)} missed bridges:')
print('  ID1    ID2    ID3     Span  Gap1  Gap2')
for tri in sorted(missed_bridges, key=lambda t: t[2] - t[0]):
    span = tri[2] - tri[0]
    gap1 = tri[1] - tri[0]
    gap2 = tri[2] - tri[1]
    print(f'  {tri[0]:4d}  {tri[1]:4d}  {tri[2]:4d}  {span:5d}  {gap1:4d}  {gap2:4d}')

# Check if there's a pattern in intermediate gaps
print(f'\nGap pattern analysis for missed bridges:')
gaps1 = [t[1] - t[0] for t in missed_bridges]
gaps2 = [t[2] - t[1] for t in missed_bridges]
all_gaps = gaps1 + gaps2

print(f'  Gap1 (tri[1]-tri[0]): min={min(gaps1)}, max={max(gaps1)}, avg={sum(gaps1)/len(gaps1):.1f}')
print(f'  Gap2 (tri[2]-tri[1]): min={min(gaps2)}, max={max(gaps2)}, avg={sum(gaps2)/len(gaps2):.1f}')
print(f'  All gaps: min={min(all_gaps)}, max={max(all_gaps)}, avg={sum(all_gaps)/len(all_gaps):.1f}')

# Check balance between gaps
print(f'\nGap balance (gap1 vs gap2):')
balanced = sum(1 for i in range(len(gaps1)) if abs(gaps1[i] - gaps2[i]) < 100)
unbalanced = len(gaps1) - balanced
print(f'  Balanced (|gap1-gap2| < 100): {balanced}')
print(f'  Unbalanced: {unbalanced}')

# For reference, compare with matched bridges
if common_bridges:
    print(f'\n\\nFor comparison - MATCHED bridges:')
    matched_gaps1 = [t[1] - t[0] for t in common_bridges]
    matched_gaps2 = [t[2] - t[1] for t in common_bridges]
    print(f'  Gap1: avg={sum(matched_gaps1)/len(matched_gaps1):.1f}')
    print(f'  Gap2: avg={sum(matched_gaps2)/len(matched_gaps2):.1f}')

print('\n' + '='*80)
