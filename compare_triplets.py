#!/usr/bin/env python3
# Compare reference and generated triplets

def parse_sections(filepath):
    """Parse pairs.txt format and return (dense_pairs, refine_pairs, triplets)"""
    dense = []
    refine = []
    triplets = []

    with open(filepath, 'r') as f:
        lines = [l.strip() for l in f.readlines()]

    i = 0
    section_idx = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('#') or line == '':
            i += 1
            continue
        if line.isdigit():
            count = int(line)
            section_data = []
            for j in range(i+1, min(i+1+count, len(lines))):
                if lines[j]:
                    section_data.append(lines[j])
            if section_idx == 0:
                dense = section_data
            elif section_idx == 1:
                refine = section_data
            elif section_idx == 2:
                triplets = section_data
            section_idx += 1
            i += count + 1
        else:
            i += 1

    return dense, refine, triplets

# Parse reference
ref_dense, ref_refine, ref_triplets = parse_sections('Datas/1/pairs.txt')
print(f"Reference: {len(ref_dense)} dense, {len(ref_refine)} refine, {len(ref_triplets)} triplets")

# Parse generated
gen_dense, gen_refine, gen_triplets = parse_sections('Datas/1/pairs_generated.txt')
print(f"Generated: {len(gen_dense)} dense, {len(gen_refine)} refine, {len(gen_triplets)} triplets")

# Convert triplets to sets for comparison
def triplet_to_tuple(line):
    parts = sorted(map(int, line.split()))
    return tuple(parts)

ref_tri_set = set(triplet_to_tuple(t) for t in ref_triplets)
gen_tri_set = set(triplet_to_tuple(t) for t in gen_triplets)

matching_tri = ref_tri_set & gen_tri_set
only_ref = ref_tri_set - gen_tri_set
only_gen = gen_tri_set - ref_tri_set

print(f"\nTriplet comparison:")
print(f"  Matching: {len(matching_tri)}")
print(f"  Only in reference: {len(only_ref)}")
print(f"  Only in generated: {len(only_gen)}")
print(f"  Precision: {len(matching_tri)/len(gen_tri_set)*100:.1f}%")
print(f"  Recall: {len(matching_tri)/len(ref_tri_set)*100:.1f}%")

# Convert dense pairs to sets
def pair_to_tuple(line):
    parts = sorted(map(int, line.split()))
    return tuple(parts)

ref_dense_set = set(pair_to_tuple(p) for p in ref_dense)
gen_dense_set = set(pair_to_tuple(p) for p in gen_dense)

matching_dense = ref_dense_set & gen_dense_set
print(f"\nDense pair comparison:")
print(f"  Matching: {len(matching_dense)}")
print(f"  Precision: {len(matching_dense)/len(gen_dense_set)*100:.1f}%")
print(f"  Recall: {len(matching_dense)/len(ref_dense_set)*100:.1f}%")

# Show some triplets only in reference
print("\nSample triplets ONLY in reference:")
for t in list(only_ref)[:5]:
    print(f"  {t}")

print("\nSample triplets ONLY in generated:")
for t in list(only_gen)[:5]:
    print(f"  {t}")
