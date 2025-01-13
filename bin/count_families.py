from utils import *
import sys

if len(sys.argv) != 3:
    sys.exit("invalid arguments\nUsage: python count_families.py <oblong_data_dir> <family_cnt>")

OBLONG_DATA_DIR, family_cnt = sys.argv[1], sys.argv[2]

buckets, sizes = count_families(OBLONG_DATA_DIR)

for i, sz in enumerate(sizes):
    print("Number of families <=", sz, ":", len(buckets[i]))

print("Total number of families:", sum([len(b) for b in buckets]))

f = open(family_cnt, 'w')
for i in range(len(buckets)):
    for j in buckets[i]:
        dir = str(j).split('/')[-2]
        f.write(f'{dir},{i}\n')

f.close()
print(f'results written to {family_cnt}')