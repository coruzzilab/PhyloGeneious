from utils import *
import sys

if len(sys.argv) != 3:
    sys.exit("invalid arguments\nUsage: python reverse_translate_families.py <input_data_dir> <oblong_data_dir>")

INPUT_DATA_DIR, OBLONG_DATA_DIR = sys.argv[1], sys.argv[2]

num_converted = reverse_translate_families(INPUT_DATA_DIR, OBLONG_DATA_DIR)

print("Number of families converted:", num_converted)
