# collection of helper functions to help fun Oblong for Bigplant
# Phylogeny Pipeline

import dnachisel
import numpy as np
import re
import shutil

from Bio import Phylo
from Bio import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path


def get_backtranslation_table_for_oblong(table='Standard'):
    backtranslation_table = dnachisel.biotools.get_backtranslation_table(
        table_name=table)
    backtranslation_table['-'] = ['---']

    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1976440/
    # Selenocysteine (Sec) and pyrrolysine (Pyl) are rare amino acids that are
    # cotranslationally inserted into proteins and known as the 21st and 22nd
    # amino acids in the genetic code. Sec and Pyl are encoded by UGA and UAG
    # codons, respectively, which normally serve as stop signals.

    # U	Selenocysteine (rare)
    # UGA but oblong don't deal with U
    backtranslation_table['U'] = ['TGA']
    # O	Pyrrolysine (rare)
    # UAG but oblong don't deal with U
    backtranslation_table['O'] = ['TAG']

    # B: [GA]A[TC], D or N; D slightly more common
    backtranslation_table['B'] = backtranslation_table['D'] + \
        backtranslation_table['N']
    # Z: [CG]A[AG], E or Q; E more common
    backtranslation_table['Z'] = backtranslation_table['E'] + \
        backtranslation_table['Q']
    # J: [CA]T[TCA] L or I; have been replacing with L
    backtranslation_table['J'] = backtranslation_table['L'] + \
        backtranslation_table['I']

    # X	any
    backtranslation_table['X'] = ['NNN']

    return backtranslation_table


def reverse_translate(protein_sequence: str,
                      randomize_codons=False,
                      table="Standard"):
    backtranslation_table = get_backtranslation_table_for_oblong(table=table)

    if randomize_codons:
        random_numbers = np.random.randint(0, 1000, len(protein_sequence))
        random_indices = [
            random_number % len(backtranslation_table[aa])
            for aa, random_number in zip(protein_sequence, random_numbers)
        ]
        return "".join([
            backtranslation_table[aa][random_indice]
            for aa, random_indice in zip(protein_sequence, random_indices)
        ])
    return "".join([backtranslation_table[aa][0] for aa in protein_sequence])


# For each aligned family file in TNT_DATA_DIR's subdir,
# reverse translate!
def reverse_translate_families(src_data_dir,
                               dst_data_dir,
                               family_filename='FAMILY.aligned'):
    """Reverse translate protein family files from src data directory and save
    resulting dna family files into dst data dir, preserving the subdirectory
    structures.
    """
    cnt = 0

    for src_child_dir in Path(src_data_dir).iterdir():
        if not Path(src_child_dir).is_dir():
            continue
        family_file = src_child_dir.joinpath(src_child_dir, family_filename)
        # print(f'family_file is {family_file}')
        if not Path(family_file).exists():
            continue
        if Path(src_child_dir.joinpath(family_filename + '.revfasta')).exists(): #don't remake files
            continue
        
        folder_num = src_child_dir.parts[-1]
        dst_child_dir = Path(dst_data_dir).joinpath(folder_num)
        dst_child_dir.mkdir(parents=True, exist_ok=True)
        output_file = dst_child_dir.joinpath(family_filename + '.revfasta')
        # output_handle = open(output_file, 'w')

        seq_records = SeqIO.parse(family_file, 'fasta')
        # generate iterator for seq records so that only one record will be in memory at any one time
        seq_records_rev = (SeqRecord(Seq.Seq(reverse_translate(seq_record.seq)), seq_record.id, seq_record.name, seq_record.description) 
        for seq_record in seq_records)
        SeqIO.write(seq_records_rev, output_file, 'fasta')

        # print(f'{folder_num} done')
        cnt += 1

    return cnt


def copy_tnt_trees(src_data_dir, dst_data_dir, oid_filename='oid.tre'):
    cnt = 0
    for src_child_dir in Path(src_data_dir).iterdir():
        if not Path(src_child_dir).is_dir():
            continue

        folder_num = src_child_dir.parts[-1]
        dst_child_dir = Path(dst_data_dir).joinpath(folder_num)
        dst_child_dir.mkdir(parents=True, exist_ok=True)

        oid_tre = Path(src_child_dir).joinpath(oid_filename)
        if oid_tre.exists():
            shutil.copy(oid_tre, dst_child_dir)
            cnt += 1

    return cnt


def count_families(data_dir,
                   family_filename='FAMILY.aligned',
                   seq_id_pattern=r'>',
                   sizes=[20, 500, 1500, 2000, 2500, 5000, 10000, float('inf')]): #procfiles, TNTB, .., .., .., .., TNTD, rest

    buckets = [[] for _ in range(len(sizes))]

    for family_file in Path(data_dir).glob('*/' + family_filename):
        if_seq_id = ( 1 if line.startswith(">") else 0 for line in open(family_file, "r"))
        count = sum(if_seq_id)
        for i, n in enumerate(sizes):
            if count <= n:
                buckets[i].append(family_file)
                break

    return buckets, sizes


def nexus_to_newick(inputfile, outputfile):
    Phylo.convert(inputfile, 'nexus', outputfile, 'newick')
