import pandas as pd
from Bio import SeqIO
import pickle
import warnings
warnings.filterwarnings("ignore")

def complete_sequence(ent_file):
    seq = ""
    for record in SeqIO.parse(ent_file, 'pdb-atom'):
        seq+=record.seq
    return seq


print("extracting data...")
FILENAME = snakemake.input[0]

with open(FILENAME, 'rb') as f:
    protein_files = pickle.load(f)

sequences = []
for file in protein_files:
    sequences.append(complete_sequence((r"data/" + file)))

df = pd.DataFrame({"Label": protein_files, "Sequence": sequences})
df.to_csv(snakemake.output[0])
