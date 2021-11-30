import pickle
import sys
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")


def get_aa_embedding(sequences):
    
    #Open saved embeddings for individual amino acids
    with open(r"amino_acid_embeddings.pickle", "rb") as f:
        emb = pickle.load(f)
    df = pd.DataFrame(emb[0], columns = ["Name", "X"])
    df["Tensor"] = emb[1]
    
    embeddings = []
    
    #Get embedding respective to protein as avg of amino acid
    for sequence in sequences:
        n = len(sequence)
        embedding = np.zeros(shape = (df["Tensor"].values[0]).shape)
        for i in sequence:
            if i == "X": 
                n-=1
                continue
            row = df[df["X"]==i]
            embedding += (row["Tensor"].values[0])
        embeddings.append(embedding/n)
            
    return embeddings 


print("calculating embeddings...")

FILENAME = snakemake.input[0]

df = pd.read_csv(FILENAME)

sequences = df["Sequence"].values
emb = get_aa_embedding(sequences)

    
with open(snakemake.output[0], 'wb') as f:
    pickle.dump(emb, f)
