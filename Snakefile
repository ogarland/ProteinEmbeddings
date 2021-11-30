import pickle
import os 
files = os.listdir("data/")
with open('proteins.pkl', 'wb') as f:
    pickle.dump(files, f)


rule sequence_extraction:
    input:
        "proteins.pkl"
    output:
        temp("results/sequences.csv")
    script:
        "sequence_extraction.py"

rule get_embeddings:
    input:
        "results/sequences.csv"
    output:
        "results/embeddings.pkl"
    script:
        "get_embeddings.py"

rule tSNE:
    input:
        "results/embeddings.pkl",
        "proteins.pkl"
    output:
        "results/tSNE.png"
    script:
        "tSNE.py"