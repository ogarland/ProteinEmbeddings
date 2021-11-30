import warnings
warnings.filterwarnings("ignore")
import torch
import esm
import pickle
import sys


model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")

#Prepare dataset
data = [
    ("alanine", "A"), ("arginine ", "R"), ("asparagine","N"),("asparagine","N"),
    ("aspartic acid", "D"), ("cysteine", "C"), ("glutamine","Q"),("glutamic acid","E"),
    ("glycine", "G"), ("histidine", "H"), ("isoleucine","I"),("leucine","L"),
    ("lysine ", "K"), ("methionine", "M"), ("phenylalanine","F"),("proline","P"),
    ("serine",  "S"), ("threonine", "T"), ("tryptophan","W"),("tyrosine","Y"),
    ("valine", "V")
]


# Load ESM-1b model
model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
batch_converter = alphabet.get_batch_converter()
batch_labels, batch_strs, batch_tokens = batch_converter(data)

# Extract per-residue representations (on CPU)
with torch.no_grad():
    results = model(batch_tokens, repr_layers=[33], return_contacts=True)
token_representations = results["representations"][33]


# Generate per-sequence representations via averaging
# NOTE: token 0 is always a beginning-of-sequence token, so the first residue is token 1.
sequence_representations = []
for i, (_, seq) in enumerate(data):
    sequence_representations.append(token_representations[i, 1 : len(seq) + 1].mean(0).numpy())

with open('amino_acid_embeddings.pickle', 'wb') as f:
    pickle.dump([data, sequence_representations], f)
