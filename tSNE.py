from sklearn.manifold import TSNE
from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import warnings
warnings.filterwarnings("ignore")

print("performing tSNE...")

FILENAME1 = snakemake.input[0]
with open(FILENAME1, 'rb') as f:
    emb = np.array(pickle.load(f))
    
FILENAME2 = snakemake.input[1]
with open(FILENAME2, 'rb') as f:
    labels = np.array(pickle.load(f))

colors = cm.rainbow(np.linspace(0, 1, len(labels)))
tsne_aa = TSNE(n_components=2, random_state=1).fit_transform(emb)
for i, pt in enumerate(tsne_aa):
    plt.scatter(tsne_aa[i,0] , tsne_aa[i,1], color = colors[i], label = labels[i])
plt.legend(loc="upper right")

plt.savefig(snakemake.output[0])
