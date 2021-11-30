# Protein Sequence Embeddings

## Submitted: November 29, 2021.

### Background and Rationale
Over the past decade, natural language processing (NLP) techniques based on self-supervised learning, have used the context of a text to predict missing words in a sentence. Upon further investigation, these models have been found to also learn representations of these words which often capture their meaning. [1] What makes this exciting it that it allows for useful information to be represented within a word’s vector, which is the mapping of a word to a vector of real numbers, which can then be used as a valuable input to a machine learning model. 

With the development of exciting new NLP techniques also allows for these methods to be tested with different types of biological data which is also represented in text form. In 2019, Facebook’s AI Research published the paper ‘Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences’, where they used unsupervised deep learning to train a model on 250 million protein sequences. From this model’s layers representations of the protein sequences can be extracted, which contains information about biological properties. In their work, they found that metric structure in the representation space accords with organizing principles at scales from physicochemical to remote homology. They also found that secondary and tertiary protein structure can be identified in representations. [2]

The model from this paper is publicly available and they provide guidance on how to extract representations for your own protein sequences from their model. This allows for the opportunity to use them as input for new machine learning models for other tasks. Giving a model a vector input which you know contains relationships to biological data gives you a step up in training a successful model.

I began to use their model for my own work but found I was running into issues with the memory requirement of loading both the model, my protein sequences, and then also storing the representations for the sequences. I found this problem occurred because in their code, they require you to input each protein into the model. But really what is happening is that each amino acid in the sequence is fed separately and then the output tensor is averaged over the entire sequence.  This made me realize that you really don’t need to input the entire sequence to the model, but rather you can just get the representation of the individual amino acids, map them to your own sequence and find the average representation for the entire sequence.

This workflow was built to allow users to create a labelled training set for their protein data sets, using the representations for learned from the model, without actually needing to load the model at all. I have stored the representations per amino acid in a pickle file which serialized the amino acid tensor structures. The workflow is as follows:

- Sequence Extraction: Protein sequences are extracted from pdb files
- Embedding Mapping: Representations for each sequence are mapped and calculated based on the stored amino acids. They are labelled with their PDB ID and are output to a file which makes it easy to later import for training of an ML model.
- t-SNE Plotting: Proteins embeddings are mapped using the dimentionality reduction method t-SNE. This allows users to investigate potential trends captured in the representations for data exploration purposes.

<img src="https://github.com/ogarland/ProteinEmbeddings/blob/main/workflow.png?raw=true" width="300" height="300">


### Usage
The execution of the following files *requires* a **Linux OS** and **MiniConda for Python 3**. If you are using Windows, please set up using a Virtual Machine.

1. Clone this file in a proper directory. This will download all the files and datasets necessary for execution. The output will be a folder named "ProteinEmbeddings".
```
git clone https://github.com/ogarland/ProteinEmbeddings.git
```
2. Create an environment with all the dependencies. This will manage all the packages under the same collection for execution.
```
conda env create --name protemb --file environment.yml
```
3. Activate the environment.
```
conda activate protemb
```
4. Now run pipeline by executing the following code. This will create two files in the results folder: `embeddings.pkl`, `tSNE.png`.
```
snakemake --cores 1 results/tSNE.png
```


### Input
ProteinEmbeddings takes PDB .ent files as input. Any number of files can be used, and they all should be placed in thefolder `data/`. 

These .ent files can be retrived from the Protein Data Bank (PDB) and contain sequence, experimental and structural information.

### Output
ProteinEmbeddings outputs will appear in the `results/` folder. The first output is a pickle file called `embeddings.pkl` which contains the labels for the proteins as well as their corresponding embedding. The second output, is t-SNE plot of the embeddings.

t-SNE Output example
------------- 
<img src="https://github.com/ogarland/ProteinEmbeddings/blob/main/results/tSNE.png?raw=true" width="500" height="300"> 

- t-SNE plot: This is a visualization for dimension reduction of high dimensional data based on Stochastic Neighbor Embedding. Axes do not refer to spatial coordinates. As the paper describing this model suggests that biological properties are captured in the embeddings, using a tSNE plot allows the users to investigate the kinds of trends that may be captured in their embeddings. 

### References

[1] Collobert, R., &amp; Weston, J. (2008). A unified architecture for natural language processing. Proceedings of the 25th International Conference on Machine Learning - ICML '08. https://doi.org/10.1145/1390156.1390177 

[2] Rives, A., Meier, J., Sercu, T., Goyal, S., Lin, Z., Liu, J., Guo, D., Ott, M., Zitnick, C. L., Ma, J., &amp; Fergus, R. (2019). Biological structure and function emerge from scaling unsupervised learning to 250 million protein sequences. BioRxiv. https://doi.org/10.1101/622803 
