### The netcore environment was made and will be used for saving all data related to our networks as well as being the environment which python works inside it.
Conda activate netcore   
### Two Databases String and ConsensusPATHDB (from Human, only) were downloaded from their corresponding websites. 
### After download (inside the built conda environment netcore), they were decompressed by gzip.

### Receiving database address in conda environment utilizing realpath command
realpath 9606.protein.links.v11.5.txt  ###   '/home/mehdi/9606.protein.links.v11.5.txt'
realpath 9606.protein.info.v11.5.txt  ### ''/home/mehdi/9606.protein.info.v11.5.txt'
realpath ConsensusPathDB_human_PPI   ###    '/home/mehdi/ConsensusPathDB_human_PPI'

### Removing first row from the file and then save the changes by means of nano
nano ConsensusPathDB_human_PPI 

 ### Running python inside netcore environment and importing needed python libraries
python                                                           
import pandas as pd
import numpy as np
import scipy as sp
import networkx as np
import matplotlib.pyplot as plt

### Loading ConsensusPATHDB
conpath=pd.read_csv('/home/mehdi/ConsensusPathDB_human_PPI', sep='\t')
### Checking database
conpath.head() 
type(conpath)

### Choosing 2 columns from it
conpath2=conpath[['interaction_participants__genename','interaction_confidence']]
conpath2.head()
type(conpath2)

### Deleting less than 0.7 rows
conpath3=conpath2[(conpath2.interaction_confidence > 0.7 )]
conpath3.head()
conpath3.interaction_confidence.mean()           ###  0.8921013785832059
conpath3.interaction_confidence.max()              ### 1.0
conpath3.interaction_confidence.min()              ### 0.700003

### Splitting interaction_participants_genename column into two new columns gene1 and gene2 using string_split function
conpath3[['gene1','gene2']]=conpath3.interaction_participants__genename.str.split(pat=',', expand=True)
conpath3.head()
type(conpath3)

### Generating two-column conpath database including gene1 and gene2
conpath4=conpath3[['gene1','gene2']]
conpath4.head()
type(conpath4)

### Building a networkx graph object
conpath=nx.from_pandas_edgelist(conpath4, source='gene1', target='gene2')
type(conpath) 

### Extracting characteristic properties of conpath graph 
print(nx.info(conpath))
#Name:
#Type: Graph
#Number of nodes: 12919
#Number of edges: 258171
#Average degree:  39.9676

### Loading String Database
string=pd.read_csv('/home/mehdi/9606.protein.links.v11.5.txt',sep=' ')
string.head()
type(string)

### Normalizing combined_score between 0 and 1 (and adding another column normalized_score to string
string['normalized_score']=(string['combined_score']-string['combined_score'].min())/(string['combined_score'].max()-string['combined_score'].min())
string.head()

### Deleting less than 0.7 rows
string2=string[(string.normalized_score > 0.7)]
string2.head()
type(string2)

### Generating edge-list dataframe for string database
string3=string2[['protein1','protein2']]
string3.head()
type(string3)

### Importing another dataframe for creating a dictionary from it
dataframe=pd.read_csv('/home/mehdi/9606.protein.info.v11.5.txt', sep='\t')
dataframe.head()
type(dataframe)
### Renaming columns in dataframe
dataframe2=dataframe.rename(columns={'#string_protein_id':'protein_names','preferred_name':'gene_names'})
dataframe2.head()
type(dataframe2)
dataframe3=pd.DataFrame(dataframe2, columns=['protein_names','gene_names'])
dataframe23.head()
type(dataframe3)
### Building a Dictionary from dataframe3
test_keys=list(dataframe3.protein_names)
test_values=list(dataframe3.gene_names)
len(test_keys)  ### 19566
len(test_values) ### 19566
pgdict=dict(zip(test_keys,test_values))

### Utilizing Dictionary for replacing protein names with gene names on string3 dataframe
string4=string3.applymap(lambda x: pgdict[x])
string4.head()
type(string4)
string=string4.rename(columns={'protein1':'gene1','protein2':'gene2'})
string.head()
type(string)

###  Building a networkx graph object from string
string=nx.from_pandas_edgelist(string, source='gene1', target='gene2')
type(string)

### Extracting characteristic properties of string graph 
print(nx.info(string))
#Name:
#Type: Graph
#Number of nodes: 15710
#Number of edges: 213834
#Average degree:  27.2227

### Saving two Databases (in gml format)
nx.write_gml(conpath,'/home/mehdi/conpath.gml')
nx.write_gml(string,'/home/mehdi/string.gml')
