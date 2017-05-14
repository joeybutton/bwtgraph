# Smith Waterman on A variant Graph (SWAG)
#### By: leomcelroy and joeybutton
A project created for CS321 - Bioinformatics Algorithms

##### What is it?
An implementation of a graph-based Smith-Waterman algorithm for aligning reads to a reference genome, when the reference is stored as a graph with variants.

##### How does it work?
By first constructing a graph of the reference, the scoring portion of Smith-Waterman creates the potential edges in the graph. 
Known variants are stored as nodes with annotations on them, so that they can be referenced when reconstructing the path after scoring of the paths has completed.
The algorithm stores the maximum possible alignment score up to each given node, while using up some portion of the read it is aligning. 
When backtracking, the annotated nodes are used to create portions of the path that differ from the reference, but match known variants.

##### How can I test it?
``` 
python swag.py variants.csv reads-NA12878.sam
```


