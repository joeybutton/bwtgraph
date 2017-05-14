# -*- coding: utf-8 -*-
"""
joeybutton and leomcelroy Bioinformatics final project

Usage on with test data: python final.py variants.csv reads-NA12878.sam

Implementation of Smith-Waterman algorithm on a variation graph for demonstrating local alignment of reads to the reference with known variations

Some notes:
--- Topologically sorted graph is just used as a lookup table between our scores matrix and our actual scoring matrix
--- Scoring matrix is just a holding data structure for values in the graph
--- Smith-Waterman algorithm is implemented using adjacencies in the graph - not adjacent cells in the matrix
--- A Header line was added to the SAM file for easy use
--- Variants.csv can be created from any standard VCF, 
"""
from Bio.SubsMat.MatrixInfo import pam250
import argparse
import csv 


REF = "TTGTTTTCTTTTGAAACAGAATCTCACTCTGCAGTCCAGGCTGGAGTGCAGCGGTGCAATCTTGGCTCACTGCAACCTCTGCCTTGTGAGTTCAAGCGATTCTCCTGCCTCAGCCTCTAGACTAGCTGGGATTACAGGTGCATGCCACCATGTCCAGCTAACTTTTTTTGTTTGTTTATTTGTTTGTTTGTTTGTTTTGAGACGGAGTCTTGCTCTATTGCCCAGGCTGGAGTGCAGTGGTGCAATCTCGGCTCACTGCAAGCTTCACCTCCCGGGTTCATGCCATTCTTCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAAGTGCCCGCCACCA"
START = 1520505
END = 1520840
PENALTY = -15

def ref_to_graph(ref):
    """
    ref <Str>: reference genome chunk represented as a string of bases

    ret <List: <Dict, node>>: the graph represented as a series of nodes with pointers to other nodes 
    """
    graph = list()
    for i in range(len(ref)):
        base = REF[i]
        out = [i+1] 
        if out[0] >= len(ref): out = []
        iinn = [i-1]
        if iinn[0] <= -1: iinn = []
        graph.append({
            'base': base,
            'out': out,
            'in': iinn
            })
    return graph

def variants_onto_graph(graph,variants):
    """
    graph <List: <Dict, node>>: our graph represented as a list of nodes with adjacency pointers
    variants <List: <Dict, variants>>: each variant represented as a dictionary with necessary info

    ret: <List, <Dict, node>>: each node including variations now present in the graph
    """
    def vs(var):
        """ Helper func for finding start of variant"""
        return int(var['start'])-START
    def ve(var):
        """ Helper func for finding end of variant"""
        return int(var['end'])-START
    for var in variants:
        if var['type'] == 'm':
            graph.append({
            'base': var['reads'],
            'out': [],
            'in': [],
            'MUT': True
            })

            for inputs in graph[int(var['start'])-START]['in']: # for all nodes into the current
                graph[inputs]['out'].append(len(graph)-1) # for all previous nodes
                graph[-1]['in'].append(inputs)

            for outputs in graph[int(var['start'])-START]['out']:
                graph[outputs]['in'].append(len(graph)-1) # for all next nodes
                graph[-1]['out'].append(outputs)                

        if var['type'] == 'del':
            graph[vs(var)]['del_out'] = ve(var)-vs(var)
            graph[vs(var)]['out'].append(ve(var))

            graph[ve(var)]['del_in'] = [vs(var)]
            graph[ve(var)]['in'].append(vs(var))
    return graph
        

def fit(topo_map, read, penalty, graph_dict):
    """
    topo_map <List>: Topological ordering of the nodes in graph_dict - used as a lookup table
    read <Str>:  Raw string of the read
    penalty <Int>: Penalty for skip
    graph_dict <Dict>: Variation Graph of the reference genome with variations endcoded as edges and nodes (key: nodeID, val: node)

    ret <List<List<Tup>>>: score matrix that encodes best possible score for elements of the graph up to that point
    """
    matrix = [[(0,(-1,-1)) for x in range(len(topo_map)+1)] for y in range(len(read)+1)]
    for i in range(1, len(read)+1):
        for j in range(1, len(topo_map)+1):
            # score matrix following smith-waterman criteria
            # Each cell represents a node in the graph, relative to a location in a read
            scores = []
            for in_edge in graph_dict[topo_map[j-1]]["in"]:
                # For all of the edges into the current node, check potential scores
                a = (read[i-1], graph_dict[topo_map[j-1]]["base"])
                if a not in pam250:
                    a = (graph_dict[topo_map[j-1]]["base"], read[i-1])

                # if a score includes both bases
                scores.append((matrix[i-1][topo_map.index(in_edge)+1][0] + pam250[a], (i-1, topo_map.index(in_edge)+1) ))
                # if a score includes only the graph
                scores.append((matrix[i][topo_map.index(in_edge)+1][0] + penalty, (i, topo_map.index(in_edge)+1)))

            # if a score includes only the specific base from the read, or 0 (Smith-Waterman edge)
            scores.extend([ (matrix[i-1][j][0] + penalty, (i-1,j)) , (0,(-1,-1))])  
            # Store the best score we've found at as a way to reach this particular node, at this point in the read
            matrix[i][j] = max(scores,key=lambda p: p[0])

    return matrix

def backtrack(scores, max_score, topo_map, read, graph_dict):
    """
    scores <List<List<Tup>>>: Score matrix that maps to the graph for traversal
    max_score <Int>:  maximum value in scores
    topo_map <List: <Int>>: Topological ordering of the nodes in reference
    read <Str>: the raw read being fit to the backtrack
    graph_dict <Dict>: Variation Graph of the reference genome with variations endcoded as edges and nodes (key: nodeID, val: node)

    ret <Tup: <Str>, <Str>>: The two alignments as Reference,read
    """
    max_i = [max(j,key=lambda p: p[0]) for j in scores].index(max_score)
    max_j = scores[max_i].index(max_score)
    
    i, j = max_i, max_j
    s1, s2 = '', ''
    
    while (i > 0 and j > 0):
        cur = scores[i][j]
        if cur[0] == 0:
            break

        pred = cur[1]
        # using the annotations from the graph, determine what the traversal looks like
        if pred[0] == i:
            # came from same column
            if 'del_out' in graph_dict[topo_map[j-1]]:
                # The reads follow a known deletion in our variation graph
                s2 = (graph_dict[topo_map[j-1]]['del_out']*'-')+s2
                s1 = (graph_dict[topo_map[j-1]]['del_out']*'-') + s1

            else:
                # Unknown Deletion
                s2 = ((j-pred[1])*'-')+s2
                s1 = str(graph_dict[topo_map[j-1]]['base']) + s1
        elif pred[1] == j:
            # came from same row
            s2 = read[i-1]+s2     
            s1 = "-" + s1
        else:
            # use both strings
            if 'MUT' in graph_dict[topo_map[j-1]]:
                # Follow a known mutation path
                s2 = ' |'+read[i-1]+'| '+s2
                s1 = ' |'+str(graph_dict[topo_map[j-1]]['base'])+'| '+ s1  
            else:
                # This is an unknown mutation/variation         
                s2 = read[i-1]+s2
                s1 = str(graph_dict[topo_map[j-1]]['base']) + s1      
        i,j = pred

    return s1, s2

def graph_to_dict(graph):
    """
    graph <List: <Dict>>: the initial variation graph, converted to a dictionary for traversal

    ret <Dict: <Dict, node>>: the graph represented in a python dictionary
    """
    dic_graph = {}
    for i in range(len(graph)):
        dic_graph[i] = graph[i]

    return dic_graph


def toposort(graph_dict):
    """
    Function for taking in the graph dictionary and returning the topological ordering as a list
    """
    L = []
    marked = []
    toVisit = list(graph_dict.keys())

    def visit(node):
        # recursive helper function
        if node in marked:
            return("NOT DAG ALERT")
        if node in toVisit:
            marked.append(node)
            for entry in sorted(graph_dict[node]["out"])[::-1]:
                visit(entry)
            toVisit.remove(node)
            marked.remove(node)
            L.insert(0,node)

    # Toposort our graph
    while len(toVisit) > 0:
        visit(toVisit[0])
    return L     

def traversals(s,e,graph):
    """
    Debug function for printing chunks of the varian graph
    """
    if s == e:
        return [e]
    travs = []
    print(s)
    print(graph[s])
    for n in graph[s]['out']:
        travs.append(traversals(n,e,graph))
    return travs

def read_file(filename,delim=","):
    """
    simple function for reading in csv/tab-delimited files
    """
    with open(filename, 'rb') as csvfile:
        reader = csv.DictReader(csvfile,delimiter=delim)
        variants = list()
        for row in reader:
            variants.append(row)
        return variants

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='read in input file')
    parser.add_argument('filename', metavar='f', type=str, nargs='+',
                        help='filename of input')
    parser.add_argument('reads', metavar='r', type=str, nargs='+',
                        help='filename of input')

    args = parser.parse_args()
    variants = read_file(args.filename[0])
    reads = read_file(args.reads[0],'\t')

    # construct the reference graph, and add in the known variants
    graph = ref_to_graph(REF)
    graph_varied = variants_onto_graph(graph,variants)
    graph_dict = graph_to_dict(graph)

    # topologically sort the graph
    L = toposort(graph_dict)

    print('REF: '+REF)
    for read in reads:
        # fit each of our reads relative to the variant graph using smith-waterman
        mat = fit(L, read['RAW'], PENALTY, graph_dict)

        # backtrack in our scoring matrix, to generate the fit and show the variations/deletions
        topo,align = backtrack(mat,max([max(l,key=lambda p:p[0]) for l in mat]), L, read['RAW'], graph_dict)
        print((int(read['start'])-START)*" "+"Rad: " +align)



