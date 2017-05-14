# -*- coding: utf-8 -*-
from Bio.SubsMat.MatrixInfo import pam250
import argparse
import csv 


REF = "TTGTTTTCTTTTGAAACAGAATCTCACTCTGCAGTCCAGGCTGGAGTGCAGCGGTGCAATCTTGGCTCACTGCAACCTCTGCCTTGTGAGTTCAAGCGATTCTCCTGCCTCAGCCTCTAGACTAGCTGGGATTACAGGTGCATGCCACCATGTCCAGCTAACTTTTTTTGTTTGTTTATTTGTTTGTTTGTTTGTTTTGAGACGGAGTCTTGCTCTATTGCCCAGGCTGGAGTGCAGTGGTGCAATCTCGGCTCACTGCAAGCTTCACCTCCCGGGTTCATGCCATTCTTCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAAGTGCCCGCCACCA"
# REF = "PENALTY"
START = 1520505
END = 1520840
PENALTY = -5

def read_file(filename,delim=","):
    with open(filename, 'rb') as csvfile:
        reader = csv.DictReader(csvfile,delimiter=delim)
        variants = list()
        for row in reader:
            variants.append(row)
        return variants

def ref_to_graph(ref):
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

def vs(var):
    return int(var['start'])-START
def ve(var):
    return int(var['end'])-START

def variants_onto_graph(graph,variants):
    for var in variants:
        # print graph[int(var['start'])-1520505], var
        if var['type'] == 'm':
            graph.append({
            'base': var['reads'],
            'out': [],
            'in': []
            })

            for inputs in graph[int(var['start'])-START]['in']: # for all nodes into the current
                graph[inputs]['out'].append(len(graph)-1) # for all previous nodes
                graph[-1]['in'].append(inputs)

            for outputs in graph[int(var['start'])-START]['out']:
                graph[outputs]['in'].append(len(graph)-1) # for all next nodes
                graph[-1]['out'].append(outputs)                

        if var['type'] == 'del':
            graph[vs(var)-1]['out'].append(ve(var)-1)
            graph[ve(var)-1]['in'].append(vs(var)-1)
    return graph
    # print(graph)
        
def traversals(s,e,graph):
    if s == e:
        return [e]
    travs = []
    print(s)
    print(graph[s])
    for n in graph[s]['out']:
        travs.append(traversals(n,e,graph))
    return travs

parser = argparse.ArgumentParser(description='read in input file')
parser.add_argument('filename', metavar='f', type=str, nargs='+',
                    help='filename of input')
parser.add_argument('reads', metavar='r', type=str, nargs='+',
                    help='filename of input')

args = parser.parse_args()
variants = read_file(args.filename[0])
reads = read_file(args.reads[0],'\t')

def graph_to_dict(graph):
    dic_graph = {}
    for i in range(len(graph)):
        dic_graph[i] = graph[i]

    return dic_graph

graph = ref_to_graph(REF)
# graph_varied = variants_onto_graph(graph,variants)
graph_dict = graph_to_dict(graph)

L = []
marked = []
toVisit = list(graph_dict.keys())

def visit(node):
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

def fit(graph, read, penalty, graph_dict):
    matrix = [[(0,(-1,-1)) for x in range(len(graph)+1)] for y in range(len(read)+1)]
    for i in range(1, len(read)+1):
        for j in range(1, len(graph)+1):
            #score
            scores = []
            for in_edge in graph_dict[graph[j-1]]["in"]:
                #pams
                a = (read[i-1], graph_dict[graph[j-1]]["base"])
                if a not in pam250:
                    a = (graph_dict[graph[j-1]]["base"], read[i-1])

                scores.append((matrix[i-1][graph.index(in_edge)+1][0] + pam250[a], (i-1, graph.index(in_edge)+1) ))
                scores.append((matrix[i][graph.index(in_edge)+1][0] + penalty, (i, graph.index(in_edge)+1)))

            scores.extend([ (matrix[i-1][j][0] + penalty, (i-1,j)) , (0,(-1,-1))])  

            matrix[i][j] = max(scores,key=lambda p: p[0])

    return matrix

def backtrack(scores, max_score, penalty, topo, read, graph_dict):
    max_i = [max(j,key=lambda p: p[0]) for j in scores].index(max_score)
    max_j = scores[max_i].index(max_score)
    # print(max_i,max_j)
    # print([i[:max_j] for i in scores[:max_i]])
    i, j = max_i, max_j
    s1, s2 = '', ''
    # print(scores)

    # print([l[max_j-10:max_j] for l in scores[max_i-10:max_i] ])
    
    while (i > 0 and j > 0):
        cur = scores[i][j]
        if cur[0] == 0:
            break

        pred = cur[1]
        if pred[0] == i:
            # came from same column
            s2 = "-"+s2     
            s1 = str(graph_dict[topo[j-1]]['base']) + s1
        elif pred[1] == j:
            # came from same row
            s2 = read[i-1]+s2     
            s1 = "-" + s1
        else:
            # use both
            s2 = read[i-1]+s2
            s1 = str(graph_dict[topo[j-1]]['base']) + s1
        
        i,j = pred
    return s1, s2
            
reads2 = [
    {"RAW":'MEANLY'}
    ]
for read in reads:
    mat = fit(L, read['RAW'], PENALTY, graph_dict)
    # for row in mat:
    #     print([l[0] for l in row])
    topo,align = backtrack(mat,max([max(l,key=lambda p:p[0]) for l in mat]), PENALTY, L, read['RAW'], graph_dict)
    print(topo)
    print(align)
    print('∆'*10)
