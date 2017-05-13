# -*- coding: utf-8 -*-
from Bio.SubsMat.MatrixInfo import pam250
import argparse
import csv 


REF = "TTGTTTTCTTTTGAAACAGAATCTCACTCTGCAGTCCAGGCTGGAGTGCAGCGGTGCAATCTTGGCTCACTGCAACCTCTGCCTTGTGAGTTCAAGCGATTCTCCTGCCTCAGCCTCTAGACTAGCTGGGATTACAGGTGCATGCCACCATGTCCAGCTAACTTTTTTTGTTTGTTTATTTGTTTGTTTGTTTGTTTTGAGACGGAGTCTTGCTCTATTGCCCAGGCTGGAGTGCAGTGGTGCAATCTCGGCTCACTGCAAGCTTCACCTCCCGGGTTCATGCCATTCTTCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAAGTGCCCGCCACCA"
START = 1520505
END = 1520840

def read_file(filename,delim=","):
    with open(filename, 'rb') as csvfile:
        reader = csv.DictReader(csvfile,delimiter=delim)
        variants = list()
        for row in reader:
            variants.append(row)
        return variants

def ref_to_graph():
    graph = list()
    for i in range(len(REF)):
        base = REF[i]
        graph.append({
            'base': base,
            'out': [i+1],
            'in': [i-1]
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

graph = ref_to_graph()
graph_varied = variants_onto_graph(graph,variants)
graph_dict = graph_to_dict(graph_varied)

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
    matrix = [[0 for x in range(len(graph)+1)] for y in range(len(read)+1)]

    for i in range(0,len(read)+1):
        matrix[i][0] = 0
    for j in range(0,len(graph)+1):
        matrix[0][j] = 0
        
    for i in range(1, len(read)+1):
        for j in range(1, len(graph)+1):
            #score
            scores = []
            for in_edge in graph_dict[graph[j-1]]["in"]:
            	if in_edge == -1:
            		continue
            	#pams
            	a = (read[i-1], graph_dict[in_edge]["base"])
            	if a not in pam250:
                	a = (graph_dict[in_edge]["base"], read[i-1])

                scores.append(matrix[i-1][in_edge] + pam250[a])
                scores.append(matrix[i][in_edge] + penalty)

            scores.extend([matrix[i-1][j], 0])
                     
            matrix[i][j] = max(scores)

    return matrix

print(fit(L, 'GTGCAATCTTGGCTCACTGCAACCTCTGCCTTGTGAGTTCAAGCGATTCTCCTGCCTCAGCCTCTAGACTAGCTGGGATTACAGGTGCATGCCACCATGT', -5, graph_dict))
