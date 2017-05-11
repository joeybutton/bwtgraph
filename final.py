# -*- coding: utf-8 -*-
from Bio.SubsMat.MatrixInfo import blosum62
import argparse
import csv 


REF = "TTGTTTTCTTTTGAAACAGAATCTCACTCTGCAGTCCAGGCTGGAGTGCAGCGGTGCAATCTTGGCTCACTGCAACCTCTGCCTTGTGAGTTCAAGCGATTCTCCTGCCTCAGCCTCTAGACTAGCTGGGATTACAGGTGCATGCCACCATGTCCAGCTAACTTTTTTTGTTTGTTTATTTGTTTGTTTGTTTGTTTTGAGACGGAGTCTTGCTCTATTGCCCAGGCTGGAGTGCAGTGGTGCAATCTCGGCTCACTGCAAGCTTCACCTCCCGGGTTCATGCCATTCTTCTGCCTCAGCCTCCCAAGTAGCTGGGACTACAAGTGCCCGCCACCA"
START = 1520505
END = 1520840

def read_file(filename):
    with open(filename, 'rb') as csvfile:
        reader = csv.DictReader(csvfile)
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
    print(graph[s])
    for n in graph[s]['out']:
        travs.append(traversals(n,e,graph))
    return travs




parser = argparse.ArgumentParser(description='read in input file')
parser.add_argument('filename', metavar='f', type=str, nargs='+',
                    help='filename of input')

args = parser.parse_args()
variants = read_file(args.filename[0])


graph = ref_to_graph()
graph_varied = variants_onto_graph(graph,variants)

traversals(170,176,graph_varied)
