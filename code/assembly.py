def substr(st, start, length):
    return st[start: start + length]

def assign_verteces(sequences, k):
    verteces = []
    seen = {}
    for seq in sequences:
        for index, letter in enumerate(seq):
            sub = substr(seq, index, k)
            if len(sub) < k:
                break
            if not seen.get(sub):
                seen[sub] = True
                verteces.append(sub)
    return verteces

def assign_edges(sequences, k):
    edge_defined = {}
    out_edges = {}
    for seq in sequences:
        for index, letter in enumerate(seq):
            sub1 = substr(seq, index, k)
            sub2 = substr(seq, index + 1, k)
            if len(sub2) < k:
                break
            if not edge_defined.get((sub1, sub2)):
                edge_defined[(sub1, sub2)] = True
                if out_edges.get(sub1):
                    out_edges[sub1].append(sub2)
                else:
                    out_edges[sub1] = [sub2]
    return out_edges

def make_cycle(start_vertex, out_edges):
    vertex = start_vertex
    path = [vertex]
    edges = out_edges.get(vertex)
    while edges:
        if len(edges) == 0:
            break
        vertex = edges.pop()
        path.append(vertex)
        edges = out_edges.get(vertex)
    return path
        

def reconstruct_sequence(path):
    sequence = path[0]
    kmer_len = len(sequence)
    for seq in path[1:]:
        sequence = sequence + substr(seq, kmer_len-1, 1)
    return sequence

def make_contig(vertex, out_edges, verteces):
    cycle = make_cycle(verteces[vertex], out_edges)
    path_growth = 1
    while path_growth > 0:
        size0 = len(cycle)
        for index, seq in enumerate(cycle):
            vertex = cycle[index]
            edges = out_edges.get(vertex)
            if not edges:
                continue
            if len(edges) > 0:
                new_cycle = make_cycle(vertex, out_edges)
                cycle = cycle[:c+1] + new_cycle + cycle[c+1:]
        path_growth = len(cycle) - size0
    contig = reconstruct_sequence(cycle)
    return contig

#k-mer length
k = 10

#counters
c = 0
ct = 1
contig_number = 1


import sys
filename = sys.argv[0]

# open sequences file
input = open('seq_frags/make_seq.out.txt' , 'r')

full_seq = input.readline()

sequences = []
for line in input:
    sequences.append(line.strip())

sequences[0]= full_seq[0:100]


print('Original sequence:')
print(f'{full_seq}')


verteces = assign_verteces(sequences, k)

out_edges = assign_edges(sequences, k)

contig = make_contig(0, out_edges, verteces)

print('After assembly, the contiguous sequence is:')
print(f'{contig}')

print('Original length: %s    After assembly length: %s' % (len(full_seq), len(contig)))
