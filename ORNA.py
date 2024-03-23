from gatb import Bank, BankFasta, Kmer, Graph, System, Sequence, Node
from gatb.core import Iterator, ProgressIterator, PairedIterator, GraphIterator, Dispatcher
import math
import os

STR_BASE = "-base"
STR_INPUT = "-input"
STR_OUTPUT = "-output"
STR_KMER = "-kmer"
STR_PAIR1 = "-pair1"
STR_PAIR2 = "-pair2"
STR_SORTING = "-sorting"
STR_KSORTING = "-ksorting"
STR_BSIZE = "-binsize"
COLL_STAT = "-cs"
STR_TYPE = "-type"

def get_length(bank):
    it = bank.iterator()
    max_length = 0
    for sequence in it:
        max_length = sequence.getDataSize()
        break
    return max_length

def get_number(bank):
    it = bank.iterator()
    sequence_count = 0
    for _ in it:
        sequence_count += 1
    return sequence_count

def get_nodes_number(graph_iterator):
    count = 0
    for _ in graph_iterator:
        count += 1
    return count

def get_median(graph, seq, kmer, nodeit):
    read_abundance = 0
    i = 0
    length = seq.getDataSize()
    abdce = []
    model = Kmer(kmer).ModelCanonical
    itKmer = model.Iterator(model)
    itKmer.setData(seq.getData())
    data = seq.getDataBuffer()
    for _ in itKmer:
        s = str(model.toString(itKmer.value()))
        node = graph.buildNode(s)
        index = graph.nodeMPHFIndex(node)
        abund = nodeit[index]
        abdce.append(abund)
    abdce.sort()
    i = len(abdce)
    if i % 2 == 0:
        median = (abdce[i // 2] + abdce[i // 2 - 1]) // 2
    else:
        median = abdce[i // 2]
    return median

def calculate_quality(qual):
    score = 0
    for char in qual:
        score += ord(char)
    return score

def calculate_abundance(graph, seq, itKmer, model):
    score = 0
    for _ in itKmer:
        s = str(model.toString(itKmer.value()))
        node = graph.buildNode(s)
        threshold = 0
        index = graph.nodeMPHFIndex(node)
        abund = int(graph.queryAbundance(node))
    return score

def read_acceptance(graph, itKmer, model, counter, base, nodeit):
    acceptance = 0
    for _ in itKmer:
        s = str(model.toString(itKmer.value()))
        node = graph.buildNode(s)
        threshold = 0
        index = graph.nodeMPHFIndex(node)
        abund = nodeit[index]
        threshold = math.ceil((math.log(abund) / math.log(base)))
        if threshold > abund:
            threshold = abund
        if threshold < 1:
            threshold = 1
        if counter[index] < threshold:
            acceptance += 1
            break
    return acceptance

def single_end(graph, filename, out_file, base, kmer, nb_cores, nodeit, seq_type):
    count = 0
    dispatcher = Dispatcher(nb_cores)
    sync = System.thread().newSynchronizer()
    it = graph.iterator()
    input_dataset = Bank.open(filename)
    it_seq = ProgressIterator(input_dataset)
    out_dataset = None
    tmp_file = out_file
    out_file1 = ""
    extension = ""
    final_out = ""

    if seq_type.lower() in ["fastq", "fq"]:
        extension = ".fq"
    else:
        extension = ".fa"

    out_file1 = tmp_file + extension
    final_out = out_file1

    if seq_type.lower() in ["fastq", "fq"]:
        out_dataset = BankFasta(final_out, True)
    else:
        out_dataset = BankFasta(final_out)

    node_size = it.size()
    counter = [0] * node_size

    for sequence in it_seq:
        length = sequence.getDataSize()
        flag = 1
        gb = 1
        acceptance = 0
        model = Kmer(kmer).ModelCanonical
        it_kmer = model.Iterator(model)
        it_kmer.setData(sequence.getData())

        if flag == 1:
            acceptance = read_acceptance(graph, it_kmer, model, counter, base, nodeit)

        if acceptance > 0:
            for _ in it_kmer:
                s = str(model.toString(it_kmer.value()))
