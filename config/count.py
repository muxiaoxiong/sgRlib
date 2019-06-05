#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-04-22 17:26:11
# @Author  : muxiaoxiong
# @email   : niexiongwei@qq.com

# count.py
import sys
sys.path.append("..")
from config.seqmap import seqmap
from config.sgrna import sgrna
import os
from itertools import islice


class count:
    def __init__(self, inputfile, genome, N):
        self.inputfile = inputfile
        self.genome = genome
        self.N = N

    def select_N(self):
        sgrnafile = seqmap(self.inputfile, self.genome).seqmap_5M()
        N_dict = {}
        outname = "sgrna_Top_N.txt"
        if os.path.exists(outname):
            os.remove(outname)
        out = open(outname, "a+")
        with open(sgrnafile) as ff:
            for line in ff:
                line = line.strip()
                test = line.split()
                gene = test[0].split("@")[0]
                if gene not in N_dict.keys():
                    N_dict[gene] = 1
                else:
                    N_dict[gene] += 1
                if N_dict[gene] <= int(self.N):
                    out.write(line + "\n")
        out.close()
        return outname

    def count_total(self):
        sgrnafile = "sgRNA_NGG.txt"
        gene_2_sgrna = {}
        with open(sgrnafile) as ff:
            for line in ff:
                line = line.strip()
                if line.startswith(">"):
                    sgrna = line.replace(">", "")
                    gene = sgrna.split("@")[0]
                    if gene not in gene_2_sgrna.keys():
                        gene_2_sgrna[gene] = [0, 0, 0]
                        gene_2_sgrna[gene][0] += 1
                    else:
                        gene_2_sgrna[gene][0] += 1
        return gene_2_sgrna

    def count_candidata(self):
        gene_2_sgrna = self.count_total()
        with open("sgrna_candidata.txt") as ff:
            for line in ff:
                line = line.strip()
                if line.startswith(">"):
                    sgrna = line.replace(">", "")
                    gene = sgrna.split("@")[0]
                    if gene not in gene_2_sgrna.keys():
                        gene_2_sgrna[gene] = [0, 0, 0]
                        gene_2_sgrna[gene][1] += 1
                    else:
                        gene_2_sgrna[gene][1] += 1
        return gene_2_sgrna

    def count_N(self):
        sgrna_N = self.select_N()
        gene_2_sgrna = self.count_candidata()
        with open(sgrna_N) as ff:
            for line in islice(ff, 1, None):
                line = line.strip()
                test = line.split()
                gene = test[0].split("@")[0]
                gene_2_sgrna[gene][2] += 1
        outname = "sgrna_info.txt"
        if os.path.exists(outname):
            os.remove(outname)
        out = open(outname, "a+")
        out.write("Gene\tTotal sgRNA\tCandidata sgRNA\tTOP N sgRNA\n")
        for gene, info in gene_2_sgrna.items():
            out.write(gene + "\t")
            for i in info:
                out.write(str(i) + "\t")
            out.write("\n")
        out.close()
