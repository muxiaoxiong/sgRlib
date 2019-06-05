#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-04-19 13:52:05
# @Author  : muxiaoxiong
# @email   : niexiongwei@qq.com

import re
import os
import sys
sys.path.append("..")
from config.gene import gene


class sgrna():
    def __init__(self, inputfile, pam="GG", lenth=21):
        self.id_seq = gene(inputfile).gene_2_dict()
        self.pam = pam
        self.lenth = lenth

    def gene_2_sgRNA(self):
        # 清楚TTTT，清除重复
        id_sgrna = {}
        for genid, seq in self.id_seq.items():
            count = 0
            pam_list = find_sgrna(seq, self.pam)
            for i in range(len(pam_list)):  # 正义链
                end = pam_list[i] - 1
                start = pam_list[i] - self.lenth
                if start >= 0 and end <= len(seq):
                    count += 1
                    sequence = seq[start:end]  # 设计的sgRNA序列（正义链）
                    sgrna_name = genid + "@S_" + str(count)  # sgRNA的id
                    if sequence in id_sgrna.keys():
                        id_sgrna[sequence].append(sgrna_name)
                    else:
                        id_sgrna[sequence] = []
                        id_sgrna[sequence].append(sgrna_name)
            count = 0
            rever_seq = fasta_reverse(seq)  # 反向互补序列
            repam_list = find_sgrna(rever_seq, self.pam)
            for i in range(len(repam_list)):  # 正义链
                end = repam_list[i] - 1
                start = repam_list[i] - self.lenth
                if start >= 0 and end <= len(rever_seq):
                    count += 1
                    sequence = rever_seq[start:end]  # 设计的sgRNA序列（反义链）
                    sgrna_name = genid + "@A_" + str(count)  # sgRNA的id
                    if sequence in id_sgrna.keys():
                        id_sgrna[sequence].append(sgrna_name)
                    else:
                        id_sgrna[sequence] = []
                        id_sgrna[sequence].append(sgrna_name)
        for sequence in list(id_sgrna.keys()):
            if len(id_sgrna[sequence]) != 1:
                gene_list = set()
                for i in id_sgrna[sequence]:
                    i = i.split("@")[0]
                    gene_list.add(i)
                if len(i) > 1:
                    del id_sgrna[sequence]
        id_sgrna = remove_id_sgr_TTTT(id_sgrna)  # 去除了序列中的TTTT
        id_sgrna = add_NGG(id_sgrna, self.pam)  # 为sgRNA加上PAM
        return id_sgrna

    def sgrna_2_file(self):
        id_sgrna = self.gene_2_sgRNA()
        id_sgrna_file_name = "sgRNA_NGG.txt"
        if os.path.exists(id_sgrna_file_name):
            os.remove(id_sgrna_file_name)
        out = open(id_sgrna_file_name, "a+")
        for sgrid, seq in id_sgrna.items():
            out.write(">" + sgrid + "\n")
            out.write(seq + "\n")
        out.close()
        return id_sgrna_file_name


def fasta_reverse(sequence):
    '''获得反向互补序列'''
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    sequence = sequence.upper()
    return sequence[::-1]


def remove_id_sgr_TTTT(id_sgrna):
    '''删除id_sgrna中活性过低的sgrna'''
    sgrna = list(id_sgrna.keys())
    for i in sgrna:
        if re.search('TTTT', i):
            del id_sgrna[i]
    return id_sgrna


def add_NGG(id_sgrna, pam="GG"):
    new_id_sgrna = {}
    for sequence, seq_list in id_sgrna.items():
        for i in seq_list:
            new_id_sgrna[i] = sequence + str("N") + pam
    return new_id_sgrna


def find_sgrna(seq, pam):
    test = set()
    for i in range(len(seq)):
        x = seq[i:]
        y = re.search(pam, x)
        try:
            test.add(y.start() + i)
        except:
            pass
    return list(test)
