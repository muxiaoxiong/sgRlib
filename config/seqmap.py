#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-04-22 09:50:34
# @Author  : muxiaoxiong
# @email   : niexiongwei@qq.com

# seqmap.py
import os
import sys
sys.path.append("..")
import platform
from itertools import islice
from collections import Counter
from config.gene import gene
from config.sgrna import sgrna


class seqmap():
    '''用来对生成的sgRNA文件进行比对'''

    def __init__(self, inputfile, genome="", pam="GG", lenth=21):
        # clean_file文件为inputfile文件中的clean文件 用来做1M文件的比对
        self.clean_file = gene(inputfile).gene_2_cleanfile()
        self.genome = genome
        self.pam = pam
        self.lenth = lenth
        self.inputfile = inputfile

    def seqmap_1M(self):
        # 用来将增加NGG的文件与clean_file文件进行不比对
        seqmapfile = "NGG_vs_clean.txt"
        NGG_sgrna = sgrna(self.inputfile).sgrna_2_file()
        if os.path.exists(seqmapfile):
            os.remove(seqmapfile)
        choose_seqmap(1, NGG_sgrna, self.clean_file, seqmapfile)
        os.remove(self.clean_file)  # 删除clean_file
        try:
            open(seqmapfile)
        except:
            print("###")
            print("The software is limited by the Windows operating system. Please use other operating system to run the file.")
            sys.exit(2)
        sgrna_gene = {}
        no_repeat_sgrna = {}
        with open(seqmapfile) as ff:
            for line in islice(ff, 1, None):
                line = line.strip()
                test = line.split()
                sgrid = test[-4]
                gene = test[0]
                if sgrid not in sgrna_gene.keys():
                    sgrna_gene[sgrid] = set()
                    sgrna_gene[sgrid].add(gene)
                else:
                    sgrna_gene[sgrid].add(gene)
        for sgrid in list(sgrna_gene.keys()):
            if len(sgrna_gene[sgrid]) != 1:
                del sgrna_gene[sgrid]

        sgrna_output = '1M_sgrna.txt'
        out = open(sgrna_output, "a+")
        with open(seqmapfile) as ff:
            for line in islice(ff, 1, None):
                line = line.strip()
                test = line.split()
                sgr = test[3]
                sequence = test[2]
                if sgr in sgrna_gene.keys():
                    no_repeat_sgrna[sgr] = sequence
                    out.write(">" + sgr + "\n")
                    out.write(sequence + "\n")
        out.close()
        return sgrna_output, no_repeat_sgrna

    def seqmap_0M(self):
        seqmapfile = "sgrna_vs_genome.txt"
        if os.path.exists(seqmapfile):
            os.remove(seqmapfile)
        (sgrna_output, no_repeat_sgrna) = self.seqmap_1M()
        count_list = []
        choose_seqmap(0, sgrna_output, self.genome, seqmapfile)
        with open(seqmapfile) as ff:
            for line in islice(ff, 1, None):
                line = line.strip()
                test = line.split()
                sgrna_seq = test[-3]
                count_list.append(sgrna_seq)
        count_list = Counter(count_list)
        unqiue = [k for k, v in count_list.items()if v == 1]  # 选择0M=1的sgRNA
        candidata_name = "sgrna_candidata.txt"
        if os.path.exists(candidata_name):
            os.remove(candidata_name)
        with open(candidata_name, "a+") as out:
            for sgrid, sequence in no_repeat_sgrna.items():
                if sequence in unqiue:
                    out.write(">" + sgrid + "\n" + sequence + "\n")
        return candidata_name

    def seqmap_5M(self):
      # 用来与基因组进行比对 通常选择5M进行通过脱靶评估
        seqmapfile = "candidata_vs_genome.txt"
        sgrna_dict = {}
        l_pam = -int(len(self.pam))
        num = 5  # 脱靶评估的数目
        if os.path.exists(seqmapfile):
            os.remove(seqmapfile)
        candidata = self.seqmap_0M()
        choose_seqmap(5, candidata, self.genome, seqmapfile)
        with open(seqmapfile) as ff:
            for line in islice(ff, 1, None):
                line = line.strip()
                test = line.split()
                sgrna = test[-5]
                if sgrna[l_pam:] != self.pam:
                    pass
                else:
                    sgrnaid = test[-4]
                    if sgrnaid not in sgrna_dict.keys():
                        sgrna_dict[sgrnaid] = []
                        sgrna_dict[sgrnaid].append(test[-3])
                        sgrna_dict[sgrnaid].append(test[-2])
                    else:
                        sgrna_dict[sgrnaid].append(test[-2])
        outname = 'sgrna_off_target.txt'
        out = open(outname, "a+")
        out.write("随机sgRNAID\t序列\t")
        for i in range(num + 1):
            out.write(str(i) + "M\t")
        out.write("Total\n")
        for sgid, sginfo in sgrna_dict.items():
            totalcount = 0
            out.write(sgid + "\t")
            out.write(sginfo[0] + "\t")
            sgrnacount = Counter(sginfo[1:])
            for i in range(num + 1):
                i = str(i)
                if i not in sgrnacount.keys():
                    sgrnacount[i] = 0
                sgrna_sort = sorted(sgrnacount)
            for sort_sg in sgrna_sort:
                sgcount = int(sgrnacount[sort_sg])
                totalcount += sgcount
                out.write(str(sgrnacount[sort_sg]) + "\t")
            out.write(str(totalcount) + "\n")
        out.close()
        return outname


def win_seqmap(num, inputfile, genome, seqmapfile):
    path = os.getcwd()
    os.chdir(path + "\\config")
    inputfile = path + "\\" + inputfile
    genome = path + "\\" + genome
    seqmapfile = path + "\\" + seqmapfile
    os.system(("seqmap.exe {0} {1} {2} {3} /output_all_matches").format(num,
                                                                        inputfile, genome, seqmapfile))
    os.chdir(path)


def linux_seqmap(num, inputfile, genome, seqmapfile):
    path = os.getcwd()
    os.chdir(path + "/config")
    os.system("chmod +x ./seqmap")
    inputfile = path + "/" + inputfile
    genome = path + "/" + genome
    seqmapfile = path + "/" + seqmapfile
    os.system(("./seqmap {0} {1} {2} {3} /output_all_matches").format(num,
                                                                      inputfile, genome, seqmapfile))
    os.chdir(path)


def choose_seqmap(num, inputfile, genome, seqmapfile):
        # 用来对seqmap外部程序调用函数
    sysstr = platform.system()
    if(sysstr == "Windows"):
        win_seqmap(num, inputfile, genome, seqmapfile)
    elif(sysstr == "Linux"):
        linux_seqmap(num, inputfile, genome, seqmapfile)
    else:
        print("当前系统暂时不支持")
        sys.exit(2)
