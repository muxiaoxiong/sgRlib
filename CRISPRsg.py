#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-04-19 13:26:50
# @Author  : muxiaoxiong
# @email   : niexiongwei@qq.com
# @file    : CRISPRsg.py

'''CRISPERsg v_1.2
对v1.0版本中的问题进行相应的优化
程序调用几个模块
gene 1.将文件读取并且制作成dict文件
     2.将输入的文件进行clean再输出

sgRNA 1.接受gene_dict 并根据此设计sgRNA
      2.设计的sgRNA去掉重复及TTTT
      3.设计的sgRNA增加NGG
      4.输出NGG文件

seqmap 1.输出的NGG文件与 输入clean文件进行比较 1M 选择只靶向一个基因的sgRNA
     2.选择的sgRNA与全基因组基因比较 0M=1的sgRNA
     3.0M=1的设计RNA与全基因组比较 5M 并按照脱靶的数目选择相应的sgRNA 并进行脱靶评估的统计

count 1.统计每个基因设计的总sgRNA数目
      2.统计每个基因candidata_sgRNA数目
      3.统计5M的脱靶数目
      4.统计选择前N个 并输出GC含量 起始位置 等

'''
import os
import getopt
import sys
sys.path.append("..")
from config.count import count
from config.sgrna import sgrna


def usage():  # 使用信息
    print(
        """
Usage:
  python CRISPRsg.py -i <genefile> -g <reference_genome> -n <num sgrna of 1 gene>

Authors:
    Mu X.X
Email:
    xiongweinie@qq.com
Copyright:
    Free for research and educational users.
Support platform:
	Windows & Linux

In order to ensure a good experience of the program
This version is only suitable for reference genomes less than 100M.
If your reference genome is very large, please use the sgRcas9_lib
"""
    )


def main():
    inputfile = ''
    all_genome = ''
    num = 0
    if len(sys.argv) == 1:  # 无参数时 直接显示使用信息
        usage()
        sys.exit()
    else:
        opts, args = getopt.getopt(sys.argv[1:], "hvi:g:n:", [
                                   "help"])  # 定义输入的- &--的情况
        for opt, value in opts:
            if opt in ("-h", "--help"):  # 显示帮助信息
                usage()
                sys.exit()
            elif opt == '-i':
                inputfile = value  # 输入的文件信息
                try:
                    open(inputfile)
                except:
                    print("文件损坏或者文件名不正确，请参考使用手册\n")
                    sys.exit(2)
            elif opt == '-g':
                all_genome = value
            elif opt == '-n':
                num = value
    old_file = ["sgrna_candidata.txt", "sgrna_Top_N.txt", "sgr_info.txt",
                "sgrna_off_target.txt", "sgrna_vs_genome", "sgRNA_NGG.txt", "NGG_vs_clean.txt"]

    for i in old_file:
        try:
            if os.path.exists(i):
                os.remove(i)
        except:
            pass

    count(inputfile, all_genome, num).count_N()
    os.remove("1M_sgrna.txt")
    os.remove("candidata_vs_genome.txt")
    os.remove("NGG_vs_clean.txt")
    os.remove("sgrna_vs_genome.txt")
    os.remove("sgRNA_NGG.txt")


if __name__ == '__main__':
    main()
