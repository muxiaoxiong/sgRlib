#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Date    : 2019-04-19 13:44:15
# @Author  : muxiaoxiong
# @email   : niexiongwei@qq.com

import re
import os
import sys


class gene():
    def __init__(self, inputfile):
        self.inputfile = inputfile

    def gene_2_dict(self):
        id_seq = {}
        with open(self.inputfile)as fh:             # 文件打开的位置
            id, seq = '', ''
            for line in fh:
                if line.startswith('#'):
                    continue
                line = line.strip()
                if line.startswith('>'):
                    if id:
                        if id_seq.get(id):
                            if len(seq) > len(id_seq[id]):
                                id_seq[id] = seq
                                seq = ''  # 清除多余的字典
                        else:
                            id_seq[id] = seq
                            seq = ''
                    gene = re.findall('gene=[\S]+', line)  # 寻找里面是否有gene
                    name = line.replace('>', '').split()[0]  # 如果没有就以第一个字符串命名
                    if len(gene):
                        id = re.sub('[^a-zA-Z0-9_-]', '',
                                    gene[0]).replace("gene", "")
                    else:
                        id = name
                else:
                    seq += line.upper().replace("U", 'T')
            id_seq[id] = seq                            # 将序列做成字典键值对的值
        return id_seq

    def gene_2_cleanfile(self):
        outname = self.inputfile.split(".")[0] + "_clean.txt"
        if os.path.exists(outname):
            return outname
        else:
            id_seq = self.gene_2_dict()
            with open(outname, "w") as out:
                for geneid, seq in id_seq.items():
                    out.write(">" + geneid + "\n")
                    out.write(seq.upper() + "\n")
            return outname
