# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 13:12:47 2020

@author: lixin
"""


import tkinter as tk
from tkinter import filedialog


def readResultPep(fn):
    global pepdic
    with open(fn, 'r') as file:
        linels = file.readlines()[1:]
        for line in linels:
            pepIL = line.split(',')[0]
            if pepIL not in pepdic.keys():
                pepdic[pepIL] = {}
                pepdic[pepIL]['pep'] = line.split(',')[1]
                pepdic[pepIL]['HQ-ORF'] = line.split(',')[2]
                pepdic[pepIL]['score'] = line.split(',')[3]
                pepdic[pepIL]['fdr'] = line.split(',')[4]
                pepdic[pepIL]['minFdr'] = line.split(',')[5]
                pepdic[pepIL]['i114'] = line.split(',')[6]
                pepdic[pepIL]['i115'] = line.split(',')[7]
                pepdic[pepIL]['i116'] = line.split(',')[8]
                pepdic[pepIL]['i117'] = line.replace('\n','').split(',')[9]
                pepdic[pepIL]['fnls'] = [fn]
                if eval(line.split(',')[5]) < 0.01:
                    pepdic[pepIL]['dcount'] = 1
                else:
                    pepdic[pepIL]['dcount'] = 0
            else:
                pepdic[pepIL]['fnls'].append(fn)
                pepdic[pepIL]['score'] = pepdic[pepIL]['score'] + ';' + line.split(',')[3]
                pepdic[pepIL]['fdr'] = pepdic[pepIL]['fdr'] + ';' + line.split(',')[4]
                if eval(line.split(',')[5]) < eval(pepdic[pepIL]['minFdr']):
                    pepdic[pepIL]['minFdr'] = line.split(',')[5]
                if eval(line.split(',')[5]) < 0.01:
                    pepdic[pepIL]['dcount'] += 1
                pepdic[pepIL]['i114'] = str(eval(pepdic[pepIL]['i114']) + eval(line.split(',')[6]))
                pepdic[pepIL]['i115'] = str(eval(pepdic[pepIL]['i115']) + eval(line.split(',')[7]))
                pepdic[pepIL]['i116'] = str(eval(pepdic[pepIL]['i116']) + eval(line.split(',')[8]))
                pepdic[pepIL]['i117'] = str(eval(pepdic[pepIL]['i117']) + eval(line.replace('\n','').split(',')[9]))


def readTargetGenedic(fn):
    genedic = {}
    with open(fn, 'r') as file:
        for line in file:
            target = line.split('\t')[3]
            if target != '':
                gene = line.split('\t')[1]
                if gene not in genedic.keys():
                    genedic[gene] = target
    return genedic


genedic = readTargetGenedic('C:/data/linzhen/proteogenomics/pacbio/HQ_all.txt')
pepdic = {}
root = tk.Tk()
root.withdraw()
fns = filedialog.askopenfilenames()
headline = 'pepIL,pep,HQ-ORF,gene,score,fdr,minFdr,dcount,count,i114,i115,i116,i117'
for fn in fns:
    readResultPep(fn)
    headline = headline + ',' + fn.split('/')[-1].split('_')[0]
with open('C:/data/linzhen/proteogenomics/pacbio/combine/sp/result_ORFpep-sp.csv', 'w') as file:
    file.write(headline + '\n')
    for (pep, value) in pepdic.items():
        line = pep + ',' + value['pep'] + ',' + value['HQ-ORF'] + ','
        
        HQorfls = value['HQ-ORF'].split(';')
        genels = []
        for HQorf in HQorfls:
            gene = HQorf.split('.')[0]
            if gene in genedic.keys():
                if genedic[gene] not in genels:
                    genels.append(genedic[gene])
        
        line = line + ';'.join(genels) + ',' + value['score'] + ',' + value['fdr']
        line = line + ',' + value['minFdr'] + ',' + str(value['dcount']) + ',' + str(len(value['fnls']))
        line = line + ',' + value['i114'] + ',' + value['i115'] + ',' + value['i116'] + ',' + value['i117']
        
        for fn in fns:
            if fn in value['fnls']:
                line = line + ',Y'
            else:
                line = line + ','
        file.write(line + '\n')

