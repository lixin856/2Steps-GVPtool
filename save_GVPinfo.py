# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 16:41:24 2021

@author: lixin
"""


import re


def KR_digestion(seq):
    pepls = []
    KR = sorted([i.start() for i in re.finditer('K', seq)] + [i.start() for i in re.finditer('R', seq)])
    pre = 0
    for i in KR:
        pep = seq[pre:i+1]
        if pep not in pepls:
            pepls.append(pep)
        pre = i+1
    pep = seq[pre:]
    if pep not in pepls:
        pepls.append(pep)
    return pepls


def get_dif_peptide(rawSeq, snpSeq):
    rawpepls = KR_digestion(rawSeq)
    snppepls = KR_digestion(snpSeq)
    
    same = [x for x in rawpepls if x in snppepls]
    dif_rawpepls = [y for y in rawpepls if y not in same]
    dif_snppepls = [z for z in snppepls if z not in same]
    return dif_rawpepls, dif_snppepls


def get_table_line(infoline, dif_rawpepls, dif_snppepls):
    global peplenCutoff
    tableline = infoline.split(' ')[0][1:] #line Num
    tableline = tableline + '\t' + infoline.split(' ')[1] # NM num
    tableline = tableline + '\t' + infoline.split(' ')[2] # SNP in nucleic acid
    tableline = tableline + '\t' + infoline.split(' ')[3] # SNP in protein
    
    # eg: line5938665, too long to write to csv file
    changeInfo = ' '.join(infoline.split(' - ')[0].split(' ')[4:])
    if len(changeInfo) > 500:
        tableline = tableline + '\t' + changeInfo.split(' ')[0] + ' ...'
    else:
        tableline = tableline + '\t' + changeInfo
    
    geneinfols = infoline.split(' - ')[-1].split('\t')[1].split(',')
    for geneinfo in geneinfols:
        if infoline.split(' ')[1] in geneinfo:
            gene = geneinfo.split(':')[0]
    tableline = tableline + '\t' + gene + '\t' + infoline.split(' - ')[-1].replace(',',';').replace('\n','')
    
    # eg: line5938665, too long to write to csv file
    if len(dif_rawpepls) > 10 and len(dif_snppepls) > 10:
        tableline = tableline + '\t' + '/'.join(dif_rawpepls[:10]) + '/...\t' + '/'.join(dif_snppepls[:10]) + '/...'
    elif len(dif_rawpepls) > 10:
        tableline = tableline + '\t' + '/'.join(dif_rawpepls[:10]) + '/...\t' + '/'.join(dif_snppepls)
    elif len(dif_snppepls) > 10:
        tableline = tableline + '\t' + '/'.join(dif_rawpepls) + '\t' + '/'.join(dif_snppepls[:10]) + '/...'
    else:
        tableline = tableline + '\t' + '/'.join(dif_rawpepls) + '\t' + '/'.join(dif_snppepls)
    
    difIL_rawpep = '/'.join(dif_rawpepls).replace('I','L')
    difIL_snppep = '/'.join(dif_snppepls).replace('I','L')
    if len(dif_snppepls) >= 1:
        if len(dif_snppepls) == 1 and len(dif_snppepls[0]) >= peplenCutoff:
            if difIL_snppep not in difIL_rawpep:
                snptype = 'yes' # useful, only 1 peptide, and replace IL
            else:
                snptype = 'yes_but_ILsame' # only 1 peptide, but after replace IL, snp and raw are the same
        else:
            snptype = 'not1pep_OR_length-min' # not only 1 peptide, or peptide length < cutoff
    else:
        snptype = '0snp-pep' # 0 snp peptide
    tableline = tableline + '\t' + snptype + '\n'
    
    return tableline


peplenCutoff = 4
rawfastafn = 'c:/data/hair/paper/dbSNP/common_all_20180418_hg38_gene.fasta'
tablefn = 'c:/data/hair/paper/commonSNP_all/commonSNP_all_info.csv'
fn = open(tablefn, 'w')
with open(rawfastafn, 'r') as file:
    for line in file:
        if line.startswith('>'):
            infoline = line            
        else:
            rawSeq = line.split('\t')[0]
            snpSeq = line.replace('\n','').split('\t')[1]
            dif_rawpepls, dif_snppepls = get_dif_peptide(rawSeq, snpSeq)
            tableline = get_table_line(infoline, dif_rawpepls, dif_snppepls)
            fn.write(tableline.replace('\t',','))
fn.close()

