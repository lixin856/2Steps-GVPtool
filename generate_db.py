# -*- coding: utf-8 -*-
"""
Created on Tue Mar  2 14:45:03 2021

@author: lixin
"""


def read_refDB(fn):
    refDBseq = ''
    geneProdic = {}
    with open(fn, 'r') as file:
        for line in file:
            if line.startswith('>'):
                refDBseq = refDBseq + ' '
                gene1 = line.replace('\n',' ').split('GN=')[-1].split(' ')[0]
                gene2 = line.split('_')[0].split('|')[-1]
                protein = line.split('|')[1]
                geneProdic[gene1] = [gene2, protein]
            else:
                refDBseq = refDBseq + line.replace('\n', '').replace('I', 'L')
    return refDBseq, geneProdic


refDBfn = 'c:/data/hair/paper/hair_sp/SwissProt_iso_human.fasta'
refDBseq, geneProdic = read_refDB(refDBfn)

rawtablefn = 'C:/data/hair/paper/commonSNP_all/commonSNP_all_info.csv'
sameCount = 0
difCount = 0
genedic = {}
with open(rawtablefn, 'r') as file:
    for line in file:
        if line.replace('\n','').split(',')[-1] == 'yes':
            snppepIL = line.split(',')[-2].replace('I', 'L')
            if snppepIL not in refDBseq:
                difCount += 1
                gene = line.split(',')[5]
                snppep = line.split(',')[-2]
                infols = [snppepIL] + line.split(',')[1:4] + line.split(',')[7:11] + line.split(',')[-3:-2]
                if gene not in genedic.keys():
                    genedic[gene] = {}
                    genedic[gene][snppep] = [infols]
                else:
                    if snppep not in genedic[gene].keys():
                        genedic[gene][snppep] = [infols]
                    else:
                        genedic[gene][snppep].append(infols)
            else:
                sameCount += 1


newtablefn = 'C:/data/hair/paper/commonSNP_all/commonSNP_pep_info.csv'
with open(newtablefn, 'w') as file:
    for (gene, value1) in genedic.items():
        file.write(gene + ',')
        n = 0
        for (pep, infols) in value1.items():
            if n == 0:
                file.write(pep + ',')
                n += 1
            else:
                file.write(',' + pep + ',')
            
            m = 0
            for info in infols:
                if m == 0:
                    file.write(','.join(info) + '\n')
                    m += 1
                else:
                    file.write(',,' + ','.join(info) + '\n')


snpDBfn = 'C:/data/hair/paper/commonSNP_all/commonSNP_snppepDB.fasta'
with open(snpDBfn, 'w') as file:
    for (gene, value1) in genedic.items():
        pepls = list(genedic[gene].keys())
        endPepNumls = []
        krPepNumls = []
        for i in range(len(pepls)):
            if pepls[i][-1] != 'K' and pepls[i][-1] != 'R':
                endPepNumls.append(i)
            else:
                krPepNumls.append(i)
        
        if gene in geneProdic.keys():
            protein = geneProdic[gene][1]
            gene2 = geneProdic[gene][0]
        else:
            protein = gene
            gene2 = gene
            
        if len(endPepNumls) == 0:
            head = '>sp|' + protein + '.m0|' + gene2 + '.m0_HUMAN OS=Homo sapiens OX=9606 GN=' + gene + '.m0\n'
            file.write(head + ''.join(pepls) + '\n')
        else:
            if len(krPepNumls) != 0:
                seq = ''
                for num in krPepNumls:
                    seq = seq + pepls[num]
                head = '>sp|' + protein + '.m0|' + gene2 + '.m0_HUMAN OS=Homo sapiens OX=9606 GN=' + gene + '.m0\n'
                file.write(head + seq + pepls[endPepNumls[0]] + '\n')
            
            if len(endPepNumls) > 1:
                j = 1
                for num in endPepNumls[1:]:
                    head = '>sp|' + protein + '.m' + str(j) + '|' + gene2 + '.m' + str(j) + '_HUMAN OS=Homo sapiens OX=9606 GN=' + gene + '.m' + str(j)
                    file.write(head + '\n' + pepls[num] + '\n')
                    j += 1

