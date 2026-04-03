# _*_ encding: utf-8 _*_
'''
@File : Format_QuIBL_final.py
@TIME : 2022/10/20 10::51
@Original Author: Chengran Zhou & Fang Li 
@Modified by: Jinjiazheng & Jiangchuan
@VERSION : 1.0
@Contact : jinjiazhengxiao@163.com
@LISCENCE : None
'''

#### This script is used to deal quilbal outfile
import sys,os,re  
import math,glob 
import ete3
import subprocess

def generate_all_triplet(files):
    file = files[0]
    print(file)
    with open(file,'r') as f:
        f.readline()
        triplet_lst = [line.split(',')[0] for line in f]
        triplet_lst = [i.strip().split('_') for i in set(triplet_lst)]
        return triplet_lst

def checkT1():
    toplis = [(sp1,sp2,sp3),(sp2,sp1,sp3),(sp3,sp1,sp2)]
    merdict = {}
    for top in toplis:
        lca = "(%s,(%s,%s));"%(top)
        ta = ete3.Tree(lca)
        if Atopology.compare(ta)['rf'] == 0:
            merdict["T1"] = top
        else:
            if "T2" in merdict:
                merdict["T3"] = top
            else:
                merdict["T2"] = top
    return merdict

def ReadOut(file):
    Mergedic = {}
    total_tree = 0
    with open(file,'r') as h1:
        h1.readline()
        for line in h1:
            iterms = line.strip().split(',')
            spl = iterms[0].strip().split('_')
            if len(set(spl+[sp1,sp2,sp3])) !=3:
                continue
            total_tree += float(iterms[10])
            if iterms[1] == outgroup["T1"][0]:
                Toplog = "T1"
                or1,or2,or3 = outgroup["T1"]
            elif iterms[1] == outgroup["T2"][0]:
                Toplog = "T2"
                or1,or2,or3 = outgroup["T2"]
            else:
                Toplog = "T3"
                or1,or2,or3 = outgroup["T3"]
            struc = "(%s,(%s,%s))"%(or1,or2,or3)
            if float(iterms[8]) -float(iterms[9]) < -10:
                Mergedic[Toplog] = {'Introgression':float(iterms[5])*float(iterms[10]),'ILS':float(iterms[4])*float(iterms[10]),"C2":float(iterms[3]),"Topo":struc}
            else:
                Mergedic[Toplog] = {'Introgression':0,'ILS':float(iterms[10]),"C2":float(iterms[3]),"Topo":struc}
    for k in ["T1","T2","T3"]:
        for k1 in ['Introgression','ILS']:
            Mergedic[k][k1] = Mergedic[k][k1]/float(total_tree)
    return Mergedic

def SubDeal(inputls):
        line1,line2 = inputls
        li1ls = line1.strip().split()
        li2ls = line2.strip().split()
        mark = "raw"
        ILSv = float(li1ls[4])
        INTv = float(li2ls[4])
        TopILS = li1ls[-1]
        TopINT = li2ls[-1]
        ILSvk = 0
        INTvk = 0
        if ILSv > 0.005:
                ILSvk = ILSv
                if INTv > 0.005:
                        mark = "%s%sILS%sINTR"%(li1ls[-1],ILSv,INTv)
                        INTvk = INTv
                else:
                        mark = "%s%sILS"%(li1ls[-1],ILSv)
        else:
                if (INTv > 0.005):
                        mark = "%s%sINTR"%(li2ls[-1],INTv)
                        INTvk = INTv
        return mark,ILSvk,INTvk,re.findall(r'\w+',li1ls[-1])

def addict(mark,ILSvkT,INTvkT,out,sp):
        if mark == "raw":
                return 0
        if Dichas[out][sp] == "NONE":
                Dichas[out][sp] = mark
        else:
                Dichas[out][sp] += ",%s"%(mark)
        if Dichas[sp][out] == "NONE":
                Dichas[sp][out] = mark
        else:
                Dichas[sp][out] += ",%s"%(mark)
        diclis[out][sp]['num'] += 1
        diclis[out][sp]['val'] += ILSvkT
        dicintr[out][sp]['num'] += 1
        dicintr[out][sp]['val'] += INTvkT
        diclis[sp][out]['num'] += 1
        diclis[sp][out]['val'] += ILSvkT
        dicintr[sp][out]['num'] += 1
        dicintr[sp][out]['val'] += INTvkT

def Readout_2(files):
    dics = {}
    ref = None
    with open(files,'r') as h1:
        h1.readline()
        for line1 in h1:
            iterms = line1.strip().split()
            key = iterms[0]
            if iterms[1] == "REF":
                ref = iterms[0]
            valls = iterms[2].split("|")
            dics[key] = [round(float(valls[0])/100,2),round(float(valls[1])/100,2)]
    return dics,ref

def main():
    global sp1,sp2,sp3,outgroup,Atopology
    rawpath,treefile = sys.argv[:2]
    rootspciesls = sys.argv[2:]
    OUTfileLs = glob.glob(os.path.join(rawpath,'Out.csv*'))
    all_triplet = generate_all_triplet(OUTfileLs)
    os.mkdir('./out1')
    for spn_lst in all_triplet:
        sp1,sp2,sp3 = spn_lst
        with open(treefile,'r') as h1:
            Trestr = h1.readline().strip()
            Atopology = ete3.Tree(Trestr)
            node = Atopology.get_common_ancestor(rootspciesls[0],rootspciesls[1]) if len(rootspciesls) >1 else rootspciesls[0]
            Atopology.set_outgroup(node)
            Atopology.prune([sp1,sp2,sp3])
        outgroup = checkT1()
        outfile = open("./out1/%s_%s_%s.xls"%(sp1,sp2,sp3),'w')
        outfile.write('Triplet\tIndex\tTopology\tType\tratio\ttimes\tTopo_struct\n')
        for Outfil in OUTfileLs:
            idx = Outfil.strip().split('.')[-1]
            indx = "Random_%s"%(idx)
            oudict = ReadOut(Outfil)
            for Topy in ["T1","T2","T3"]:
                Thedic = oudict[Topy]
                writeL1 = ["_".join([sp1,sp2,sp3]),indx,Topy,"ILS",str(Thedic['ILS']),str(Thedic['C2']),Thedic['Topo']]
                writeL2 = ["_".join([sp1,sp2,sp3]),indx,Topy,"Introgression",str(Thedic['Introgression']),str(Thedic['C2']),Thedic['Topo']]
                outfile.write('%s\n%s\n'%('\t'.join(writeL1),'\t'.join(writeL2)))
        outfile.close()
# Merge and average
    merge_dic,Number_dic,Time_dic,treety_dic = {},{},{},{}
    for file in os.listdir('./out1'):
        file = os.path.join('./out1',file)
        with open(file,'r') as h1:
            h1.readline()
            for line1 in h1:
                iterms = line1.strip().split()
                Key,index,Topology,Type,ratio,times,Topo_struct = iterms
                if Key not in treety_dic:
                    treety_dic[Key] = {}
                if Topology not in treety_dic[Key]:
                    treety_dic[Key][Topology] = Topo_struct
                try:
                    merge_dic[Key][Topology][Type] += float(ratio)
                    Number_dic[Key][Topology][Type] += 1
                    Time_dic[Key][Topology][Type] += float(times)
                except:
                    if Key not in merge_dic:
                        merge_dic[Key] = {}
                        Number_dic[Key] = {}
                        Time_dic[Key] = {}
                    if Topology not in merge_dic[Key]:
                        merge_dic[Key][Topology] = {}
                        Number_dic[Key][Topology] = {}
                        Time_dic[Key][Topology] = {}
                    merge_dic[Key][Topology][Type] = float(ratio)
                    Number_dic[Key][Topology][Type] = 1
                    Time_dic[Key][Topology][Type] = float(times)
    outf = open('./out1/allstat.xls','w')
    outf.write("ID\tTriplet\tTopology\tType\tAverage\tTreesNumber\ttimes\ttree\n")
    n = 0
    for KeyN,Kdic in merge_dic.items():
        for topology in ['T1','T2','T3']:
            trees = treety_dic[KeyN][topology]
            for ty in ["ILS",'Introgression']:
                n += 1
                aver2 = round(Kdic[topology][ty]/Number_dic[Key][Topology][Type],6)
                time2 = round(Time_dic[KeyN][topology][ty]/Number_dic[Key][Topology][Type],6)
                num = aver2*0
                outf.write('%s\n'%('\t'.join([str(n),KeyN,topology,ty,str(aver2),str(num),str(time2),trees])))
    outf.close()

# Calculate each pairwise species
    global Dichas,dictchek,diclis,dicintr
    Spnls = list(set(sum(all_triplet,[])))
    Dichas,dictchek,diclis,dicintr = {},{},{},{}
    num = len(Spnls)
    for i in range(num):
        spn = Spnls[i]
        Dichas[spn] = {}
        dictchek[spn] = {}
        diclis[spn] = {}
        dicintr[spn] = {} 
        for i1 in range(num):
            spn1 = Spnls[i1]
            dictchek[spn][spn1] = 0
            diclis[spn][spn1] = {'val':0,'num':0}
            dicintr[spn][spn1] = {'val':0,'num':0}
            if spn == spn1:
                Dichas[spn][spn1] = "OUT"
            else:
                Dichas[spn][spn1] = "NONE"
    with open('./out1/allstat.xls','r') as f:
        f.readline()
        while True:
            ls = [[f.readline() for i in range(2)]for i1 in range(3)]
            if not(ls[0][0]):
                break
            T1ls,T2ls,T3ls = ls
            out,sp1,sp2 = re.findall(r'\w+',T1ls[0].strip().split()[-1])
            markT2,ILSvkT2,INTvkT2,T2spls = SubDeal(T2ls)
            markT3,ILSvkT3,INTvkT3,T3spls = SubDeal(T3ls)
            SP0 = out
            if sp1 == T2spls[0]:
                SP1 = sp2
                SP2 = sp1
            else:
                SP1 = sp1
                SP2 = sp2
            dictchek[out][SP1] += 0
            dictchek[out][SP2] += 0
            addict(markT2,ILSvkT2,INTvkT2,out,SP1)
            addict(markT3,ILSvkT3,INTvkT3,out,SP2)
    os.mkdir('./out2')
    for i1 in range(num):
        spn1 = Spnls[i1]
        outf = open('./out2/%s.txt'%(spn1),'w')
        outf.write("Id\tType\tILSratio|INTRratio\tILSnum|INTRnum\tToponumber\tnote\n")
        for i2 in range(num):
            spn2 = Spnls[i2]
            Mars = Dichas[spn1][spn2]
            rils = 0
            rintr = 0
            if Mars == "NONE":
                types = "NONE"
            elif Mars == "OUT":
                types = "REF"
            elif "ILS" in Mars and "INTR"  in Mars:
                types = "BOTH"
            elif "ILS" in Mars:
                types = "ILS"
            else:
                types = "INTR"
            if types in ["ILS",'BOTH']:
                rils = round(float(diclis[spn1][spn2]['val']*100)/diclis[spn1][spn2]['num'],2)
                rintr = round(float(dicintr[spn1][spn2]['val']*100)/dicintr[spn1][spn2]['num'],2)
            wirtels = [spn2,types,"%s|%s"%(rils,rintr),"%s|%s"%(diclis[spn1][spn2]['num'],dicintr[spn1][spn2]['num']),str(dictchek[spn1][spn2]),str(Dichas[spn1][spn2])]
            outf.write('%s\n'%('\t'.join(wirtels)))
        outf.close()

#Generate matrix
    outfile = open('ILS_introng.txt','w')
    outfile.write('speciesA\tspeciesB\tratio\tType\n')
    Speciesoder = {}; n=0
    for spn in Spnls:
        n += 1
        Speciesoder[spn] = n
    MerDict = {}
    for file in os.listdir('./out2'):
        file = os.path.join('./out2',file)
        dim,ref = Readout_2(file)
        if ref == None:
            continue
        for k,v in dim.items():
            k1 = "%s_%s"%(ref,k)
            k2 = "%s_%s"%(k,ref)
            if k1 in MerDict or k2 in MerDict:
                continue
            MerDict[k1] = [str(i) for i in v]
    for keys,vlis in MerDict.items():
        sp1,sp2 = keys.split("_")
        if Speciesoder[sp1] >Speciesoder[sp2]:
            ILSty = '%s\t%s'%(sp1,sp2)
            INTty = '%s\t%s'%(sp2,sp1)
        else:
            ILSty = '%s\t%s'%(sp2,sp1)
            INTty = '%s\t%s'%(sp1,sp2)
        ILSty = '%s\t%s'%(sp1,sp2) if Speciesoder[sp1] >Speciesoder[sp2] else '%s\t%s'%(sp2,sp1)
        INTty = '%s\t%s'%(sp1,sp2) if Speciesoder[sp1] <Speciesoder[sp2] else '%s\t%s'%(sp2,sp1)
        if sp1 == sp2:
            outfile.write('%s\t0.0\tBOTH\n'%(ILSty))
            continue
        outfile.write('%s\t%s\tILS\n%s\t%s\tIntrongression\n'%(ILSty,vlis[0],INTty,vlis[1]))
    outfile.close()    

if __name__ == "__main__":
    if len(sys.argv) <3:
        sys.exit("python script.py QuIBL_out_path Species_tree || PS:The output files in the raw_path directory should be in the format Out.csv.1, Out.csv.2, ...")
    main()
