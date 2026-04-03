import sys,os
import pysam
import glob
def ReadPeP(files):
	h1 = open(files,'r')
	dics = {}
	for line in h1:
		if line.startswith('>'):
			ids = line.strip().split('>')[1]
			dics[ids] = ""
			continue
		if line.startswith('\n'):
			continue
		dics[ids] += line.strip()
	h1.close()
	return dics 

maffpaht = "/mafft_output"
signcopy = "/path/to/OrthoFinder_Output/Single_Copy_Orthologue_Sequences"
genome_path = "/path/to/species_genome"
pepfiledict = {}
cdsfiledict = {}
spnfile = open('/path/to/spn.txt','r')
for line in spnfile:
	spn = line.strip().split()[0]
	spnpepfile = "%s/%s/%s.pep.fa"%(genome_path,spn,spn)
	spncdsfile = "%s/%s/%s.cds.fa"%(genome_path,spn,spn)
	pepfiledict[spn] = pysam.Fastafile(spnpepfile)
	cdsfiledict[spn] = pysam.Fastafile(spncdsfile)
spnfile.close()

algnmentfilels  = glob.glob('%s/*_aligned.fa'%(maffpaht))
Checkresult = open('check.txt','w')
for aligfile in algnmentfilels:
	names = os.path.basename(aligfile).split('_aligned.fa')[0]
	rawpepfile = "%s/%s.fa"%(signcopy,names)
	rawdict = ReadPeP(rawpepfile)
	spnoder = []
	recdsseq = {}
	checkwrit = []
	for genid,seq in rawdict.items():
		spns = genid.split('_')[0]
		pepf = pepfiledict[spns]
		cdsf = cdsfiledict[spns]
		spnoder.append(spns)
		beckid = genid[:-3] + "." + genid[-2:] if genid[-3] == "X" else genid[:-2] + "." + genid[-1:]
		if seq[-1] in ["*","."]:
			seq = seq[:-1]
		if genid in pepf.references or beckid in pepf.references:
			pep1 = pepf.fetch(genid) if genid in pepf.references else  pepf.fetch(beckid)
			cds1 = cdsf.fetch(genid) if genid in cdsf.references else  cdsf.fetch(beckid)
			if pep1[-1] in ["*","."]:
				pep1 = pep1[:-1]
		else:
			pep1 = ""
			cds1 = ""
		if seq == pep1:
			recdsseq[spns] = cds1
		else:
			checkwrit.append(spns)
	if checkwrit:
		Checkresult.write('%s\t%s\n'%(names,';'.join(checkwrit)))
	else:
		outcdsfile = open('CDS_modify/%s.cds.fa'%(names),'w')
		print('CDS_modify/%s.cds.fa'%(names))
		for wspn in spnoder:
			outcdsfile.write('>%s\n%s\n'%(wspn,recdsseq[wspn]))
		outcdsfile.close()
Checkresult.close()
