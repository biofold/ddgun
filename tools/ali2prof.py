#!/usr/bin/env python
import sys, os, tempfile, numpy
from commands import getstatusoutput
from optparse import OptionParser
#from Bio.Blast import NCBIStandalone
from Bio.SearchIO._legacy import NCBIStandalone
from Bio import SeqIO
aalist='ARNDCQEGHILKMFPSTWYV'
aaprof='VLIMFWYGAPSTCHRKQENDBUZX'


def readconfs(confile='0-configure'):
	global confs
	confs={}
	progfile = os.path.dirname(__file__)
	progdir  = os.path.abspath(progfile)
	fileconf=progdir+'/'+confile
	if os.path.isfile(fileconf):
		lines=open(fileconf,'r').readlines()
		for line in lines:
			v=line.split('=')
			confs[v[0].strip()]=v[1].replace('\n','').strip()
	else:				
		sys.stderr.write('ERROR: Configuration file '+fileconf+' not found\n')
	return


# Functions for working with muscle multiple sequence alignments	

def get_sequence(fastafile):
        seq=''
        scode=''
        lines=open(fastafile,'r').readlines()
        for line in lines:
                if line[0]!='>':
                        seq=seq+line.replace('\n','')
                else:
                        #scode=line[1:].split()[0].split('|')[0]
                        scode=line[1:].split()[0]
        return scode,seq


def readfasta_align(protid,alifile):
	vali=[]
	protseq=''
	handle = open(alifile, "rU")
	alis=list(SeqIO.parse(handle, "fasta"))
	handle.close()
	for ali in alis:
		pid,seq=ali.id,str(ali.seq)
		if pid==protid or protid.find(pid)==0: 
			protseq=seq
		else:
			vali.append([pid,seq])
	if protseq!='': 
		vali=[[protid,protseq]]+vali
		protseq=protseq.replace('-','')
	else:
		print >> sys.stderr,'ERROR: Protein sequence ID not found in the alignment file',alifile
		sys.exit(1)
	return protseq,vali	
	
	
def init_profile(seq,amino=aaprof):
	# 20 standard AAs + 4 other AAs
	# 4 deletion, insertion, conservation, total
	na=len(amino)
	extra=4
	n=len(seq)
	vprof=numpy.zeros((n,na+extra))
	for i in range(n):
		apos=amino.find(seq[i])
		# assign position of X
		if apos==-1: apos=na-1
		vprof[i][apos]=1
	return vprof


def global_profile(seq,malign,amino):
	vprof=init_profile(seq,amino)
	for salign in malign[1:]:
		vprof=update_globalprofile(vprof,malign[0][1],salign[1],amino)
	return vprof
	

def update_globalprofile(profile,query_seq,subjct_seq,amino=aaprof):
	seq=query_seq.replace('-','')
	seq1=query_seq
	seq2=subjct_seq
	na=len(amino)
	start=1
	end=len(seq)
	alen=len(query_seq)
	c=start-1
	for i in range(alen):
		if seq1[i]!='-':
			#check position	
			if seq1[i]==seq[c] or seq1[i]=='X':
				if seq2[i]=='-':
					# Deletion position na
					profile[c][na]+=1
				else:
					apos=amino.find(seq2[i])
					if apos==-1: apos=na-1
					profile[c][apos]+=1
					# Total position na+3
					profile[c][na+3]+=1
					#print alignment[0],len(seq),len(profile.tolist()),len(seq1),len(seq2),c,i,alignment
					# Conservation position na+2
					if seq[c]==seq2[i]: profile[c][na+2]+=1
				c=c+1
		else:
			# Insertion position na+1
			profile[c-1][na+1]+=1
	#print c,end,len(alignment[1][2]),alignment
	return profile		


def write_plainhssp(seq,alignments,profile,amino=aaprof,norm=True,nums=[]):
	na=len(amino)
	nali=len(alignments)
	nres=len(seq)
	outprof=''
	#for ali in alignments:
	#	print ali[1]
	profline=''
	header='ID Pseudo hssp file.\nSEQ\n'
	aaline=' SeqNo PDBNo     V     L     I     M     F     W     Y     G     A     P     S     T     C     H     R     K     Q     E     N     D     X   DEL   INS  CONS  PTOT    WT\n'
	nalign='SEQ\nNALIGN  '+str(nali)+'\n## SEQUENCE PROFILE AND INFOS\n'
	for alignment in alignments:
		header=header+'SEQ %-64s...'%(alignment[0])+'\n'
	for i in range(nres):
		prof=profile[i]
		line='%5d %5d '%(i+1,i+1)
		if norm:
			tot=numpy.sum(prof[:20])
			if tot>0:
				nval=100./tot
			else:
				nval=1
		else:
			nval=1
		for j in range(20):
			line=line+'%6.1f' %(nval*prof[j])
		xres=numpy.sum(prof[20:na])
		aa=seq[i]
		rdel,rins,cons,tot=(prof[na],prof[na+1],prof[na+2],prof[na+3])
		#line=line+'%4d %4d %4d %5d %5d %4s' %(xres,rdel,rins,cons,tot,aa)
		line=line+'%6d %5d %5d %5d %5d %5s' %(xres,rdel,rins,cons,tot,aa)
		profline=profline+line+'\n'
	if  profline!='': outprof=header+nalign+aaline+profline+'//'
	return	outprof


def getfasta_profile(seqfile,alignfile,outfile,amino=aaprof):
	scode,seq=get_sequence(seqfile)
	protseq,alignments=readfasta_align(scode,alignfile)
	profile=global_profile(protseq,alignments,amino)
	outprof=write_plainhssp(protseq,alignments,profile,amino,norm=True,nums=[])
	if outprof!='': 
		open(outfile,'w').write(outprof)
		print 'ProtId:',scode
		print 'Alignment:',alignfile
		print 'Output:',outfile
	return	

	
def check_input(options,args):
	num=3
	err=0
	seqfile=''
	alignfile=''
	outfile=''
	if len(args)==num:
		seqfile=args[0]
        	alignfile=args[1]
        	outfile=args[2]
        	#readconfs()
		if not(os.path.isfile(seqfile)):
                        sys.stderr.write('ERROR: Alignment file '+seqfile+' not found \n')
                        err=1
        	#Check blastout file
        	if not(os.path.isfile(alignfile)): 
        		sys.stderr.write('ERROR: Alignment file '+alignfile+' not found \n')
        		err=1
        	#Check output directory
        	if not(os.path.isdir(os.path.abspath(os.path.dirname(outfile)))):
        		sys.stderr.write('ERROR: Output directory '+outfile+' not exist \n')
        		err=1
	else:
		err=1
		cmd='python program.py seqfile alignment outfile [-t blast_type]\n'
		sys.stderr.write('ERROR: Incorrect input\n'+cmd)	
	return seqfile,alignfile,outfile,err		


if (__name__  == '__main__'):
        parser = OptionParser('\nUsage : python %prog  seqfile alignment outfiles \n')
        (options, args) = parser.parse_args()
        # Defined default input
        seqfile,alignfile,outfile,err=check_input(options,args)
        if err==0:
                	getfasta_profile(seqfile,alignfile,outfile,amino=aaprof)
