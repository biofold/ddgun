#!/usr/bin/env python
##  Copyright (C) 2019  Ludovica Montanucci, Emidio Capriotti and Piero Fariselli
##
##  This program and all program in this package are free software;
##  you can redistribute it and/or modify it under the terms of the
##  GNU General Public License as published by the Free Software
##  Foundation; either version 2 of the License, or (at your option)
##  any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA
##
##  REFERENCES
##
##  Montanucci L, Capriotti E, Frank Y, Ben-Tal N, Fariselli P. (2019).
##  DDGun: an untrained method for the prediction of protein stability
##  changes upon single and multiple point variations.
##  BMC Bioinformatics. 20 (Suppl 14): 335. PMID:31266447
##


import os, sys, pickle, tempfile, argparse
from commands import getstatusoutput
import numpy as np

global pblast, uniref90, pprof, prog_path, data_path, aalist
prog_path=os.path.dirname(os.path.abspath(__file__))
data_path=prog_path+'/data'
tool_path=prog_path+'/tools'
util_path=prog_path+'/utils'
data_path=prog_path+'/data'
sys.path.append(tool_path)

from hsspTools import readHSSP, hssp2dic

aalist='ARNDCQEGHILKMFPSTWYV'
pprof=tool_path+'/ali2prof.py'
pblast=util_path+'/hh-suite/bin/hhblits'
uniref90=data_path+'/uniclust30_2018_08/uniclust30_2018_08'



def get_options():
	global uniref90, pblast
	parser = argparse.ArgumentParser(description='Program for generating protein mutations features.')
	parser.add_argument ('seqfile', type=str)
	parser.add_argument ('mutations')
	parser.add_argument ("--aa1", "--aaindex1", action="store",type=str, dest="aa1", help="Potential aaindex1")
	parser.add_argument ("--aa2", "--aaindex2", action="store",type=str, dest="aa2", help="Potential aaindex2")
	parser.add_argument ("--aa3", "--aaindex3", action="store",type=str, dest="aa3", help="Potential aaindex3")
	parser.add_argument ("-w", "--win", action="store",type=int, dest="win", help="Windows around mutation")
	parser.add_argument ("-o", "--out-file", action="store",type=str, dest="outfile", help="Output file")
	parser.add_argument ("-m", "--mutation-list",action="store_true",dest="ml", help="Mutation is parsed as a comma separated list")
	parser.add_argument ("-v", "--verbose", action="store",type=int, dest="verb", help="Verbose output")
	parser.add_argument ("--outdir", "--out-dir", action="store",type=str, dest="outdir", help="Output directory")
	parser.add_argument ("-d", "--db", action="store",type=str, dest="dbfile", help="DB file for hhblits")
	parser.add_argument ("-s", "--search-prog", action="store",type=str, dest="hhblits", help="hhblits")
	args = parser.parse_args()
	aa1='KYTJ820101'
	aa2='HENS920102'
	aa3='SKOJ970101'
	win=2
	muts={}
	sep=','
	verb=0
	outdir=None
	outfile=None
	seqfile=args.seqfile
	if os.path.isfile(seqfile)==False:	
		print >> sys.stderr,'ERROR: Incorrect sequence file '+seqfile+'.'
		sys.exit(1)
	if args.aa1: aa1=args.aa1
	if args.aa2: aa2=args.aa2
	if args.aa3: aa1=args.aa3
	if args.win: win=args.win
	if args.verb in [1,2]: verb=args.verb
	if args.outdir: outdir=args.outdir
	if args.outfile: outfile=args.outfile
	if args.dbfile: uniref90=args.dbfile
	if args.hhblits: pblast=args.hhblits
	if not os.path.isfile(pblast):
		print >> sys.stderr,'ERROR: hhblits program not found in',pblast
		sys.exit(4)
	if not os.path.isfile(uniref90+'_a3m_db.index'):
		print >> sys.stderr,'ERROR: DB file clust30_2018_08 not found in',uniref90
		sys.exit(5)
	if args.ml:
		if parse_mut(args.mutations): 
			muts[sort_mut(args.mutations)]=[mut for mut in sort_mut(args.mutations).split(sep) if mut!='']
	else:
		if os.path.isfile(args.mutations):
			lmut=open(args.mutations).read()
			muts=dict((sort_mut(mut),sort_mut(mut).split(sep)) for mut in lmut.replace(' ','').split('\n') if parse_mut(mut))
	if len(muts)==0:
		print >> sys.stderr,'ERROR: Incorrect mutation list.'
		sys.exit(2)
	return seqfile,muts,[aa1,aa2,aa3],win,verb,outfile,outdir


def parse_mut(imut,sep=','):
	v_mut=imut.split(sep)
	v_pos=[]
	for mut in v_mut:
		c=True
		try:
			pos=int(mut[1:-1])
			if pos in v_pos:
				c=False
			else:
				v_pos.append(pos)
			if aalist.index(mut[0])==-1 or aalist.index(mut[-1])==-1: c=False
		except:
			c=False
		if not c:
			if mut!='': print >> sys.stderr,'WARNING: Incorrect mutation',imut 
			break
	return c


def sort_mut(imut,sep=','):
	v_mut=imut.split(sep)
	t_mut=[(int(j[1:-1]),j) for j in v_mut]
	t_mut.sort()
	return ','.join([i[1] for i in t_mut])


def ali2fasta(filein,fileout):
	fb=open(filein)
	vd=''
	for line in fb:
		v=line.rstrip().split()
		vd=vd+'>'+v[0]+'\n'+v[1]+'\n'
	fb.close()
	fb=open(fileout,'w')
	fb.write(vd)
	fb.close()
		
		
def get_hssp(hsspfile):
	hssp=readHSSP(hsspfile)
	dhssp=hssp2dic(hssp)
	return dhssp


def get_pot_res(res,pot='KYTJ820101'):
	with open(data_path+'/aaindex1.pkl') as f:
		dpot=pickle.load(f).get(pot,{})
	return dpot.get(res,0.0)


def get_pot_prof(hsspfile,l_mut,pot='KYTJ820101'):
	l_score={}
	l_hssp={}
	with open(data_path+'/aaindex1.pkl') as f:
		dpot=pickle.load(f).get(pot,{})
	if len(dpot.keys())==0:
		print >> sys.stderr,'Incorrect potential',pot
		return l_score
	dhssp=get_hssp(hsspfile)
	for i in l_mut.keys():
		v_mut=l_mut[i]
		v_dat=[]
		v_prof=[]
		for mut in v_mut:
			swt=0.0
			snew=0.0
			wt=mut[0]
			pos=int(mut[1:-1])
			new=mut[-1]
			prof=dhssp.get(pos,{})
			if len(prof.keys())==0 or prof['WT']!=wt: 
				print >> sys.stderr,'WARNING: Profile position',pos,'not found or incorrect residue.'
				v_dat=[]
				v_prof=[]
				break
			v_prof.append((prof[wt],prof[new]))
			swt=dpot.get(wt,0.0)*prof.get(wt,0.0)*0.01
			snew=dpot.get(new,0.0)*prof.get(new,0.0)*0.01
			v_dat.append((swt,snew))	
		l_score[i]=v_dat
		l_hssp[i]=v_prof
	return l_score,l_hssp
	

def get_subs_prof(hsspfile,l_mut,pot='HENS920102'):
	l_score={}
	with open(data_path+'/aaindex2.pkl') as f:
		dpot=pickle.load(f).get(pot,{})
	if len(dpot.keys())==0:
		print >> sys.stderr,'Incorrect potential',pot
		return l_score
	dhssp=get_hssp(hsspfile)
	for i in l_mut.keys():
		v_mut=l_mut[i]
		v_dat=[]
		for mut in v_mut:
			swt=0.0
			snew=0.0
			wt=mut[0]
			pos=int(mut[1:-1])
			new=mut[-1]
			prof=dhssp.get(pos,{})
			if len(prof.keys())==0 or prof['WT']!=wt: 
				print >> sys.stderr,'WARNING: Profile position',pos,'not found or incorrect residue.'
				v_dat=[]
				break
			for aa in aalist:
				swt=swt+dpot.get((wt,aa),0.0)*prof.get(aa,0.0)*0.01
				snew=snew+dpot.get((new,aa),0.0)*prof.get(aa,0.0)*0.01
			v_dat.append((swt,snew))
		l_score[i]=v_dat
	return l_score
	

def get_seq_prof(hsspfile,l_mut,w=2,pot='SKOJ970101'):
	l_score={}
	with open(data_path+'/aaindex3.pkl') as f:
		dpot=pickle.load(f).get(pot,{})
	if len(dpot.keys())==0:
		print >> sys.stderr,'Incorrect potential',pot
		return l_score
	dhssp=get_hssp(hsspfile)
	n=len(dhssp.keys())
	for i in l_mut.keys():
		v_mut=l_mut[i]
		v_dat=[]
		for mut in v_mut:
			swt=0.0
			snew=0.0
			wt=mut[0]
			pos=int(mut[1:-1])
			new=mut[-1]
			prof=dhssp.get(pos,{})
			if len(prof.keys())==0 or prof['WT']!=wt: 
				print >> sys.stderr,'WARNING: Profile position',pos,'not found or incorrect residue.'
				v_dat=[]
				break
			s=max(1,pos-w)
			e=min(n,pos+w)
			for j in range(s,pos)+range(pos+1,e+1):
				iprof=dhssp.get(j,{})
				for aa in aalist:
					swt=swt+dpot.get((aa,wt),0.0)*iprof.get(aa,0.0)*0.01
					snew=snew+dpot.get((aa,new),0.0)*iprof.get(aa,0.0)*0.01
			v_dat.append((swt,snew))
		l_score[i]=v_dat
	return l_score


def run_seq_pipeline(seqfile,blast_prog=pblast,db=uniref90,outdir=None,e=1e-9):
	if outdir:
		tmpdir=outdir
		rd=''
	else:
		tmpdir=tempfile.mkdtemp()
		rd=tmpdir
	seqname=seqfile.split('/')[-1]
	blastfile=tmpdir+'/'+seqname+'.blast'
	hsspfile=tmpdir+'/'+seqname+'.hssp'
	if os.path.isfile(blastfile)==False:
		#cmd=blast_prog+' -i '+seqfile+' -d '+db+' -e '+str(e)+' -j 1 -b 1000 -v 1000 -o '+blastfile
		cmd=blast_prog+' -d  '+db+'  -i '+seqfile+' -opsi '+blastfile+'x  -n 2 -cpu 4 '
		print >> sys.stderr,'1) Run HHBLITS Search'
		print >> sys.stderr,cmd
		out=getstatusoutput(cmd)
		if out[0]!=0:
			blastfile=''
			print >> sys.stderr,'HHBLITS_ERROR:'+out[1]
			sys.exit(1)
		ali2fasta(blastfile+'x',blastfile)
		getstatusoutput('rm '+blastfile+'x')
	if os.path.isfile(hsspfile)==False:
		cmd=pprof+' '+seqfile+' '+blastfile+' '+hsspfile
		print >> sys.stderr,'2) Generate HSSP File'
		print >> sys.stderr,cmd
		out=getstatusoutput(cmd)
		if out[0]!=0:
			print >> sys.stderr,'HSSP_ERROR:'+out[1]
			#getstatusoutput('rm -r '+blastfile+' '+rd)
			sys.exit(1)
	return hsspfile
	

def get_muts_score(seqfile,hsspfile,muts,pots,win,outdir=None):
	l_data={}
	l_mut=muts
	s_hyd,l_hssp=get_pot_prof(hsspfile,l_mut,pots[0])
	s_sub=get_subs_prof(hsspfile,l_mut,pots[1])
	s_pro=get_seq_prof(hsspfile,l_mut,win,pots[2])
	for i in l_mut.keys():
		v_mut=l_mut[i]
		n=len(v_mut)
		hs=s_hyd.get(i,[])
		ss=s_sub.get(i,[])
		ps=s_pro.get(i,[])
		if len(hs)==0 or len(ss)==0 or len(ps)==0:
			print >> sys.stderr,'WARNING: Incorrect profile calculation for mutation',i
			continue
		l_data[i]=[]
		for j in range(n):
			v_score=[hs[j][1]-hs[j][0],ss[j][1]-ss[j][0],ps[j][1]-ps[j][0]]
			l_data[i].append(v_score)
	if not outdir:
		odir=os.path.dirname(hsspfile)
		#getstatusoutput('rm -r '+odir)
	return l_data,l_hssp
		
		
def print_data(seqfile,l_data,l_hssp,verb,sep=','):
	# Coefficients
	#0.331 0.267 -0.328
	nfile=seqfile.split('/')[-1]
	s_mut=[]
	out_data=[]
	header='#SEQFILE\tVARIANT\tS_DDG\tT_DDG\n'
	if verb==1: header='#SEQFILE\tVARIANT\t\tS_KD\tS_BL\tS_PROF\tDDG\tT_DDG\n'
	if verb==2: header='#SEQFILE\tVARIANT\tCONSERVATION\tS_KD\tS_BL\tS_PROF\tDDG\tT_DDG\n'
	for mut in l_data.keys():
		s_mut.append([[int(i[1:-1]) for i in mut.split(sep)],mut])
	s_mut.sort()
	
	for lpos,mut in s_mut:
		pm=[]
		n=len(lpos)
		v=[[],[],[],[]]
		line='\t'.join([nfile,mut])
		for j in range(n):
			vm=l_data[mut][j]
		 	pred=0.27*vm[0]+0.30*vm[1]-0.43*vm[2]
			pm.append(pred)
			for k in range(3): v[k].append('%.3f' %vm[k])
			v[3].append('|'.join([str(f) for f in l_hssp[mut][j]]))
			#print line+'\t'+str(j+1)+'\t'+'\t'.join([str(round(i,3)) for i in vm])+'\t'+str(round(pred,1))
		if len(pm)==1:
			mpred=pm[0]
		else:
			mpred=max(pm)+min(pm)-sum(pm)/float(n)
		sdata=','.join(v[0])+'\t'+','.join(v[1])+'\t'+','.join(v[2])
		sext=','.join(v[3])
		spred=','.join([str(round(i,1)) for i in pm ])+'\t'+str(round(mpred,1))
		if verb==1: line=line+'\t'+sdata
		if verb==2: line=line+'\t'+sext+'\t'+sdata
		out_data.append(line+'\t'+spred+'\n')
	if len(out_data)>0: out_data=[header]+out_data
	return out_data



if __name__ == '__main__':
	seqfile,muts,pots,win,verb,outfile,outdir=get_options()
	hsspfile=run_seq_pipeline(seqfile,pblast,uniref90,outdir)
	l_data,l_hssp=get_muts_score(seqfile,hsspfile,muts,pots,win,outdir)
	if len(l_data)==0:
		print >> sys.stderr,'ERROR: Incorrect mutation list.'
		sys.exit()
	out_data=print_data(seqfile,l_data,l_hssp,verb)
	if len(out_data)==0:
		print >> sys.stderr,'ERROR: No predictions returned. Check your input.'
		sys.exit(1)
	if not outfile:
		for line in out_data: print line.rstrip()
	else:
		try:
			fout=open(outfile,'w')
			fout.writelines(out_data)
			fout.close()
		except:
			print >> sys.stderr,'ERROR: File',outfile,'can not be saved.'
	
