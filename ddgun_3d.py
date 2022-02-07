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

from __future__ import print_function
import os, sys, pickle, tempfile, argparse
try:
        from subprocess import getstatusoutput
except:
        from commands import getstatusoutput
import numpy as np

global pblast, uniref90, at, pdssp, pprof, prog_path, data_path, aalist
prog_path=os.path.dirname(os.path.abspath(__file__))
data_path=prog_path+'/data'
tool_path=prog_path+'/tools'
util_path=prog_path+'/utils'
data_path=prog_path+'/data'
sys.path.append(tool_path)


import pdbTools as pdb
import dsspTools as dssp
from hsspTools import readHSSP, hssp2dic

aalist='ARNDCQEGHILKMFPSTWYV'
pprof=tool_path+'/ali2prof.py'
pblast=util_path+'/hh-suite/bin/hhblits'
pdssp=util_path+'/dssp/dsspcmbi'
uniref90=data_path+'/uniclust30_2018_08/uniclust30_2018_08'
at=pdb.HeavyAtomType


def get_options():
        global uniref90, pblast
        parser = argparse.ArgumentParser(description='Program for generating protein mutations features.')
        #parser = OptionParser('\nUsage : python %prog  pdbfile chain mutation []')
        parser.add_argument ('pdbfile', type=str)
        parser.add_argument ('chain')
        parser.add_argument ('mutations')
        parser.add_argument ("-o", "--out-file", action="store",type=str, dest="outfile", help="Output file")
        parser.add_argument ("--aa1", "--aaindex1", action="store",type=str, dest="aa1", help="Potential aaindex1")
        parser.add_argument ("--aa2", "--aaindex2", action="store",type=str, dest="aa2", help="Potential aaindex2")
        parser.add_argument ("--aa3", "--aaindex3", action="store",type=str, dest="aa3", help="Potential aaindex3")
        parser.add_argument ("--raa3", "--raaindex3", action="store",type=str, dest="raa3", help="Potential 3D aaindex3")
        parser.add_argument ("-w", "--win", action="store",type=int, dest="win", help="Windows around mutation")
        parser.add_argument ("-m", "--mutation-list",action="store_true",dest="ml", help="Mutation is parsed as a comma separated list")
        parser.add_argument ("--dmin", action="store",type=float, dest="dmin", help="Minimum radius around mutation")
        parser.add_argument ("--dmax", action="store",type=float, dest="dmax", help="Maximum radius around mutation")
        parser.add_argument ("-v", "--verbose", action="store",type=int, dest="verb", help="Verbose output")
        parser.add_argument ("--outdir", "--out-dir", action="store",type=str, dest="outdir", help="Output directory")
        parser.add_argument ("-d", "--db", action="store",type=str, dest="dbfile", help="DB file for hhblits")
        parser.add_argument ("-s", "--search-prog", action="store",type=str, dest="hhblits", help="hhblits")
        args = parser.parse_args()
        muts={}
        sep=','
        outdir=None
        outfile=None
        aa1='KYTJ820101'
        aa2='HENS920102'
        aa3='SKOJ970101'
        raa3='BASU010101'
        win=2
        verb=0
        dmin=0.0
        dmax=5.0
        pdbfile=args.pdbfile
        chain=args.chain
        if os.path.isfile(pdbfile)==False:        
                print('ERROR: Incorrect PDB file '+pdbfile+'.', file=sys.stderr)
                sys.exit(1)
        if args.aa1: aa1=args.aa1
        if args.aa2: aa2=args.aa2
        if args.aa3: aa3=args.aa3
        if args.raa3: raa3=args.raa3
        if args.win: win=args.win
        if args.dmin: dmin=args.dmin
        if args.dmax: dmax=args.dmax
        if args.verb in [1,2]: verb=args.verb
        if args.outdir: outdir=args.outdir
        if args.outfile: outfile=args.outfile
        if args.dbfile: uniref90=args.dbfile
        if args.hhblits: pblast=args.hhblits
        if not os.path.isfile(pblast):
                print('ERROR: hhblits program not found in',pblast, file=sys.stderr)
                sys.exit(4)
        if not os.path.isfile(uniref90+'_a3m_db.index'):
                print('ERROR: DB file clust30_2018_08 not found in',uniref90, file=sys.stderr)
                sys.exit(5)
        if args.ml:
                if args.mutations.lower() == 'all':
                        cmuts=get_all(pdbfile,chain)
                elif args.mutations.lower() == 'alascan':
                        cmuts=get_alascan(pdbfile,chain)
                else:
                        cmuts=expand_muts([sort_mut(mut) for mut in args.mutations.split(sep) if mut!=""])
                muts=dict((sort_mut(mut),sort_mut(mut).split(sep)) for mut in cmuts if parse_mut(mut))
                if len(muts)==0:
                        print('ERROR: Incorrect mutation list.', file=sys.stderr)
                        sys.exit(2)
                #if parse_mut(args.mutations):
                #        cmuts=expand_muts(sort_mut(args.mutations).split(sep)) 
                #        muts[sort_mut(args.mutations)]=[mut for mut in cmuts if mut!='']
        else:
                if os.path.isfile(args.mutations):
                        lmut=open(args.mutations).read()
                        cmuts=expand_muts(lmut.replace(' ','').split('\n'))
                        muts=dict((sort_mut(mut),sort_mut(mut).split(sep)) for mut in cmuts if parse_mut(mut))
        if len(muts)==0:
                print('ERROR: Incorrect mutation list.', file=sys.stderr)
                sys.exit(2)
        if dmin>=dmax:
                print('ERROR: Incorrect radius distance range ',dmin,dmax, file=sys.stderr)
                sys.exit(3)
        return pdbfile,chain,muts,[aa1,aa2,aa3,raa3],[dmin,dmax],win,verb,outfile,outdir
    

def get_pdbchain(pdbfile,chain,atoms=at):
        och=pdb.readPDB(pdbfile,chain,atoms)
        return och 


def parse_mut(imut,sep=','):
        v_mut=imut.split(sep)
        v_pos=[]
        for mut in v_mut:
                c=True
                try:
                        #pos=int(mut[1:-1])
                        pos=mut[1:-1]
                        if pos in v_pos:
                                c=False
                        else:
                                v_pos.append(pos)
                        if aalist.find(mut[0])==-1 or aalist.find(mut[-1])==-1: c=False
                except:
                        c=False
                if not c:
                        if mut!='': print('WARNING: Incorrect mutation',imut, file=sys.stderr) 
                        break
        return c


def sort_mut(imut,sep=','):
        v_mut=imut.split(sep)
        try:
                t_mut=[(int(j[1:-1]),j) for j in v_mut]
        except:
                t_mut=[(j[1:-1],j) for j in v_mut]
        t_mut.sort()
        return sep.join([i[1] for i in t_mut])
        

def expand_muts(lmut):
        vmut=[]
        for mut in lmut:
                if mut=='': continue
                if aalist.find(mut[-1])!=-1: vmut.append(mut)
                if mut[-1]=="*":
                        for aa in aalist:
                                if aa!=mut[0]: vmut.append(mut[:-1]+aa)
                else:
                        pass
        return vmut


def get_alascan(pdbfile,chain):
        cmut=[]
        och=pdb.readPDB(pdbfile,chain,['CA'])
        lres=och.getPDBNumbers()
        seq=och.getSequence()
        n=len(seq)
        for i in range(n):
                if seq[i]=='A':
                        continue
                else:
                        for aa in aalist:
                                if aa!=seq[i]: cmut.append(seq[i]+lres[i]+aa)
        return cmut
        

def get_all(pdbfile,chain):
        cmut=[]
        och=pdb.readPDB(pdbfile,chain,['CA'])
        lres=och.getPDBNumbers()
        seq=och.getSequence()
        n=len(seq)
        for i in range(n):
                for aa in aalist:
                        if aa!=seq[i]: cmut.append(seq[i]+lres[i]+aa)
        return cmut

        
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


def get_dssp_rsa(dsspfile,chain,l_mut):
        l_rsa={}
        odssp=dssp.getAccfromNumbers(dsspfile,chain)
        for i in list(l_mut.keys()):
                v_mut=l_mut[i]
                i_rsa=[]
                for mut in v_mut:
                        rsa=None
                        wt=mut[0]
                        pos=mut[1:-1]
                        rdssp=odssp.get(pos,[])
                        if len(rdssp)==3 and rdssp[0]==wt: rsa=rdssp[2]
                        i_rsa.append(rsa)
                if None not in i_rsa:
                        l_rsa[i]=i_rsa
        return l_rsa
                
                
def get_hssp(hsspfile):
        hssp=readHSSP(hsspfile)
        dhssp=hssp2dic(hssp)
        return dhssp


def get_env_residue(pdbfile,chain,l_mut,d=[0.0,5.0],atoms=at):
        lres_env={}
        lnres={}
        och=pdb.readPDB(pdbfile,chain,atoms)
        lres=och.getPDBNumbers()
        for i in list(l_mut.keys()):
                v_mut=l_mut[i]
                v_dat=[]
                v_res=[]
                for mut in v_mut:
                        res_env=[]
                        nres=None
                        wt=mut[0]
                        pos=mut[1:-1]
                        new=mut[-1]
                        r_i=och.getResidueByNumber(pos)
                        if not r_i or r_i.NameOne!=wt:
                                print('WARNING: Incorrect residue',mut[:-1], file=sys.stderr)
                                v_dat=[]
                                v_res=[]
                                break
                        nres=lres.index(pos)+1
                        for j in lres:
                                if j==pos: continue
                                r_j=och.getResidueByNumber(j)        
                                dist=r_i.residueDistance(r_j)
                                if dist>=d[0] and dist<=d[1]:
                                        res_env.append((r_j.NameOne,j))
                        v_res.append(wt+str(nres)+new)
                        v_dat.append(res_env)
                if len(v_dat)>0: lres_env[i]=v_dat
                if len(v_res)>0: lnres[i]=v_res
        return lres,lnres,lres_env


def get_pot_res(res,pot='KYTJ820101'):
        dpot=pickle.load(open(data_path+'/aaindex1.pkl','rb')).get(pot,{})
        return dpot.get(res,0.0)


def get_pot_3d(pdbfile,chain,hsspfile,l_mut,d=[0.0,5.0],pot='BASU010101'):
        l_score={}
        l_nres={}
        n=len(l_mut)
        dpot=pickle.load(open(data_path+'/aaindex3.pkl','rb')).get(pot,{})
        if len(list(dpot.keys()))==0:
                print('Incorrect potential',pot, file=sys.stderr)
                return l_nres,l_score
        dhssp=get_hssp(hsspfile)
        l_pos,l_nres,lres_env=get_env_residue(pdbfile,chain,l_mut,d)
        if len(l_nres)==0:
                print('ERROR: Incorrect mutation list.', file=sys.stderr)
                sys.exit(0)
        for i in list(l_mut.keys()):
                v_mut=l_mut[i]
                v_dat=[]
                res_env=lres_env.get(i,[])
                if len(res_env)==0:
                        print('WARNING: No residues around',i, file=sys.stderr)
                        continue
                for j in range(len(v_mut)):
                        mut=v_mut[j]
                        swt=0.0
                        snew=0.0
                        wt=mut[0]
                        pos=mut[1:-1]
                        new=mut[-1]
                        if not res_env[j]:
                                print('WARNING: No residues around',wt+pos, file=sys.stderr)
                                v_dat=[]
                                break
                        #print mut,res_env[j]
                        for l in res_env[j]:
                                try:
                                        prof=dhssp.get(l_pos.index(l[1])+1,{})
                                except:
                                        continue
                                for k in aalist:
                                        swt=swt+0.01*prof.get(k,0.0)*dpot.get((wt,k),0.0)
                                        snew=snew+0.01*prof.get(k,0.0)*dpot.get((new,k),0.0)
                        v_dat.append((swt,snew))
                l_score[i]=v_dat
        return l_nres,l_score,lres_env        



def get_pot_prof(hsspfile,l_mut,pot='KYTJ820101'):
        l_score={}
        l_hssp={}
        dpot=pickle.load(open(data_path+'/aaindex1.pkl','rb')).get(pot,{})
        if len(list(dpot.keys()))==0:
                print('Incorrect potential',pot, file=sys.stderr)
                return l_score
        dhssp=get_hssp(hsspfile)
        for i in list(l_mut.keys()):
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
                        if len(list(prof.keys()))==0 or prof['WT']!=wt: 
                                print('WARNING: Profile position',pos,'not found or incorrect residue.', file=sys.stderr)
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
        dpot=pickle.load(open(data_path+'/aaindex2.pkl','rb')).get(pot,{})
        if len(list(dpot.keys()))==0:
                print('Incorrect potential',pot, file=sys.stderr)
                return l_score
        dhssp=get_hssp(hsspfile)
        for i in list(l_mut.keys()):
                v_mut=l_mut[i]
                v_dat=[]
                for mut in v_mut:
                        swt=0.0
                        snew=0.0
                        wt=mut[0]
                        pos=int(mut[1:-1])
                        new=mut[-1]
                        prof=dhssp.get(pos,{})
                        if len(list(prof.keys()))==0 or prof['WT']!=wt: 
                                print('WARNING: Profile position',pos,'not found or incorrect residue.', file=sys.stderr)
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
        dpot=pickle.load(open(data_path+'/aaindex3.pkl','rb')).get(pot,{})
        if len(list(dpot.keys()))==0:
                print('Incorrect potential',pot, file=sys.stderr)
                return l_score
        dhssp=get_hssp(hsspfile)
        n=len(list(dhssp.keys()))
        for i in list(l_mut.keys()):
                v_mut=l_mut[i]
                v_dat=[]
                for mut in v_mut:
                        swt=0.0
                        snew=0.0
                        wt=mut[0]
                        pos=int(mut[1:-1])
                        new=mut[-1]
                        prof=dhssp.get(pos,{})
                        if len(list(prof.keys()))==0 or prof['WT']!=wt: 
                                print('WARNING: Profile position',pos,'not found or incorrect residue.', file=sys.stderr)
                                v_dat.append(None)
                                continue
                        s=max(1,pos-w)
                        e=min(n,pos+w)
                        for j in list(range(s,pos))+list(range(pos+1,e+1)):
                                iprof=dhssp.get(j,{})
                                for aa in aalist:
                                        swt=swt+dpot.get((aa,wt),0.0)*iprof.get(aa,0.0)*0.01
                                        snew=snew+dpot.get((aa,new),0.0)*iprof.get(aa,0.0)*0.01
                        v_dat.append((swt,snew))
                l_score[i]=v_dat
        return l_score


def run_3d_pipeline(pdbfile,chain,blast_prog=pblast,db=uniref90,outdir=None,atoms=at,dssp_prog=pdssp,e=1e-9):
        err=0
        seq=''
        if outdir:
                tmpdir=outdir
                rd=''
        else:
                tmpdir=tempfile.mkdtemp()
                rd=tmpdir
        pdbname=os.path.abspath(pdbfile).split('/')[-1]
        chainfile=tmpdir+'/'+pdbname+'.'+chain
        dsspfile=chainfile+'.dssp'
        seqfile=chainfile+'.fasta'
        blastfile=chainfile+'.blast'
        hsspfile=chainfile+'.hssp'
        och=get_pdbchain(pdbfile,chain,atoms)
        och.writePDB(chainfile,[],'Y',chain,'Y')
        seq=och.getSequence()
        if seq=='':
                print('ERROR: Incorrect PDB file',pdbfile,'or chain',chain, file=sys.stderr)
                getstatusoutput('rm -r '+chainfile+' '+rd)
                sys.exit(1)
        f=open(seqfile,'w')
        f.write('>'+pdbname+'.'+chain+'\n'+seq)
        f.close()
        if os.path.isfile(chainfile)==False:
                print('ERROR: Chain file',chain,'not found.', file=sys.stderr)
                sys.exit(2)
        if os.path.isfile(dsspfile)==False: 
                cmd=dssp_prog+' '+chainfile+' '+dsspfile
                print('1) Generate DSSP File', file=sys.stderr)
                print(cmd, file=sys.stderr)
                out=getstatusoutput(cmd)
        if os.path.isfile(dsspfile)==False:
                print('ERROR: DSSP file',chain,'not found.', file=sys.stderr)
                getstatusoutput('rm -r '+chainfile+' '+seqfile+' '+rd)
                sys.exit(3)
        if os.path.isfile(blastfile)==False:
                #cmd=blast_prog+' -i '+seqfile+' -d '+db+' -e '+str(e)+' -j 1 -b 1000 -v 1000 -o '+blastfile
                cmd=blast_prog+' -d  '+db+'  -i '+seqfile+' -opsi '+blastfile+'x  -n 2 -cpu 4 '
                print('2) Run HHBLITS Search', file=sys.stderr)
                print(cmd, file=sys.stderr)
                out=getstatusoutput(cmd)
                if out[0]!=0:
                        print('HHBLITS_ERROR:'+out[1], file=sys.stderr)
                        getstatusoutput('rm -r '+chainfile+' '+seqfile+' '+dsspfile+' '+rd)
                        sys.exit(4)
                ali2fasta(blastfile+'x',blastfile)
                getstatusoutput('rm '+blastfile+'x')
        if os.path.isfile(hsspfile)==False:
                cmd=pprof+' '+seqfile+' '+blastfile+' '+hsspfile
                print('3) Generate HSSP File', file=sys.stderr)
                print(cmd, file=sys.stderr)
                out=getstatusoutput(cmd)
                if out[0]!=0:
                        print('HSSP_ERROR:'+out[1], file=sys.stderr)
                        getstatusoutput('rm -r '+chainfile+' '+seqfile+' '+dsspfile+' '+blastfile+' '+rd)
                        sys.exit(5)
        return chainfile,dsspfile,hsspfile


def get_muts_score(chainfile,chain,dsspfile,hsspfile,muts,pots,d,win=2,outdir=None):
        l_data={}
        #l_mut=list(set(muts))
        l_mut=muts
        l_nmut,s_3d,lres_env=get_pot_3d(chainfile,chain,hsspfile,l_mut,d,pots[3])
        d_rsa=get_dssp_rsa(dsspfile,chain,l_mut)
        s_hyd,l_hssp=get_pot_prof(hsspfile,l_nmut,pots[0])
        s_sub=get_subs_prof(hsspfile,l_nmut,pots[1])
        s_pro=get_seq_prof(hsspfile,l_nmut,win,pots[2])
        for i in list(l_mut.keys()):
                v_mut=l_mut[i]
                n=len(v_mut)
                hs=s_hyd.get(i,[])
                ss=s_sub.get(i,[])
                ps=s_pro.get(i,[])
                s3=s_3d.get(i,[])
                rsa=d_rsa.get(i,[])
                if len(hs)==0 or len(ss)==0 or len(ps)==0:
                        print('WARNING: Incorrect profile calculation for mutation',i, file=sys.stderr)
                        continue
                if len(s3)==0 or len(rsa)==0:
                        print('WARNING: Incorrect structural features for mutation',i, file=sys.stderr)
                        continue
                l_data[i]=[]
                for j in range(n):
                        v_score=[hs[j][1]-hs[j][0],ss[j][1]-ss[j][0],ps[j][1]-ps[j][0],\
                                s3[j][1]-s3[j][0],0.0,rsa[j],0]
                        l_data[i].append(v_score)
        if not outdir: 
                odir=os.path.dirname(chainfile)        
                getstatusoutput('rm -r '+odir)
        return l_data,l_hssp,lres_env
                
                
def print_data(pdbfile,chain,l_data,l_hssp,lres_env,verb,sep=','):
        # Coefficients
        # 0.18 0.20 0.29 0.33
        nfile=pdbfile.split('/')[-1]
        s_mut=[]
        out_data=[]
        for mut in list(l_data.keys()):
                try:
                        s_mut.append([[int(i[1:-1]) for i in mut.split(sep)],mut])
                except:
                        s_mut.append([[i[1:-1] for i in mut.split(sep)],mut])
        s_mut.sort()
        header='#PDBFILE\tCHAIN\tVARIANT\tS_DDG\tT_DDG\n'
        if verb==1: header='#PDBFILE\tCHAIN\tVARIANT\t\tS_KD\tS_BL\tS_PROF\tS_3D[WT]\tRSA[WT]\tDDG\tT_DDG\n'
        if verb==2: header='#PDBFILE\tCHAIN\tVARIANT\tCONSERVATION\tCONTACTS\tS_KD\tS_BL\tS_PROF\tS_3D[WT]\tRSA[WT]\tDDG\tT_DDG\n'
        for lpos,mut in s_mut:
                pm=[]
                n=len(lpos)
                v=[[],[],[],[],[],[],[],[],[]]
                line='\t'.join([nfile,chain,mut])
                for j in range(n):
                        vm=l_data[mut][j]
                        f1=(1.1-vm[5]*0.01)
                        pred=(0.18*vm[0]+0.20*vm[1]-0.29*vm[2]-0.33*vm[3])*f1
                        pm.append(pred)
                        for k in list(range(4))+[5]: v[k].append('%.3f' %vm[k])
                        v[7].append('|'.join([str(f) for f in l_hssp[mut][j]]))
                        v[8].append('|'.join([r+p for r,p in lres_env[mut][j]]))
                        #print line+'\t'+str(j+1)+'\t'+'\t'.join([str(round(i,3)) for i in vm])+'\t'+str(round(pred,1))
                if len(pm)==1:
                        mpred=pm[0]
                else:
                        mpred=max(pm)+min(pm)-sum(pm)/float(n)
                sdata=','.join(v[0])+'\t'+','.join(v[1])+'\t'+','.join(v[2])+'\t'+\
                        ','.join(v[3])+'\t'+','.join(v[5])
                sext=','.join(v[7])+'\t'+','.join(v[8])
                spred=','.join([str(round(i,1)) for i in pm ])+'\t'+str(round(mpred,1))
                if verb==1: line=line+'\t'+sdata
                if verb==2: line=line+'\t'+sext+'\t'+sdata
                out_data.append(line+'\t'+spred+'\n')
        if len(out_data)>0: out_data=[header]+out_data
        return out_data



        
                

if __name__ == '__main__':
        pdbfile,chain,muts,pots,d,win,verb,outfile,outdir=get_options()
        chainfile,dsspfile,hsspfile=run_3d_pipeline(pdbfile,chain,pblast,uniref90,outdir)
        l_data,l_hssp,lres_env=get_muts_score(chainfile,chain,dsspfile,hsspfile,muts,pots,d,win,outdir)
        if len(l_data)==0: 
                print('ERROR: Incorrect mutation list.', file=sys.stderr)
                sys.exit()
        out_data=print_data(pdbfile,chain,l_data,l_hssp,lres_env,verb)
        if len(out_data)==0:
                print('ERROR: No predictions returned. Check your input.', file=sys.stderr)
                sys.exit(1)
        if not outfile:
                for line in out_data: print(line.rstrip())
        else:
                try:
                        fout=open(outfile,'w')
                        fout.writelines(out_data)
                        fout.close()
                except:
                        print('ERROR: File',outfile,'can not be saved.', file=sys.stderr)
