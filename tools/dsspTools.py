#!/usr/bin/env python
#
# 

import sys

Norm_Acc={"A" :106.0,  "B" :160.0,         # D or N 
   "C" :135.0,  "D" :163.0,  "E" :194.0,
   "F" :197.0,  "G" : 84.0,  "H" :184.0,
   "I" :169.0,  "K" :205.0,  "L" :164.0,
   "M" :188.0,  "N" :157.0,  "P" :136.0,
   "Q" :198.0,  "R" :248.0,  "S" :130.0,
   "T" :142.0,  "V" :142.0,  "W" :227.0,
   "X" :180.0,         # undetermined (deliberate)
   "Y" :222.0,  "Z" :196.0}         # E or Q 


def readDSSP(file,chain='_'):
   '''
   '''
   import string

   lines = open(file).readlines()
   start_match="  #  RESIDUE AA STRUCTURE BP1 BP2  ACC"
   len_smatch=len(start_match)
   dssp=[]
   l=0
   while(lines[l][:len_smatch] != start_match):
      l=l+1  
   for line in lines[l+1:]:
      rchain=line[11]  # current chain
      if(rchain==chain or chain =='_'):
            resName=string.strip(line[13]) # residue
            if resName.islower(): #cysteines
	        resName='C'
            resNum=string.strip(line[5:11])  # current residue number
            resSS=line[16]
            resAcc=string.atoi(line[34:38])
	    resAngle=(float(string.strip(line[103:109])),float(string.strip(line[109:114])))
            dssp.append([resName,resNum,resSS,resAcc,resAngle])
   return dssp 


def getSSfromNumbers(dsspfile, chain='_'):
    ''' getSSfromNumbers(dsspfile, chain='_') 
        returns the hash pdbNumber : SecondaryStructure
    '''
    dsspcontent=readDSSP(dsspfile,chain)
    hashSS={}
    for resName,resNum,resSS,resAcc,resAngle in dsspcontent:
        if resSS == ' ':
            resSS ='C' 
        hashSS[resNum]=resSS
    return hashSS


def getSS(dsspfile, chain='_', ssUsed=['H','E']):
    ''' getSS(dsspfile, chain='_') 
        return the string of SecondaryStructure
    '''
    dsspcontent=readDSSP(dsspfile,chain)
    ss=''
    for resName,resNum,resSS,resAcc,resAngle in dsspcontent:
        if resSS in ssUsed:
            ss+=resSS
        else:
            ss+='C'
    return ss

def getAccfromNumbers(dsspfile,chain='_', th_acc=0):
	''' getAccfromNumbers(dsspfile,chain='_', th_acc=0)
	    returns an hash pdbNumber: [resName, resAcc, relAcc]
        '''
	dsspcontent=accNormList(dsspfile, chain, th_acc=0)
	hashAccNorm={}
	for resName,resNum,resSS,resAcc,relAcc in dsspcontent:
		hashAccNorm[resNum]=[resName,resAcc,relAcc]
	return hashAccNorm

def getPhiPsi(dsspfile,chain='_'):
	''' getPhiPsi(dsspfile, chain='_')
	    return a list of tuple of Phi and Psi angle
	'''
	dsspcontent=readDSSP(dsspfile,chain)
	angle=[]
	for resName,resNum,resSS,resAcc,resAngle in dsspcontent:
		angle.append(resAngle)
	return angle
		    

def accNormList(dsspfile, chain='_', th_acc=0):
    ''' accNormList(dsspfile, chain, th_acc=0)
        returns the list of the residue of chain whose accessibility is >= th_acc '''
    
    dsspcontent=readDSSP(dsspfile,chain)
    accreturn=[]
    for resName,resNum,resSS,resAcc,resAngle in dsspcontent:
	if resName!='!':
        	relAcc=int(100.0*resAcc/Norm_Acc[resName]+0.5)
        	relAcc=min(100,relAcc)
        	if relAcc>=th_acc:
	   		accreturn.append([resName,resNum,resSS,resAcc,relAcc])
    return accreturn
    

# leggi il file dssp
if __name__=='__main__':
    if(len(sys.argv) == 1):
        print "syntax :",sys.argv[0]," dsspfile [chain]"
        sys.exit()
    else:
        try:
            if(len(sys.argv)>=3):
	        chain=sys.argv[2]
            else:	
	        chain='_'     
            Protein = readDSSP(sys.argv[1],chain)
        except IOError, (errno, strerror):
            print sys.argv[1]
            print "I/O error(%s): %s" % (errno, strerror)
	
    for i in Protein:
        print i
    
#    acc= accNormList(sys.argv[1],chain,th_acc=16)
#    for i in acc:
#        print i
