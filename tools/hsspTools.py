#!/usr/bin/env python
#
# 

import sys
import string, math

header_hssp=['SeqNo', 'PDBNo', 'V', 'L', 'I', 'M', 'F', 'W', 'Y', 'G', 'A', \
	'P', 'S', 'T', 'C', 'H', 'R', 'K', 'Q', 'E', 'N', 'D', 'X', \
	'DEL', 'INS', 'CONS', 'PTOT', 'WT']

aalist=['V', 'L', 'I', 'M', 'F', 'W', 'Y', 'G', 'A', \
        'P', 'S', 'T', 'C', 'H', 'R', 'K', 'Q', 'E', 'N', 'D']


def readHSSP(file,chain='_'):
   '''
   '''
   import re

   lines = open(file).readlines()
   #start_match="SeqNo PDBNo   V   L   I   M   F   W   Y   G   A   P   S   T   C   H   R   K   Q   E   N   D"
   #start_match="SeqNo PDBNo     V     L     I     M     F     W     Y     G     A     P     S     T     C     H     R     K     Q     E     N     D"
   #start_match="SeqNo PDBNo"
   start_match="SeqNo"
   hssp=[]
   filelen=len(lines)
   l=0
   pat=re.compile(start_match)
   while not pat.search(lines[l]):
      l=l+1  
   if l >= filelen:
       sys.stderr.write('error in file \n')
       sys.exit()
   pat=re.compile('\s+\d+\s+\d+')
   l=l+1
   line=lines[l]
   while l < filelen and line[0] != '#':
      line=lines[l]
#      pat=re.compile('^Sequences\s+\((\d+)\:(\d+)\)\s+Aligned\. Score:\s+(\d+)')
      m = pat.search(line)
      if (m) and ( chain =='_' or chain==line[11] ):
	    pos=line[:5].strip()
            pdbpos=line[5:12].strip()
            lst=[pos,pdbpos]+line[12:].split()
            #lst=line.split()
            hssp.append(lst)
      l=l+1	    
   return hssp 


def hssp2dic(hssp,lkey=header_hssp):
        ldic={}
	n=len(lkey)
        for i in hssp:
                dicv={}
                for j in range(n):
			try:
                        	#dicv[lkey[j]]=int(i[j])
				dicv[lkey[j]]=float(i[j])
			except:
				dicv[lkey[j]]=i[j]
                ldic[dicv['SeqNo']]=dicv
        return ldic


def CI(dhssp,pos,reslist):
        ci=0.0
        mdic=meanprofile(dhssp,reslist)
        for res in reslist:
                ci=ci+(float(dhssp[pos][res])-mdic[res])**2
        ci=math.sqrt(ci)
        return ci


def addCI(dhssp,reslist):
        mdic=meanprofile(dhssp,reslist)
        keys=dhssp.keys()
        for i in keys:
                ci=0.0
                for res in reslist:
                        ci=ci+(float(dhssp[i][res])-mdic[res])**2
                dhssp[i]['CI']=math.sqrt(ci)
        return dhssp


def meanprofile(dhssp,reslist):
        mdic={}
        keys=dhssp.keys()
	n=0
        for i in keys:
                for res in reslist:
                        if mdic.get(res,0)==0: mdic[res]=0
			ni=int(dhssp[i]['PTOT'])
                        mdic[res]=mdic[res]+int(dhssp[i][res])*ni
			n=n+ni
        for res in reslist:
                mdic[res]=float(mdic[res])/n
        return mdic


def getProfile(hsspfile, chain='_',slicelist=[]):
    ''' getProfile(hsspfile, chain='_',slicelist=[])
        returns the profile as taken from the hssp like file  
	if slicelist=[] return everything 
	else            return only the indices in the slice
    '''
    hsspcontent=readHSSP(hsspfile,chain)
    proflen=20
    end=6
    start=proflen+6
    prof=[]
    for v in hsspcontent:
        vtmp=[]
	for i in range(-start,-end,1):
#            vtmp.append(string.atoi(v[i]))
            vtmp.append(v[i])
	prof.append(vtmp)
    if not slicelist:
        return prof
    else:
	try:
            retprof=[prof[i] for i in slicelist]
	except:
	    print "hsspTools",hsspfile,i,len(prof)
        return retprof

# leggi il file hssp
if __name__=='__main__':
    if(len(sys.argv) == 1):
        print "syntax :",sys.argv[0]," hsspfile [chain]"
        sys.exit()
    else:
        try:
            if(len(sys.argv)>=3):
	        chain=sys.argv[2]
            else:	
	        chain='_'     
        except IOError, (errno, strerror):
            print sys.argv[1]
            print "I/O error(%s): %s" % (errno, strerror)


    hssp=readHSSP(sys.argv[1])
    d=hssp2dic(hssp)
    for i in hssp:
	print i
