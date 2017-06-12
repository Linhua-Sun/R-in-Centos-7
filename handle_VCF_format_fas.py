#!/Users/sunlinhua/anaconda/bin/python

## 2017-6-12 python script to handle VCF files to fasta

import sys

if len(sys.argv)==1:
	print('Please give me a datamash transpose VCF clean file path!')
	sys.exit()
	
f = open(sys.argv[1],'r')  
result = list() 
for line in f.readlines():
    line = line.strip()
    result.append(line)

## NOTE ON REF AND ALT
REF=result[3].split()
ALT=result[4].split()
del REF[0]
del ALT[0]

for i in range(9,len(result)):
    New=result[i].split()
    CUN_1=list()
    CUN_2=list()
    
    for j in New[1:]:
        New_SPLIT=j.split("|")
        CUN_1.append(New_SPLIT[0])
        CUN_2.append(New_SPLIT[1])
    
    TCUN_1=list()
    for k in range(0,len(CUN_1)):
        if CUN_1[k]=="1":
            TCUN_1.append(ALT[k])
        else:
            TCUN_1.append(REF[k])
    TCUN_2=list()     
    for k in range(0,len(CUN_2)):
        if CUN_2[k]=="1":
            TCUN_2.append(ALT[k])
        else:
            TCUN_2.append(REF[k])
    print(str(">"+New[0]+"_1"+"\n"+"".join(TCUN_1)))
    print(str(">"+New[0]+"_2"+"\n"+"".join(TCUN_2)))

