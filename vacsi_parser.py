'''
Macro de mise en forme du tableau csv de VACSI national par age

Génére un fichier texte contenant seulement les valeurs numériques
nécessaires, et dans un format qui permet de génerer facilement un fichier ROOT
'''

import sys, getopt
import os
import six
import math
from array import array
from datetime import datetime
import time
import numpy as np

try:
  opts, args = getopt.getopt(sys.argv[1:],"hi",["input="])
except getopt.GetoptError:
  print('vacsi_parser.py -i <inputfile>')

inputfile=''

for opt, arg in opts:
  if opt == '-h':
    print('vacsi_parser.py -i <inputfile>')
    sys.exit()
  elif opt in ("-i", "--input"):
    inputfile=arg

runtime=datetime.now().isoformat()
outputfile='vacsi_parsed_'+runtime+'.txt'
op=open(outputfile, 'w')

inFile = open(inputfile,'r')

data = inFile.readlines()

treated_data=[]
data_by_date=[]
vac_status=[]
age_range=[]
dates=[]
for line in data:

    info=line.split(';')
    if info[0]=='fra': # Skip the first line
        continue
        
    thedate=datetime.strptime(info[2], "%d/%m/%Y")
    mydate=time.mktime(thedate.timetuple())
    
    if mydate not in dates:
        bunch=[]
        bunch.append(mydate)
        dates.append(mydate)
        for i in range(15):
            results=np.zeros(4)
            bunch.append(results)
        data_by_date.append(bunch)

    idx=dates.index(mydate)

    age=info[1]
    if age not in age_range:
        age_range.append(age)
    #print(age_range)
    idx_age=age_range.index(age)
    
    #print(idx,idx_age)
    
    for i in range(4):
        data_by_date[idx][idx_age+1][i]=-1
    
    ndose1=float(info[6])
    ncomplet=float(info[7])
    nrappel=float(info[8])
        
    data_by_date[idx][idx_age+1][0]=age
    data_by_date[idx][idx_age+1][1]=ndose1
    data_by_date[idx][idx_age+1][2]=ncomplet
    data_by_date[idx][idx_age+1][3]=nrappel


# A la fin on remplit le fichier texte

op.write('\n')
for data in data_by_date:

    op.write(str(data[0])+'\n')
    for i in range(15):
        for idx in range(len(data[i+1])):
            op.write(str(data[i+1][idx])+'\n')

print("Catégories d'âge détectées:")
print(age_range)
op.close()
sys.exit()


