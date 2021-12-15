'''
Macro de mise en forme du tableau csv de la Drees

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
  print('Ratio.py -i <inputfile>')

inputfile=''

for opt, arg in opts:
  if opt == '-h':
    print('Ratio.py -i <inputfile>')
    sys.exit()
  elif opt in ("-i", "--input"):
    inputfile=arg

runtime=datetime.now().isoformat()
outputfile='Ratio_'+runtime+'.txt'
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
    if info[0]=='date': # Skip the first line
        continue
        
    thedate=datetime.strptime(info[0], "%Y-%m-%d")
    mydate=time.mktime(thedate.timetuple())
    
    if mydate not in dates:
        bunch=[]
        bunch.append(mydate)
        dates.append(mydate)
        for i in range(5):
            results=np.zeros(50)
            bunch.append(results)
        data_by_date.append(bunch)

    idx=dates.index(mydate)

    type=info[1]
    age=info[2]
    if type not in vac_status:
        vac_status.append(type)
    if age not in age_range:
        age_range.append(age)
        
    idx_vax=vac_status.index(type)
    idx_age=age_range.index(age)
 
    # On garde seulement les entrées en soins critiques avec ou sans
    # PCR+
    # On garde aussi la population correspondante, pour normaliser
    
    HO=-1
    HOPCR=-1
    SC=-1.
    SCPCR=-1.
    
    for i in range(5):
        data_by_date[idx][idx_age+1][5*idx_vax+i]=-1
    
    # Normalisation par millions d'habitants
    # On supprime dès le départ les points avec mons de 10000 personnes
    if float(info[9])>=10000:
        HO=float(info[5])/(float(info[9])/1000000.)
        HOPCR=float(info[6])/(float(info[9])/1000000.)
        SC=float(info[7])/(float(info[9])/1000000.)
        SCPCR=float(info[8])/(float(info[9])/1000000.)

    data_by_date[idx][idx_age+1][5*idx_vax]=HO
    data_by_date[idx][idx_age+1][5*idx_vax+1]=HOPCR
    data_by_date[idx][idx_age+1][5*idx_vax+2]=SC
    data_by_date[idx][idx_age+1][5*idx_vax+3]=SCPCR
    data_by_date[idx][idx_age+1][5*idx_vax+4]=float(info[9])/1000000.

# A la fin on remplit le fichier texte

op.write('\n')
for data in data_by_date:

    op.write(str(data[0])+'\n')
    for i in range(5):
        for idx in range(len(data[i+1])):
            op.write(str(data[i+1][idx])+'\n')

print("Catégories de vaccination détectées:")
print(vac_status)
print("Catégories d'âge détectées:")
print(age_range)
op.close()
sys.exit()


