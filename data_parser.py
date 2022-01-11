'''
Macro de mise en forme du tableau csv de la Drees

Génére un fichier texte contenant seulement les valeurs numériques
nécessaires, et dans un format qui permet de génerer facilement un fichier ROOT

    Categories (45 catégories)
     
    --> Statut vaccinal (vac)
     
     0 = Complet de 6 mois et plus - avec rappel
     1 = Complet de 6 mois et plus - sans rappel
     2 = Complet de moins de 3 mois - avec rappel
     3 = Complet de moins de 3 mois - sans rappel
     4 = Complet entre 3 mois et 6 mois - avec rappel
     5 = Complet entre 3 mois et 6 mois - sans rappel
     6 = Non-vaccinés
     7 = Primo dose efficace  // Pas consideré
     8 = Primo dose récente   // Pas consideré
     
     --> Tranche d'age
     
     0 = [40,59]
     1 = [80+]
     2 = [0,19]               // Pas consideré
     3 = [20,39]
     4 = [60,79]
     
     Peuvent changer selon le fichier csv à mettre au point

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
vac_status=['Complet de 6 mois et plus - avec rappel','Complet de 6 mois et plus - sans rappel','Complet de moins de 3 mois - avec rappel','Complet de moins de 3 mois - sans rappel','Complet entre 3 mois et 6 mois - avec rappel','Complet entre 3 mois et 6 mois - sans rappel','Non-vaccinés','Primo dose efficace','Primo dose récente']
age_range=['[40,59]','[80,+]','[0,19]','[20,39]','[60,79]']


dates=[]
for line in data:

    info=line.split(';')
    if info[0]=='date': # Skip the first line
        continue
      
    # La Drees ne donne pas toujours la date avec le meme format...
    
    thedate=''
    if '/' not in info[0]:
        thedate=datetime.strptime(info[0], "%Y-%m-%d")
    else:
        thedate=datetime.strptime(info[0], "%d/%m/%Y")
        
    mydate=time.mktime(thedate.timetuple())
    
    if mydate not in dates:
        bunch=[]
        bunch.append(mydate)
        dates.append(mydate)
        for i in range(5):
            results=np.zeros(70)
            bunch.append(results)
        data_by_date.append(bunch)

    idx=dates.index(mydate)

    type=info[1]
    age=info[2]
    if type not in vac_status:
        vac_status.append(type)
        print('Uh, new vac status, strange'+type)
    if age not in age_range:
        age_range.append(age)
        print('Uh, new age range, strange'+age)
        
    idx_vax=vac_status.index(type)
    idx_age=age_range.index(age)
 
    # On garde seulement les entrées en soins critiques avec ou sans
    # PCR+
    # On garde aussi la population correspondante, pour normaliser
    
    HO=-1
    HOPCR=-1
    SC=-1.
    SCPCR=-1.
    DC=-1.
    DCPCR=-1.
    
    for i in range(7):
        data_by_date[idx][idx_age+1][7*idx_vax+i]=-1
    
    # Normalisation par millions d'habitants
    # On supprime dès le départ les points avec mons de 10000 personnes
    if float(info[13])>=0:
        HO=float(info[7])/(float(info[13])/1000000.)
        HOPCR=float(info[8])/(float(info[13])/1000000.)
        SC=float(info[9])/(float(info[13])/1000000.)
        SCPCR=float(info[10])/(float(info[13])/1000000.)
        DC=float(info[11])/(float(info[13])/1000000.)
        DCPCR=float(info[12])/(float(info[13])/1000000.)
        
    data_by_date[idx][idx_age+1][7*idx_vax]=HO
    data_by_date[idx][idx_age+1][7*idx_vax+1]=HOPCR
    data_by_date[idx][idx_age+1][7*idx_vax+2]=SC
    data_by_date[idx][idx_age+1][7*idx_vax+3]=SCPCR
    data_by_date[idx][idx_age+1][7*idx_vax+4]=DC
    data_by_date[idx][idx_age+1][7*idx_vax+5]=DCPCR
    data_by_date[idx][idx_age+1][7*idx_vax+6]=float(info[13])/1000000.

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


