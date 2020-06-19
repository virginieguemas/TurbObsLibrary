# This modules contain a overarching function reading all observational data 
# from this database together with its subroutines : one for each dataset.
#
# Author : Virginie Guemas - 2020
###############################################################################
import numpy as np

rootpath='/home/guemas/Obs/'
###############################################################################
def main(campaign='sheba',record='tower') : 

    maintable={}
    if campaign =='sheba':
        maintable.update(shebatower())
    return maintable 
###############################################################################
def shebatower(freq='Hourly'):
    
    lstfill=(999.,9999.,99999.)
    f=open(rootpath+'SHEBA/Tower/'+freq+'/prof_file_all6_ed_hd.txt','r')
    lines=f.readlines()
    f.close()
    for iline in range(len(lines)):
        lines[iline]=lines[iline].strip('\n') 
        lines[iline]=lines[iline].split('\t')
    table={}
    for ifld in range(len(lines[0])):
        values=[]
        for iline in range(2,len(lines)):
            values.append(float(lines[iline][ifld]))
        values=np.array(values)
        for fill in lstfill:
            values=np.where(values==fill,np.nan,values)
        table[lines[0][ifld]]={'values':np.array(values),'unit':lines[1][ifld]}
    return table
