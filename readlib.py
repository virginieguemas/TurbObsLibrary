# This modules contain an overarching function reading all observational data 
# from this database together with its subroutines : one for each dataset.
#
# Author : Virginie Guemas - 2020
###############################################################################
import numpy as np
import cdms2 as cdms
import MV2 as MV
import os

rootpath='/home/guemas/Obs/'
###############################################################################
def main(campaign='sheba',record='tower',freq='Hourly') : 

    maintable={}
    if campaign =='sheba':
        maintable.update(shebatower(freq))
    return maintable 
###############################################################################
def shebatower(freq='Hourly'):
    """
    This function loads the SHEBA tower data. It takes one argument :
    - freq = 'Daily' / 'Hourly' / 'Inthourly' / 'Intdaily'. 
    """
    lstfreq=('Hourly','Inthourly','Daily','Intdaily')
    if freq not in lstfreq: 
      print('Argument freq should be in',lstfreq)
      return

    lstfill=(999.,9999.,99999.)

    filename={'Hourly':'prof_file_all6_ed_hd.txt','Inthourly':'main_file6_hd.txt','Daily':'prof_file_davg_all6_ed.txt','Intdaily':'main_file_davg6_n4_hd.txt'}
    f=open(rootpath+'SHEBA/Tower/'+freq+'/'+filename[freq],'r')
    lines=f.readlines()
    f.close()
    
    if freq=='Daily':
        lines=lines[0].split('\r')
        lines.remove('')
        for iline in range(len(lines)):
            lines[iline]=lines[iline].split()
    else :
        for iline in range(len(lines)):
            lines[iline]=lines[iline].strip('\n') 
            lines[iline]=lines[iline].split('\t')
    
    table={}
    startval={'Hourly':2,'Inthourly':1,'Intdaily':1,'Daily':0}
    for ifld in range(len(lines[0])):
        values=[]
        for iline in range(startval[freq],len(lines)):
            values.append(float(lines[iline][ifld]))
        values=np.array(values)
        for fill in lstfill:
            values=np.where(values==fill,np.nan,values)
        if freq == 'Hourly': 
            table[lines[0][ifld]]={'values':values,'unit':lines[1][ifld]}
        elif freq == 'Intdaily' or freq == 'Inthourly': 
            table[lines[0][ifld]]={'values':values}
        elif freq == 'Daily':
            table[ifld]={'values':values}

    return table
################################################################################
def shebapam(freq='1hour'):

    lstpamtab=[]      # 4 output dicts, 1 per PAM station
    rootname=rootpath+'SHEBA/Mesonet_PAMIII/1hour/isff'
    os.system('cdscan -x cdsample.xml '+rootname+'*.nc')
    f=cdms.open('cdsample.xml')
    lstvars=f.listvariables()
    lstvars.remove('base_time')
    # Which variables to include in each station table
    for station in xrange(4): # 4 stations
      table={}
      for var in lstvars:
        table[var]=f[var][:,station]
        #f[var] is already a masked cdms variable
        table[var].name=table[var].short_name
      lstpamtab.append(table)
    f.close()
    os.remove('cdsample.xml')

    return lstpamtab      
################################################################################

