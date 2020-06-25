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

    Author : virginie.Guemas@meteo.fr - 2020
    """
    lstfreq=('Hourly','Inthourly','Daily','Intdaily')
    if freq not in lstfreq: 
      print('Argument freq should be in',lstfreq)
      return

    lstfill=(999.,9999.,99999.)

    filename={'Hourly':'prof_file_all6_ed_hd.txt','Inthourly':'main_file6_hd.txt','Daily':'prof_file_davg_all6_ed.txt','Intdaily':'main_file_davg6_n4_hd.txt'}
    splitarg={'Hourly':'\t','Inthourly':'\t','Daily':None,'Intdaily':'\t'}

    f=open(rootpath+'SHEBA/Tower/'+freq+'/'+filename[freq],'rU')
    lines=f.readlines()
    f.close()
    
    for iline in range(len(lines)):
      lines[iline]=lines[iline].strip('\n') 
      lines[iline]=lines[iline].split(splitarg[freq])
    
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
    """
    This function loads the SHEBA PAM station data. It takes one argument :
    - freq = '1hour' / '5min'

    Warning : The 5min files contain a character string 'station' that cdms can
    not read through cdscan but the station dimension still needs to be read
    even thought the station variable can not. The workaround changes the station
    variable name to stationame at first reading and create an empty STATIONRENAMED
    file. This existence of this file is tested to know whether the renaming has
    to be done or not.

    Author : virginie.Guemas@meteo.fr - 2020  
    """

    lstpamtab=[]      # 4 output dicts, 1 per PAM station
    filebase={'1hour':'isff','5min':'sheba'}
    dirname=rootpath+'SHEBA/Mesonet_PAMIII/'+freq+'/'
    rootname=dirname+filebase[freq]

    # Work around to skip the reading of the station variable which can not be read by cdscan because it is a character string whereas cdscan needs to read the station dimension.
    if freq=='5min':
      os.system('if [ ! -e '+dirname+'STATIONRENAMED ] ; then for file in `ls '+dirname+'*nc` ; do ncrename -v station,stationname $file ; done ; fi')
      os.system('touch '+dirname+'STATIONRENAMED')

    # Reading of all netcdf files at once
    os.system('cdscan --exclude var,stationname -x cdsample.xml '+rootname+'*.nc')
    f=cdms.open('cdsample.xml')
    lstvars=f.listvariables()
    lstvars.remove('base_time')
    # Which variables to include in each station table
    for station in xrange(4): # 4 stations
      table={}
      for var in lstvars:
        table[var]=f[var][:,station]
        #f[var] is already a masked cdms variable
        if 'short_name' in table[var].attributes : 
          table[var].name=table[var].short_name
      lstpamtab.append(table)
    f.close()
    os.remove('cdsample.xml')

    return lstpamtab      
################################################################################

