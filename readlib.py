# This modules contain an overarching function reading all observational data 
# from this database together with its subroutines : one for each dataset.
#
# Author : Virginie Guemas - 2020
###############################################################################
import numpy as np
import cdms2 as cdms
import MV2 as MV
import os
import sys

rootpath='/home/guemas/Obs/'
###############################################################################
def main(campaigns=['sheba'],sites=['tower'],freq='Hourly') : 
    """
    This function loads any observational data from this database.
    It takes three arguments :
    - campaigns = a list of campaign names. Default : ['sheba']
    - sites = a list of sheba sites amongst ['tower', 'Atlanta','Cleveland-Seattle-Maui','Baltimore','Florida']. Default : ['tower']   
    - freq = '5min' / 'Hourly' / 'Inthourly' / 'Daily' / 'Intdaily'. Default : 'Hourly'. Only 'Hourly' is available for all sites.

    Author : virginie.guemas@meteo.fr - 2020
    """
    if not isinstance(campaigns,list):
      sys.exit('Argument campaigns should be a list')

    if not isinstance(sites,list):
      sys.exit('Argument sites should be a list')

    if freq == 'Hourly' :
      freqpam = '1hour'
    else:
      freqpam = freq

    lstfiles=[]
    for campaign in campaigns:
      if campaign =='sheba':
        if 'tower' in sites:
          if sites != ['tower']:
            index = sites.index('tower')
            sites.remove('tower')
            tmp = shebapam(freqpam,sites)
            tmp.insert(index,shebatower(freq))
            lstfiles.extend(tmp)
          else:
            lstfiles.extend(shebatower(freq))
        else:
          lstfiles.extend(shebapam(freqpam,sites))
      else:
        sys.exit('Error : unknown campaign in the campaigns list')

    return lstfiles 
###############################################################################
def shebatower(freq='Hourly'):
    """
    This function loads the SHEBA tower data. It takes one argument :
    - freq = 'Daily' / 'Hourly' / 'Inthourly' / 'Intdaily'. Default : 'Hourly' 

    Author : virginie.Guemas@meteo.fr - 2020
    """

    lstfreq=('Hourly','Inthourly','Daily','Intdaily')
    if freq not in lstfreq: 
      sys.exit(('Argument freq should be in ',lstfreq))

    lstfill=(999.,9999.,99999.)

    filename={'Hourly':'prof_file_all6_ed_hd.txt','Inthourly':'main_file6_hd.txt','Daily':'prof_file_davg_all6_ed.txt','Intdaily':'main_file_davg6_n4_hd.txt'}
    splitarg={'Hourly':'\t','Inthourly':'\t','Daily':None,'Intdaily':'\t'}
    # Character between the columns of the ASCII file

    f=open(rootpath+'SHEBA/Tower/'+freq+'/'+filename[freq],'rU')
    lines=f.readlines()
    f.close()
    
    for iline in range(len(lines)): 
    # Toward a list of lines which are lists of column values
      lines[iline]=lines[iline].strip('\n') 
      lines[iline]=lines[iline].split(splitarg[freq])
    
    table={} # 1 dict containing all cdms variables referenced through their file ids
    startval={'Hourly':2,'Inthourly':1,'Intdaily':1,'Daily':0}
    # Some files do not have the units and one not even the variable names 
    for ifld in range(len(lines[0])):
      var=lines[0][ifld] # This will be the variable name if existing in the file
      values=[]
      for iline in range(startval[freq],len(lines)):
        values.append(float(lines[iline][ifld]))
      values=MV.array(values)
      # values along a column organised into an array
      if ifld == 0 : 
        values.units='days since 1997-01-01 00:00:00'
        values.id='time'
        time=cdms.createAxis(values)
        time.id='time'
      # JD is the time axis to be provided to all other cdms variables
      else :
        for fill in lstfill:
          values=MV.masked_where(values==fill, values)
        # the array becomes a cdms masked variable 
        if freq == 'Hourly': 
          values.units=lines[1][ifld] # Units only available in the Hourly file
        values.setAxisList((time,)) # Set the time axis
        if freq == 'Daily':
          table[ifld]=values # Not even variable name in the Daily file
        else:
          values.id=lines[0][ifld]
          values.name=lines[0][ifld]
          table[lines[0][ifld]]=values

    return table
################################################################################
def shebapam(freq='1hour',sites=['Atlanta','Cleveland-Seattle-Maui','Baltimore','Florida']):
    """
    This function loads the SHEBA PAM station data. It takes two arguments :
    - freq = '1hour' / '5min'. Default : '1hour'
    - sites = a list of sites amongst ['Atlanta','Cleveland-Seattle-Maui','Baltimore','Florida']. Default : ['Atlanta','Cleveland-Seattle-Maui','Baltimore','Florida']

    Warning : The 5min files contain a character string 'station' that cdms can
    not read through cdscan but the station dimension still needs to be read
    even thought the station variable can not. The workaround changes the station
    variable name to stationame at first reading and create an empty STATIONRENAMED
    file. This existence of this file is tested to know whether the renaming has
    to be done or not.

    Author : virginie.guemas@meteo.fr - 2020  
    """

    lstfreq=('1hour','5min')
    if freq not in lstfreq: 
      sys.exit(('Argument freq should be in',lstfreq))
    
    if not isinstance(sites,list):
      sys.exit('Argument sites should be a list')

    lstpamtab=[]      # 1 output dict per PAM station listed in sites
    stationames=('Atlanta','Cleveland-Seattle-Maui','Baltimore','Florida')
    columns=[]
    for site in sites:
      if site not in stationames:
        sys.exit('Argument sites should be a list of stations from',stationames)
      columns.append(stationames.index(site)) # Which columns to read

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
    for station in columns: # 4 stations
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

