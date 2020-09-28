# This modules contain an overarching function reading all observational data 
# from this database together with its subroutines : one for each dataset.
#
# Author : Virginie Guemas - 2020
###############################################################################
import numpy as np
#mport cdms2 as cdms
#import MV2 as MV
import xarray as xr
import datetime
import os
import sys
from glob import glob

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

    Author : virginie.Guemas@meteo.fr - June 2020
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
    This function outputs a list of Xarray Datasets, one Dataset per site.

    Author : virginie.guemas@meteo.fr - June 2020  
    Modified : Switch from cdms to xarray - Virginie Guemas - September 2020
    """

    lstfreq=('1hour','5min')
    if freq not in lstfreq: 
      sys.exit(('Argument freq should be in',lstfreq))
    
    if not isinstance(sites,list):
      sys.exit('Argument sites should be a list')

    stationames=('Atlanta','Cleveland-Seattle-Maui','Baltimore','Florida')
    for site in sites:
      if site not in stationames:
        sys.exit('Argument sites should be a list of stations from',stationames)

    filebase={'1hour':'isff','5min':'sheba'}
    dirname=rootpath+'SHEBA/Mesonet_PAMIII/'+freq+'/'
    rootname=dirname+filebase[freq]

    # Reading of all netcdf files at once
    #f=xr.open_mfdataset(rootname+'*.nc',combine='by_coords',data_vars='minimal',coords='minimal',compat='override')
    # The line above does not work properly for freq = 1hour because the Rlwdiff_in_ARM and Rlwdiff_out_ARM are present in all netcdfs but the last. The last netcdf is therefore not properly read and only nans appear. The same issue appears for freq='5min' with a few variables absent in several files. The next few lines serve as a workaround.

    paths=sorted(glob(rootname+'*.nc'))
    # list of paths of files to be loaded
    datasets = [xr.open_dataset(p) for p in paths]
    # 1 dataset per file to be all concatenated together
    if freq == '1hour':
      datasets[12]['Rlwdiff_in_ARM']=xr.DataArray(np.empty((114,4))*np.nan,dims=('time','station'))
      datasets[12]['Rlwdiff_out_ARM']=xr.DataArray(np.empty((114,4))*np.nan,dims=('time','station'))
      # Rlwdiff_in_ARM & Rlwdiff_out_ARM need to be included in all files for the concatenation to work properly. Inclusion with nan only in the files where they are missing.
    elif freq == '5min':
      for ii in range(len(datasets)) :
        if ii == paths.index(rootname+'.980207.nc'):
          lentime=287
        elif ii == paths.index(rootname+'.981005.nc'):
          lentime=209
        else:
          lentime=288
        for var in ('vdcmax_batt','vdcmin_batt','lat_sec','rainr_eti','raina_eti','ulev','vlev','Rlw_in','Rlw_out','lag_sec'):
          if var not in datasets[ii].keys():
            datasets[ii][var]=xr.DataArray(np.empty((lentime,4))*np.nan,dims=('time','station'))
        # Same for this list of variables

    f=xr.concat(datasets, dim='time')
    
    f.attrs={} # Remove the history which otherwise would be copied to each dataset and does not make sense.
    lstpamdat=[]      # 1 output dataset per PAM station listed in sites
    # Split the f dataset into one per PAM station
    for site in sites: # 4 stations 
      lstpamdat.append(f.isel(station=stationames.index(site)))

    return lstpamdat      
################################################################################
def accacia(flights=['FAAM','MASIN']):
    """
    This function loads the ACCACIA flights data. It takes one argument :
    - flights = a list of flights amongst [ 'FAAM', 'MASIN' ]. Default : flights=['FAAM','MASIN'] 
    This function outputs a list of Xarray Datasets, one Dataset per flight.

    Author : virginie.guemas@meteo.fr - July 2020
    Modified : Switch from cdms to xarray - Virginie Guemas - September 2020
    """
    sys.path.append(rootpath+'ACCACIA/')
    import accacia_info

    if not isinstance(flights,list):
      sys.exit('Argument flights should be a list')

    flightnames=('FAAM','MASIN')
    lstaccdat=[]  # 1 output dataset per flight listed in flights
    for flight in flights:
      if flight not in flightnames:
        sys.exit('Argument flights should be a list of flights from',flightnames)
      f=open(rootpath+'ACCACIA/ACCACIA_'+flight+'_flights_database_obs_lr_Virginie.txt','rU')
      lines=f.readlines()
      f.close()
      for iline in range(len(lines)): 
      # Toward a list of lines which are lists of column values
        lines[iline]=lines[iline].strip('\n') 
        lines[iline]=lines[iline].split()
    
      ds=xr.Dataset() 
      # 1 Dataset containing all DataArrays available for each flight
      for ifld in range(len(lines[0])):
        values=[]
        for iline in range(1,len(lines)):
          values.append(float(lines[iline][ifld]))
          # values along a column organised into an array
        array=xr.DataArray(values,dims=('time'),attrs={'long_name':accacia_info.accacia_names(lines[0][ifld]),'units':accacia_info.accacia_units(lines[0][ifld])})
        # the list of values becomes an Xarray DataArray
        ds[lines[0][ifld]]=array
        # the Xarray DataArray is included into the flight Dataset
        
      timecoord=[]
      if flight == 'FAAM': 
        for ii in range(len(values)):
          timecoord.append(datetime.datetime(2013,3,ds.calday[ii])+datetime.timedelta(seconds=ds.meantime.values[ii]))
      elif flight == 'MASIN':
        time_bnds=[]
        for ii in range(len(values)):
          timecoord.append(datetime.datetime(2012,12,31)+datetime.timedelta(days=ds.dayofyear.values[ii],seconds=ds.meantime.values[ii]))
          time_bnds.append(datetime.datetime(2012,12,31)+datetime.timedelta(days=ds.dayofyear.values[ii],seconds=ds.starttime.values[ii]))
          time_bnds.append(datetime.datetime(2012,12,31)+datetime.timedelta(days=ds.dayofyear.values[ii],seconds=ds.endtime.values[ii]))
        ds['time_bnds']=xr.DataArray(np.reshape(time_bnds,(int(len(time_bnds)/2),2)),dims=('time','bnds'))
      ds=ds.assign_coords(time=timecoord)

      lstaccdat.append(ds)

    return lstaccdat
################################################################################
