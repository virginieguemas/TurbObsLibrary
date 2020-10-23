# This modules contain an overarching function reading all observational data 
# from this database together with its subroutines : one for each dataset.
#
# Author : Virginie Guemas - 2020
###############################################################################
import numpy as np
import xarray as xr
import datetime
import os
import sys
from glob import glob

rootpath='/home/guemas/Obs/'
###############################################################################
def main(campaigns=['sheba'],sites=['tower'],freq='Hourly',flights=['FAAM','MASIN']) : 
    """
    This function loads any observational data from this database.
    It takes three arguments :
    - campaigns = a list of campaign names amongst ['sheba','accacia','acse','ascos']. Default : ['sheba']
    - sites = a list of sheba sites amongst ['tower', 'Atlanta','Cleveland-Seattle-Maui','Baltimore','Florida']. Default : ['tower']   
    - freq = a list of sheba output frequency amongst '5min' / 'Hourly' / 'Inthourly' / 'Daily' / 'Intdaily'. Default : 'Hourly'. 'Hourly' and '5min' are available for the PAM stations. 'Hourly', 'Inthourly', 'Daily' and 'Intdaily' are available for the tower. 
    - flights = a list of accacia flight names amongst [ 'FAAM', 'MASIN' ]. Default : flights=['FAAM','MASIN'] 
    
    It returns a list of Xarray Datasets, one per sheba site or per ACCACIA flight or for the ACSE campaign.

    Author : virginie.guemas@meteo.fr - 2020
    """
    if not isinstance(campaigns,list):
      sys.exit('Argument campaigns should be a list')

    if not isinstance(sites,list):
      sys.exit('Argument sites should be a list')

    if not isinstance(flights,list):
      sys.exit('Argument flights should be a list')

    if freq == 'Hourly' :
      freqpam = '1hour'
    else:
      freqpam = freq

    lstds=[]
    for campaign in campaigns:
      if campaign =='sheba':
        if 'tower' in sites:
          if sites != ['tower']:
            index = sites.index('tower')
            sites.remove('tower')
            tmp = shebapam(freqpam,sites)
            tmp.insert(index,shebatower(freq))
            lstds.extend(tmp)
          else:
            tmp=shebatower(freq)
            lstds.append(shebatower(freq))
        else:
          lstds.extend(shebapam(freqpam,sites))
      elif campaign == 'accacia':
        lstds.extend(accacia(flights))  
      elif campaign == 'acse':
        lstds.append(acse())
      elif campaign == 'ascos':
        lstds.append(ascos())
      else:
        sys.exit('Error : unknown campaign in the campaigns list')

    return lstds 
###############################################################################
def shebatower(freq='Hourly'):
    """
    This function loads the SHEBA tower data. It takes one argument :
    - freq = 'Daily' / 'Hourly' / 'Inthourly' / 'Intdaily'. Default : 'Hourly' 

    This function outputs an Xarray Dataset.

    Author : virginie.guemas@meteo.fr - June 2020
    Modified : Switch from cdms to xarray - Virginie Guemas - September 2020
    """

    sys.path.append(rootpath+'SHEBA/Tower/')
    import sheba_info
    # Where to find a subsidiary file where I stored the information I found on the variables

    lstfreq=('Hourly','Inthourly','Daily','Intdaily')
    if freq not in lstfreq: 
      sys.exit(('Argument freq should be in ',lstfreq))

    lstfill=(9.99,999.,9999.,99999.) # Missing values

    filename={'Hourly':'prof_file_all6_ed_hd.txt','Inthourly':'main_file6_hd.txt','Daily':'prof_file_davg_all6_ed.txt','Intdaily':'main_file_davg6_n4_hd.txt'}
    splitarg={'Hourly':'\t','Inthourly':'\t','Daily':None,'Intdaily':'\t'}
    # Character between the columns of the ASCII file
    interp={'Hourly':0,'Inthourly':1,'Daily':0,'Intdaily':1}
    # Flag because different variables are stored in 'hl' in interpolated files or not

    f=open(rootpath+'SHEBA/Tower/'+freq+'/'+filename[freq],'rU')
    lines=f.readlines()
    f.close()
    
    for iline in range(len(lines)): 
    # Toward a list of lines which are lists of column values
      lines[iline]=lines[iline].strip('\n') 
      lines[iline]=lines[iline].split(splitarg[freq])
    
    ds=xr.Dataset() 
    # One dataset containing all DataArrays variables available 
    for ifld in range(len(lines[2])): 
    # First and second lines might have additional variables for which there are no values
      values=[]
      for iline in range(2,len(lines)): # First line is the variable name, second the unit
        values.append(float(lines[iline][ifld]))
        # values along a column organised into an array

      for fill in lstfill:
        values=np.where(np.array(values)==fill, np.nan, np.array(values))
        # mask values in lstfill  

      array=xr.DataArray(values,dims=('time'),attrs={'long_name':sheba_info.shebatower_names(lines[0][ifld],interp[freq]),'units':lines[1][ifld]})
      # the list of values becomes an Xarray DataArray
      name=lines[0][ifld]
      name=name.replace('*','star') 
      # '*' in the variable name is not an accepted charater to call it back
      ds[name]=array
      # the Xarray DataArray is included into the sheba tower dataset

      if ifld == 0 :
        timecoord=[]
        for ii in range(len(values)):
          timecoord.append(datetime.datetime(1996,12,31,0,0)+datetime.timedelta(seconds=values[ii]*86400.))
        ds['time']=timecoord
        # JD (first column) is the time axis to be provided to the whole dataset
  
    return ds
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
    # Concatenation
    
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
    # Where to find a subsidiary file where I stored the information I was provided on the variables

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
      # One dataset containing all variables available for each flight
      for ifld in range(len(lines[0])):
        values=[]
        for iline in range(1,len(lines)): # First line is variable name
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
          # Creation of time axis from the calday and meantime variables
      elif flight == 'MASIN':
        time_bnds=[]
        for ii in range(len(values)):
          timecoord.append(datetime.datetime(2012,12,31)+datetime.timedelta(days=ds.dayofyear.values[ii],seconds=ds.meantime.values[ii]))
          # Creation of time axis from the dayofyear and meantime variables
          time_bnds.append(datetime.datetime(2012,12,31)+datetime.timedelta(days=ds.dayofyear.values[ii],seconds=ds.starttime.values[ii]))
          time_bnds.append(datetime.datetime(2012,12,31)+datetime.timedelta(days=ds.dayofyear.values[ii],seconds=ds.endtime.values[ii]))
          # Creation of time bounds
        ds['time_bnds']=xr.DataArray(np.reshape(time_bnds,(int(len(time_bnds)/2),2)),dims=('time','bnds'))
      ds=ds.assign_coords(time=timecoord)
      # Inclusion of time axis in the flight dataset

      lstaccdat.append(ds)

    return lstaccdat
################################################################################
def acse():
    """
    This function reads the ACSE campaign data and outputs an Xarray dataset. 

    Author : virginie.guemas@meteo.fr - October 2020
    """
    sys.path.append(rootpath+'ACSE/Version2/')
    import acse_info

    ds=xr.open_dataset(rootpath+'/ACSE/Version2/ACSE_CANDIFLOS_fluxes_Oden_20140710_v5.0.nc',drop_variables='doy')
    # I drop doy because there is a bug in its conversion to a time series and we do not need it because we have a correct time index.

    for var in ds.keys():
      if acse_info.acse_names(var) != '':
        ds[var].attrs = {'long_name':acse_info.acse_names(var),'units':ds[var].units}

    return ds
################################################################################
def ascos():
    """
    This function reads the ASCOS campaign data and outputs an Xarray dataset. 

    Author : virginie.guemas@meteo.fr - October 2020
    """
    
    ds=xr.open_dataset(rootpath+'/ASCOS/Metalley/ascos_Metalley.nc',decode_times=False)

    timecoord=[]
    # Define the time axis from year, month, day, hour, minute, second
    months=np.where(ds.day<10,9,8)
    for ii in range(len(ds['time'])):
      timecoord.append(datetime.datetime(year=2008,month=months[ii],day=ds.day[ii],hour=ds.hour[ii],minute=ds.minute[ii],second=ds.second[ii]))
    ds=ds.assign_coords(time=timecoord)

    # Define height level coordinates to define variables depending on (time, height) instead of var_height1, var_height2 ...
    height_axis1=[0.94, 4.04, 5.21, 8.19, 15.4, 30.6]
    ds=ds.assign_coords(height_axis1=height_axis1)
    height_axis2=[0.2,1.02,1.79,5.32,8.36]    
    ds=ds.assign_coords(height_axis2=height_axis2)
    height_axis3=[3.19,14.92]
    ds=ds.assign_coords(height_axis3=height_axis3)
    height_axis4=[0.05,0.15,0.4,1.]
    ds=ds.assign_coords(height_axis4=height_axis4)
    height_axis5=np.arange(11)+1
    ds=ds.assign_coords(height_axis5=height_axis5)

    # Reorganize var_height1, var_height2 ... in var
    lstvars1 = ('wspd','ts','g','ugeo','vgeo','wdirgeo','u','v','wdir','w','tke','tsts','uu','vv','ww','wts','uv','uw','vw','flag_csat','flag_grill','flag_gill','flag_licor','tc','q','tctc','wtc','wCO2','wq')
    lstvars2 = ('TC_Air',)
    lstvars3 = ('RH_Air','T_Air')
    lstvars4 = ('TC_Ice',)
    lstvars5 = ('TC_Surf',)
    for var in lstvars1+lstvars2+lstvars3+lstvars4+lstvars5:
      if var in lstvars1:
        height = height_axis1
        nameh = 'height_axis1'
      elif var in lstvars2:
        height = height_axis2
        nameh = 'height_axis2'
      elif var in lstvars3:
        height = height_axis3
        nameh = 'height_axis3'
      elif var in lstvars4:
        height = height_axis4
        nameh = 'height_axis4'
      elif var in lstvars5:
        height = height_axis5
        nameh = 'height_axis5'
      var0 = np.empty((len(ds['time']),len(height)))*np.nan
      for hh in range(len(height)):
        if var in lstvars5:
          suffix = str(int(height[hh]))
        else:
          suffix = str(int(height[hh]*100))+'_cm'
        name = var+'_'+suffix
        if name in ds.keys(): # Some levels missing for some variables
          var0[:,hh] = ds[name].values
          attrs = ds[name].attrs
      xr_var = xr.DataArray(var0, dims=['time',nameh], attrs=attrs)
      ds[var] = xr_var
    # Removing var_height1, var_height2 ...
      for hh in range(len(height)):
        if var in lstvars5:
          suffix = str(int(height[hh]))
        else:
          suffix = str(int(height[hh]*100))+'_cm'
        name = var+'_'+suffix
        if name in ds.keys():
          ds=ds.drop(name)

    ds['height_axis1'].attrs={'units':'m'}
    ds['height_axis2'].attrs={'units':'m'}
    ds['height_axis3'].attrs={'units':'m'}
    ds['height_axis4'].attrs={'units':'m'}

    return ds
################################################################################
