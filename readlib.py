# This modules contain an overarching function reading all observational data 
# from this database, i.e. the main function, together with its subroutines,
# one for each dataset:
# - sheba
# - ascos
# - accacia
# - acse
# - ao16
# - stable
# - nsidc
# as well as a few subsidiary functions (shebaaircraft, shebapam, shebatower, 
# shebatowergather for sheba, readtabtxt for shebatower, accacia and stable, 
# oden for acse and ao16, readstablefile and postprostable for stable)
#
# Author : Virginie Guemas - 2020
###############################################################################
import numpy as np
import xarray as xr
import datetime
import os
import sys
import glob
import xlrd
import readlib, meteolib
import inspect
from glob import glob

rootpath=inspect.getfile(readlib)[0:-18]
loadice = True
###############################################################################
def main(campaigns=['sheba'],sheba_sites=['tower'],mosaic_sites=['tower'],freq='Hourly',flights=['FAAM','MASIN'],startdate='2019-10-09', enddate='2019-10-31') : 
    """
    This function loads any observational data from this database.
    It takes seven arguments :
    - campaigns = a list of campaign names amongst ['sheba','accacia','acse','ascos','ao16','stable','mosaic']. Default : ['sheba']
    - sheba_sites = a list of sheba sites amongst ['tower', 'Atlanta','Cleveland-Seattle-Maui','Baltimore','Florida','aircraft']. Default : ['tower']   
    - mosaic_sites = a list of mosaic sites amongst ['tower','asf30','asf40','asf50']. Default : ['tower']
    - freq = a list of sheba output frequency amongst '5min' / 'Hourly' / 'Daily'. Default : 'Hourly'. 'Hourly' and '5min' are available for the PAM stations. 'Hourly' and 'Daily' are available for the tower. There is no choice for the aircrafts. 
    - flights = a list of accacia flight names amongst [ 'FAAM', 'MASIN' ]. Default : flights=['FAAM','MASIN'] 
    - startdate = YYYY-MM-DD. Only for MOSAiC. Example : '2019-10-09'
    - enddate = YYYY-MM-DD. Only for MOSAiC. Example : '2020-10-01'
    MOSAiC data starts on 2019-10-09 and end on 2020-10-01.
    
    It returns a list of Xarray Datasets, one per sheba site and/or per mosaic site and/or one for the sheba aircraft and/or one per ACCACIA flight and/or one for the ACSE campaign and/or one for the ASCOS campaign and/or one for the AO16 campaign and/or one for the STABLE campaign.

    Author : virginie.guemas@meteo.fr - 2020
    """
    if not isinstance(campaigns,list):
      sys.exit('Argument campaigns should be a list')

    if not isinstance(sheba_sites,list):
      sys.exit('Argument sheba_sites should be a list')

    if not isinstance(mosaic_sites,list):
      sys.exit('Argument mosaic_sites should be a list')

    if not isinstance(flights,list):
      sys.exit('Argument flights should be a list')

    if freq == 'Hourly' :
      freqpam = '1hour'
    else:
      freqpam = freq

    lstds=[]
    for campaign in campaigns:
      if campaign =='sheba':
        indext = -1
        indexa = -1
        if 'tower' in sheba_sites:
          indext = sheba_sites.index('tower')
          sheba_sites.remove('tower')
        if 'aircraft' in sheba_sites:   
          indexa = sheba_sites.index('aircraft')
          sheba_sites.remove('aircraft')
        if sheba_sites != []:
          tmp = shebapam(freqpam,sheba_sites)
        else:
          tmp = []
        if indexa >= 0:
          tmp.insert(indexa,shebaaircraft())
        if indext >= 0:
          tmp.insert(indext,shebatowergather(freq))
        lstds.extend(tmp)
      elif campaign == 'accacia':
        lstds.extend(accacia(flights))  
      elif campaign == 'acse':
        lstds.append(acse())
      elif campaign == 'ascos':
        lstds.append(ascos())
      elif campaign == 'ao16':
        lstds.append(ao16())
      elif campaign == 'stable':
        lstds.append(stable())
      elif campaign == 'mosaic':
        for site in mosaic_sites:
          lstds.append(mosaic(site,startdate,enddate))
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
               Gather all levels from each variable by defining level axis - 
                                            Virginie Guemas - October 2020
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
 
    # Define height level coordinates to define variables depending on (time, height) instead of var1, var2 ...
    if freq == 'Hourly':
      lstvars = ('z','ws','wd','T','q','rh','rhi','ustar','hs','ww','sgu','sgv','sgw','sgT','cT2','cu2','cv2','cw2','No','fl')
      height = np.arange(5)+1
      nameh = 'level'
      ds = ds.assign_coords(level=height)
    elif freq == 'Daily':
      lstvars = ('z','ws','wd','T','q','rh','rhi','ustar','hs','ww')
      height = np.arange(5)+1
      nameh = 'level'
      ds=ds.assign_coords(level=height)
    elif freq == 'Inthourly' or freq == 'Intdaily':
      lstvars = ('ws','wd','T','q','rhi','usb','hsb','hlb','zob','zotb','zoqb','zogb')
      height = [2.5, 10]
      nameh ='height'
      ds=ds.assign_coords(height=height)

    # Reorganize var1, var2 ... in var
    for var in lstvars:
      var0 = np.empty((len(ds['time']),len(height)))*np.nan
      for hh in range(len(height)):
        if var == 'usb' or var == 'hsb' or var == 'hlb'  :
          name = var + '_' + str(height[hh])
        else :
          name = var+str(height[hh])
        if name in ds.keys(): # Some missing levels
          var0[:,hh] = ds[name].values
          longname = ds[name].long_name.replace(' of fifth level','').replace(' at fifth level','').replace(' at level 5','').replace(' at 10m','')
          units = ds[name].units
          ds=ds.drop(name)
      xr_var = xr.DataArray(var0, dims=['time',nameh], attrs={'long_name':longname,'units':units})
      ds[var] = xr_var

    if freq == 'Inthourly' or freq == 'Intdaily' :
      ds['height'].attrs={'units':'m'}

    # Include NSIDC sea ice concentrations
    if loadice:
      sic = nsidc(lat=ds.lat,lon=ds.lon)
      ds = xr.Dataset.merge(ds, sic)
      ds.attrs['nsidc_g02202v3'] = sic.nsidc_g02202v3
      sic = nsidc(lat=ds.lat,lon=ds.lon,dataset='g02202v4')
      ds = xr.Dataset.merge(ds, sic)
      ds.attrs['nsidc_g02202v4'] = sic.nsidc_g02202v4

    return ds
################################################################################
def shebatowergather(freq='Hourly'):
    """
    This function gathers in the same Xarray Dataset the SHEBA tower data from 
    the 'Hourly' and 'Inthourly' files on one side or from the 'Daily' and 
    'Intdaily' files on the other side.
    - freq = 'Daily' / 'Hourly'. Default : 'Hourly' 

    This function outputs an Xarray Dataset.

    Author : virginie.guemas@meteo.fr - October 2020
    """
    
    lstfreq=('Hourly','Daily')
    if freq not in lstfreq: 
      sys.exit(('Argument freq should be in ',lstfreq))

    ds = shebatower(freq)
    ds2 = shebatower('Int'+freq.lower())

    newname = {'hl':'hlmed','ws':'wsint','wd':'wdint','T':'Tint','q':'qint','rhi':'rhint','ustar':'ustarmed','hs':'hsmed','ww':'wwbis'}
    if freq == 'Daily':
      newname['RR_org']='RR_org_int'
      newname['RR_ncr']='RR_ncr_int' # Different values in Daily and Intdaily but I do not know what those variables here.

    for name in ds2.keys():
      if name in newname:
        ds[newname[name]]=ds2[name] 
      else:
        ds[name]=ds2[name]  # Some variables are identical and will be overwritten - some variables contain more missing values in Intdaily than in Daily - Weird but fixed by this line

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
        sys.exit(('Argument sites should be a list of stations from',stationames))

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
      ds = f.isel(station=stationames.index(site))

      # Include NSIDC sea ice concentrations
      if loadice:
        sic = nsidc(lat=ds.latitude,lon=ds.longitude)
        ds = xr.Dataset.merge(ds, sic)
        ds.attrs['nsidc_g02202v3'] = sic.nsidc_g02202v3
        sic = nsidc(lat=ds.latitude,lon=ds.longitude,dataset='g02202v4')
        ds = xr.Dataset.merge(ds, sic)
        ds.attrs['nsidc_g02202v4'] = sic.nsidc_g02202v4

      # Output as a list of datasets per PAM station
      lstpamdat.append(ds)

    return lstpamdat      
################################################################################
def shebaaircraft():
    """
    This function loads the SHEBA Aircraft data (Perovich et al. 2002). It takes no argument.
    This function outputs an Xarray Dataset.

    Author : Sebastien Blein - December 2020  
    """
    #Function issues:
    # - Short names have been added to the original dataset
    # - deal with the 'N/A' values and cells with several values.
    # - deal with the Lat Long format
    # - to be optimized?
    ds=xr.Dataset()
    
    loc = (rootpath+'SHEBA/Aircraft/DonaldPerovich/JGRtables_sb.xls')
    xls_book = xlrd.open_workbook(loc, formatting_info = True)

    datecoord=[]
    datecoord2=[]
    check_first_date = True
    for isheet in np.arange(3):
        xls_sheet = xls_book.sheet_by_index(isheet)
        #
        for row in range(xls_sheet.nrows):
            for col in range(xls_sheet.ncols):
                if xls_sheet.cell_value(row, col) == 'Date':
                    rowstart, colstart = row, col
                if xls_sheet.cell_value(row, col) == 'Shortname':
                    rows_short, col_short = row, col
        if xls_sheet.cell_type(rowstart+1, colstart)==3 and xls_sheet.cell_type(rowstart, colstart+1)==1:
            time_axis = 0
        elif xls_sheet.cell_type(rowstart+1, colstart)==1 and xls_sheet.cell_type(rowstart, colstart+1)==3:
            time_axis = 1
        else:
            sys.exit('Error: check time axis in xls.sheet num. '+str(isheet))

        if time_axis == 0:
            for col in np.arange(colstart,xls_sheet.ncols):
                var_name = []
                var_val = []
                if xls_sheet.cell_value(rowstart, col) == 'Date':
                    if check_first_date == True :
                        check_first_date = False
                        date = [xls_sheet.cell_value(i, colstart) for i in np.arange(rowstart+1,xls_sheet.nrows)]
                        for d in np.arange(xls_sheet.nrows-rowstart-1):
                            datecoord.append(datetime.datetime(* xlrd.xldate.xldate_as_tuple(date[d], xls_book.datemode)))
                    if check_first_date == False :
                        date_to_check = []
                        date_to_check = [xls_sheet.cell_value(i, colstart) for i in np.arange(rowstart+1,xls_sheet.nrows)]
                        if date_to_check != date :
                            sys.exit('Error: xls sheets have different time array')
                if xls_sheet.cell_value(rowstart, col) != 'Date':
                    var_name = xls_sheet.cell_value(rows_short, col)
                    var_val = np.array([xls_sheet.cell_value(d, col) for d in np.arange(rowstart+1,xls_sheet.nrows)])
                    #if np.any([var_val[i]=='N/A' for i in np.arange(len(var_val))])==True:
                    #    var_val = np.where(np.array(var_val)=='N/A',np.nan,var_val)
                    #    var_val = np.array([np.float(var_val[i]) for i in np.arange(len(var_val))])
                    array = xr.DataArray(var_val,dims=('date'),attrs={'long_name':xls_sheet.cell_value(rowstart, col)})
                    ds[var_name]=array
        if time_axis == 1:
            for row in np.arange(rowstart,xls_sheet.nrows):
                var_name = []
                var_val = []
                if xls_sheet.cell_value(row, colstart) == 'Date':
                    if check_first_date == False :
                        date2 = [xls_sheet.cell_value(rowstart, i) for i in np.arange(colstart+1,xls_sheet.ncols)]
                        for d in np.arange(xls_sheet.ncols-colstart-1):
                            datecoord2.append(datetime.datetime(* xlrd.xldate.xldate_as_tuple(date2[d], xls_book.datemode)))
                if xls_sheet.cell_value(row, colstart) != 'Date':
                    var_name = xls_sheet.cell_value(row, col_short)
                    var_val = [xls_sheet.cell_value(row, d) for d in np.arange(colstart+1,xls_sheet.ncols)]
                    #if np.any([var_val[i]=='N/A' for i in np.arange(len(var_val))])==True:
                    #    var_val = np.where(np.array(var_val)=='N/A',np.nan,var_val)
                    #    var_val = np.array([np.float(var_val[i]) for i in np.arange(len(var_val))])
                    array = xr.DataArray(var_val,dims=('date2'),attrs={'long_name':xls_sheet.cell_value(row, colstart)})
                    ds[var_name]=array
    ds=ds.assign_coords(date=datecoord) 
    ds=ds.assign_coords(date2=datecoord2) 

    # Include NSIDC sea ice concentrations
    #sic = nsidc(lat=ds.rename({'date':'time'}).lat,lon=ds.rename({'date':'time'}).lon)
    #ds = xr.Dataset.merge(ds.rename({'time':'date'}), sic)
    #ds.attrs['nsidc_g02202v3'] = sic.nsidc_g02202v3 

    # Using the nsidc function requires changing the latitude and longitude
    # from that dataset from character to float and from minutes and seconds
    # to degrees.

    return ds
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
          timecoord.append(datetime.datetime(2013,3,int(ds.calday[ii]))+datetime.timedelta(seconds=ds.meantime.values[ii]))
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

      # Include NSIDC sea ice concentrations
      if loadice:
        sic = nsidc(lat=ds.lat,lon=ds.lon)
        ds = xr.Dataset.merge(ds, sic)
        ds.attrs['nsidc_g02202v3'] = sic.nsidc_g02202v3
        sic = nsidc(lat=ds.lat,lon=ds.lon,dataset='g02202v4')
        ds = xr.Dataset.merge(ds, sic)
        ds.attrs['nsidc_g02202v4'] = sic.nsidc_g02202v4

      lstaccdat.append(ds)

    return lstaccdat
################################################################################
def oden(filename):
    """
    This function reads the data provided by john.prytherch@misu.su.se from 
    the MISU station on the Oden icebreaker (during the ACSE and AO16 campaign)
    and outputs an Xarray dataset. It takes as argument the netcdf file name
    provided for each campaign.

    Author : virginie.guemas@meteo.fr - October 2020
    Modified : Generalized from the acse function into the oden one 
                                      - Virginie Guemas - March 2021
    """
    sys.path.append(rootpath+'ACSE/Version2/')
    import acse_info

    ds=xr.open_dataset(filename,drop_variables='doy')
    # I drop doy because there is a bug in its conversion to a time series and we do not need it because we have a correct time index.

    for var in ds.keys():
      if acse_info.acse_names(var) != '':
        ds[var].attrs = {'long_name':acse_info.acse_names(var),'units':ds[var].units}
      if acse_info.acse_units(var) != '':
        ds[var].attrs = {'long_name':ds[var].long_name,'units':acse_info.acse_units(var)}

    return ds
################################################################################
def acse():
    """
    This function reads the ACSE campaign data and outputs an Xarray dataset. 

    Author : virginie.guemas@meteo.fr - October 2020
    Modified : Relies on the oden function - Virginie Guemas - March 2021
    """

    ds = oden(rootpath+'/ACSE/Version2/ACSE_CANDIFLOS_fluxes_Oden_20140710_v5.0.nc')

    # Include NSIDC sea ice concentrations
    if loadice:
      sic = nsidc(lat=ds.lat,lon=ds.lon)
      ds = xr.Dataset.merge(ds, sic)
      ds.attrs['nsidc_g02202v3'] = sic.nsidc_g02202v3
      sic = nsidc(lat=ds.lat,lon=ds.lon,dataset='g02202v4')
      ds = xr.Dataset.merge(ds, sic)
      ds.attrs['nsidc_g02202v4'] = sic.nsidc_g02202v4

    return ds
################################################################################
def ao16():
    """
    This function reads the AO16 campaign data and outputs an Xarray dataset. 

    Author : virginie.guemas@meteo.fr - March 2021
    """

    ds = oden(rootpath+'/AO16/AO2016_foremast_30min_v4_0.nc')

    # We can not load NSIDC sea ice concentration because we don't have the 
    # ship position but we already have AMSR sea ice concentration in the
    # AO16 dataset itself.

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
      timecoord.append(datetime.datetime(year=2008,month=months[ii],day=int(ds.day[ii]),hour=int(ds.hour[ii]),minute=int(ds.minute[ii]),second=int(ds.second[ii])))
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

    # Include NSIDC sea ice concentrations
    if loadice:
      sic = nsidc(lat=ds.latitude,lon=ds.longitude)
      ds = xr.Dataset.merge(ds, sic)
      ds.attrs['nsidc_g02202v3'] = sic.nsidc_g02202v3
      sic = nsidc(lat=ds.latitude,lon=ds.longitude,dataset='g02202v4')
      ds = xr.Dataset.merge(ds, sic)
      ds.attrs['nsidc_g02202v4'] = sic.nsidc_g02202v4

    return ds
################################################################################
def stable():
    """
    This function reads the STABLE campaign data and outputs an Xarray dataset. 

    Author : virginie.guemas@meteo.fr - June 2022
    Modified : Merging of different files and sea ice concentration - Virginie Guemas - Octoboer 2022
    """

    lstdates=['20130310','20130325','20130326']

    for date in lstdates:
      fileflight=glob(rootpath+'STABLE/reduced/'+date+'*lead-parallel_avg-val_fluxes.tab')
      # File containing the fluxes as well as the wind speed and temperature at flight level.
      ds = readstablefile(fileflight[0], idcolumn=True)  

      # Additional files (one per leg) contain the surface temperature and pressure and
      # humidity at flight level
      lstlegs=sorted(glob(rootpath+'STABLE/reduced/'+date+'*p5_lead-parallel.tab'))
      for fileleg in lstlegs:
        ds1 = readstablefile(fileleg, idcolumn=False)
        # We need the average along each leg
        ds2 = postprostable(ds1)
        if fileleg==lstlegs[0]:
          dslegs = ds2
        else:
          dslegs = xr.concat((dslegs,ds2),dim='time')

      # Merge the datasets from all files from this flight
      ds=xr.merge([ds,dslegs],compat='override')

      # Final concatenation of the flights
      if date == lstdates[0]:
        dstot = ds
      else:
        dstot = xr.concat((dstot,ds),dim='time')

    # Include NSIDC sea ice concentrations
    if loadice:
      sic = nsidc(lat=dstot.Latitude,lon=dstot.Longitude)
      dstot = xr.Dataset.merge(dstot, sic)
      dstot.attrs['nsidc_g02202v3'] = sic.nsidc_g02202v3
      sic = nsidc(lat=dstot.Latitude,lon=dstot.Longitude,dataset='g02202v4')
      dstot = xr.Dataset.merge(dstot, sic)
      dstot.attrs['nsidc_g02202v4'] = sic.nsidc_g02202v4

    return dstot
################################################################################
def postprostable(dsin, fluxes = False):
    """
    This function computes the turbulent fluxes, wind speed and ...
    from raw STABLE files for each leg. It takes an Xarray Dataset storing all
    the variables from a STABLE leg file and outputs an Xarray Dataset containing
    turbulent fluxes and ... following the conventions of the flight files.

    Author : virginie.guemas@meteo.fr - June 2022
    Modified : Aditional output variables and fluxes option - Virginie Guemas - October 2022
    """
    # Compute Specific humidity at flight level
    dsin['Q']=meteolib.Q(rh=dsin.RH,T=dsin.TTT + 273.15, P=dsin.PPPP)
    dsin['Q'].attrs={'long_name':'Specific humidity from RH, TTT and PPPP at nose boom','units':'kg kg'}
    # Compute the average Longitude, Latitude, Altitude, Zonal, Meridional, Vertical
    # Wind speed, Temperature, Pressure, Relative Humidity, Specific Humidity and Surface 
    # temperature # for each leg
    dsout=xr.DataArray.mean(dsin)
    # Put back the time in the array (removed by the DataArray.mean operation)
    dsout=dsout.assign_coords(time=xr.DataArray.mean(dsin.time))
    # Restore attributes 
    for var in dsout.data_vars:
      dsout[var].attrs=dsin[var].attrs
   
    if fluxes : # Horizontal wind speed, Potential temperature and Air density are
      # also missing if fluxes have not been computed
      #
      # Compute Horizontal wind speed
      dsout['Wind_vel_horizontal']=(dsout.U**2+dsout.V**2)**0.5
      dsout['Wind_vel_horizontal'].attrs={'long_name':'Horizontal wind speed from mean U and V from five-hole probe at nose-boom','units':'m s-1'}
      # Compute Horizontal wind speed standard deviation
      dsout['Wind_vel_std_dev']=xr.DataArray.std((dsin.U**2+dsin.V**2)**0.5, ddof=1)
      dsout['Wind_vel_std_dev'].attrs={'long_name':'Horizontal wind speed standard deviation along the leg from five-hole probe at nose-boom','units':'m s-1'}
      # Not completely satisfactory - up to 0.05 difference with the one provided by Lupkes et al. Why?

      # Compute Air density
      dsout['Density_air']=meteolib.RHO(q=dsout.Q,T=dsout.TTT + 273.15, P=dsout.PPPP)
      dsout['Density_air'].attrs={'long_name':'Air density from Q, TTT and PPPP at nose boom','units':'kg m-3'}

      # Compute Potential air temperature
      dsout['THETA']=meteolib.Theta(q=dsout.Q,T=dsout.TTT + 273.15, P=dsout.PPPP)
      dsout['THETA'].attrs={'long_name':'Potential air temperature from Q, from TTT and PPPP at nose boom','units':'K'}
      # Compute Potential air temperature standard deviation
      dsout['THETA_std_dev']=xr.DataArray.std(meteolib.Theta(T=dsin.TTT + 273.15, P=dsin.PPPP))
      dsout['THETA_std_dev'].attrs={'long_name':'Potential air temperature standard deviation at nose boom','units':'K'}
      # Not completely satisfactory - up to 0.02 difference with the one provided by Lupkes et al. Why?

      # Compute Zonal momentum flux
      # dsout['Mom_flux_x']=xr.cov(dsin.U,dsin.W) Same result as below
      # dsout['Mom_flux_x']=dsout.Density_air * xr.DataArray.mean((dsin.U -xr.DataArray.mean(dsin.U))*(dsin.W-xr.DataArray.mean(dsin.W)))
      # dsout['Mom_flux_x'].attrs={'long_name':'West-East momentum flux estimated through the eddy-covariance method','units':'N m-2'}
      # Nothing to do with what is provided by Lupkes ??????????????????????????????????????????????????????????
      #
      # Total momentum flux
      dsout['Mom_flux_tot']=(xr.cov(dsin.U,dsin.W)**2 + xr.cov(dsin.U,dsin.W)**2)**0.5
      dsout['Mom_flux_tot'].attrs={'long_name':'Momentum flux estimated through the eddy-covariance method','units':'N m-2'}
      # Nothing to do with what is provided by Lupkes ??????????????????????????????????????????????????????????

    return dsout
################################################################################
def readstablefile(filename, idcolumn=True):
    """ 
    This function reads a file from the STABLE campaign and stores it in an 
    Xarray dataset. It takes two arguments:

    - filename = a complete path and filename to the text file containing the STABLE table to read. 
No default. 
    - idcolumn = True/False whether there is an ID column or not. Default : True 

    Author : virginie.guemas@meteo.fr - June 2022
    """

    if not(os.path.exists(filename)):
      sys.exit(filename+' does not exist')

    sys.path.append(rootpath+'STABLE/')
    import stable_info
    # Where to find a subsidiary file where I stored the information I found on the STABLE 
    # variables

    if idcolumn:
      firstcol = 2
    # In the STABLE files, either the first column is the leg name and the second is the date
    else:
    # or the first column is the date
      firstcol = 1

    lines = readtabtxt(filename,endline='\n',splitcolumn='\t')

    ds=xr.Dataset() 
    # One dataset containing all variables available for each flight or leg depending on the file
    for ifld in range(firstcol,len(lines[0])): # we skip the ID and date/time columns
      values=[]
      for iline in range(1,len(lines)): # First line is variable name and unit
        values.append(float(lines[iline][ifld]))
        # values along a column organised into an array
      fldid=lines[0][ifld].split("[")[0].strip().replace(' ','_')
      # we strip the unit from the name in the first line and replace spaces by _ to get a variable id 
      array=xr.DataArray(values,dims=('time'),attrs={'long_name':stable_info.stable_names(fldid),'units':stable_info.stable_units(fldid)})
      # the list of values becomes an Xarray DataArray
      ds[fldid]=array
      # the Xarray DataArray is included into the flight or leg Dataset

    timecoord=[]
    for iline in range(1,len(lines)):
      timecoord.append(datetime.datetime.fromisoformat(lines[iline][firstcol-1]))
      # date/time column provides the complete date YYYY-MM-DDTHH:MM:SS
    ds['time']=timecoord

    if idcolumn: # Flight files also have leg ids as a second time axis
      values=[]
      for iline in range(1,len(lines)): 
        values.append(float(lines[iline][0][8:10]))
      # The ids contain the date and the leg number, we keep here the leg number
      array=xr.DataArray(values,dims=('time'),attrs={'long_name':'leg number'})
      ds['leg']=array
      
      # For unknown reasons, the chronological order of flight files has been lost
      ds=ds.sortby('leg')
  
    return ds
################################################################################
def readtabtxt(filename,endline='',splitcolumn=' '):
    """
    This function reads a table stored in a text file (ASCII). It returns a list
    of lines, which are list of column values.
    It takes three arguments:

    - filename = a complete path and filename to the text file containing the table to read. No default. 
    - endline = an optional string present at the end of each line. Example : \'\\n\'. Default : empty.
    - splitcolumn = an optional string present between each column. Example : \'\\t\'. Default : space.

    Author : virginie.guemas@meteo.fr - June 2022
    """

    if not(os.path.exists(filename)):
      sys.exit(filename+' does not exist')
 
    # Store the txt tables into a list of lines
    f = open(filename,mode='r',encoding='utf-8')
    lines=f.readlines()
    f.close()

    # Turn it into a list of list of column values
    for iline in range(len(lines)):
      lines[iline]=lines[iline].strip(endline)
      lines[iline]=lines[iline].split(splitcolumn)

    return lines
################################################################################
def nsidc(lat, lon, dataset='g02202v3', hemisphere='nh'):
    """
    This functon reads the NSIDC sea ice concentration datasets and returns sea 
    ice concentration time series within an Xarray Dataset. The time series 
    follow the track of a ship pathway provided by a time series of longitudes 
    and latitudes. It takes three arguments:

    - lat = a DataArray of latitudes (in degrees).
    - lon = a DataArray of longitudes (in degrees).
    - dataset = a NSIDC dataset id : 'g02202v3'/'g02202v4'. Default : 'g02202v3'.
    - hemisphere = 'nh' or 'sh' for Arctic or Antarctic respectively. Default : 'nh'.

    Author : virginie.guemas@meteo.fr - July 2022
    Modified : November 2022 - Virginie Guemas - Add option g02202v4
    """
    
    loose_grid_step = 20

    # g02202v3 dataset corresponds to NSIDC id g02202 version 3. 
    if dataset == 'g02202v3':
      # For any dataset option, we need the basename (before the date) and the suffix (after the date)
      basename = 'seaice_conc_daily_'+hemisphere+'_*_'
      suffix = '_v03r01.nc'
      # We need the list of variables to be loaded
      lstvars = ['seaice_conc_cdr','stdev_of_seaice_conc_cdr','goddard_merged_seaice_conc','goddard_nt_seaice_conc','goddard_bt_seaice_conc']
      # and the list of variables not to be loaded
      dropvars = ['projection','melt_onset_day_seaice_conc_cdr','qa_of_seaice_conc_cdr']

    # g02202v4 dataset corresponds to NSIDC id g02202 version 4. 
    elif dataset == 'g02202v4':
      # For any dataset option, we need the basename (before the date) and the suffix (after the date)
      basename = 'seaice_conc_daily_'+hemisphere+'_'
      suffix = '_v04r00.nc'
      # We need the list of variables to be loaded
      lstvars = ['cdr_seaice_conc','stdev_of_cdr_seaice_conc','nsidc_nt_seaice_conc','nsidc_bt_seaice_conc']
      # and the list of variables not to be loaded
      dropvars = ['melt_onset_day_cdr_seaice_conc','qa_of_cdr_seaice_conc','spatial_interpolation_flag','temporal_interpolation_flag','projection']
    else:
      sys.exit('This dataset option is not coded yet')

    if dataset == 'g02202v3' or dataset == dataset == 'g02202v4':
      # Here is a general comment describing the data to be included 
      comment = 'The CDR algorithm output is a rule-based combination of ice concentration estimates from two well-established algorithms: the NASA Team (NT) algorithm (Cavalieri et al. 1984) and NASA Bootstrap (BT) algorithm (Comiso 1986). The CDR algorithm blends the NT and BT output concentrations by selecting, for each grid cell, the higher concentration value. The CDR algorithm capitalizes on the strengths of each contributing algorithm to produce ice concentration fields that should be more accurate than those from either algorithm alone. The CDR begins in 1987 with DMSP SSM/I passive microwave data, rather than in 1978 with NASA Nimbus-7 SMMR data, because the complete processing history of SMMR brightness temperatures cannot be traced and therefore the CDR program requirement for transparency is not met.'
    else:
      sys.exit('This dataset option is not coded yet')
    
    # Provide new names in the campaign dataset with the NSIDC dataset id as a suffix
    newnames = {}
    for var in lstvars+['time','latitude','longitude']:
      newnames[var]=var+'_'+dataset

    for jt in range(len(lat)):
    # We loop over the (latitude, longitude) time series providing the ship location
      if dataset == 'g02202v3':
        filename = glob(rootpath+'NSIDC/'+dataset+'/data/'+basename + lat.time[jt].dt.strftime('%Y%m%d').data + suffix)
      elif dataset == 'g02202v4':
        filename = glob(rootpath+'NSIDC/'+dataset+'/data/*/aggregate/'+basename + lat.time[jt].dt.strftime('%Y').data + suffix)
      if len(filename)!=1 :
        sys.exit(('We expect one and only one file to be loaded for a single date. Here we get ',filename))
      # For each measurement date, we open the associated NSIDC file
      seaicefld = xr.open_dataset(filename[0],drop_variables=dropvars)

      if dataset == 'g02202v4':
        seaicefld = seaicefld.isel(tdim = lat.time[jt].dt.dayofyear-1 ) 
      # For the g02202v4 dataset, all the days of a year are stored in a single file. We select the correct day.
        seaicefld = seaicefld.rename({'x':'xgrid','y':'ygrid'})
      # Rename dimensions to the same as g02202v3, easier to handle dimensions below
      
      if dataset == 'g02202v3':
        seaicefld = seaicefld.squeeze(dim='time')
      # The NSIDC time dimension is not a dimension anymore but only a variable. The time dimension will
      # rather be the ship measurement time. 

      if np.isnan(lon[jt]) or np.isnan(lat[jt]) or filename[0] == rootpath+'NSIDC/g02202v3/data/seaice_conc_daily_nh_f13_19980311_v03r01.nc' or filename[0] == rootpath+'NSIDC/g02202v3/data/seaice_conc_daily_nh_f13_19980820_v03r01.nc':
        # If the location of the ship is unknown, 
        # we use a random point from the NSIDC grid
        seaicefld = seaicefld.isel({'ygrid':0,'xgrid':0}).reset_coords(['xgrid','ygrid']).drop_vars(['xgrid','ygrid'])
        # we set every variable to np.nan
        for var in lstvars+['latitude','longitude']:
          seaicefld[var]=np.nan
      else:
        for var in ['latitude', 'longitude']:
          seaicefld[var] = seaicefld[var].where(seaicefld[lstvars[0]]<=1.)
        # For latitudes, longitudes, sea ice concentrations and associated variables, we change the 
        # continent values to np.nan. We rely on sea ice concentration values to find continents for lat/lon.
        for var in lstvars:
          seaicefld[var] = seaicefld[var].where(seaicefld[var]<=1.)
        # Distance between the ship location and each point on the NSIDC grid
        dist = distance((lon[jt].values,lat[jt].values),(seaicefld.longitude.values,seaicefld.latitude.values))
        # Loose scanning of the grid to speed up the process, i.e. every loose_grid_step points
        left_shift = int(loose_grid_step/2)
        distbis = dist[left_shift:seaicefld.latitude.shape[0]:loose_grid_step,left_shift:seaicefld.latitude.shape[1]:loose_grid_step]
        # Location of the nearest neighbour to the ship on the loose grid
        (idyb,idxb) = np.where(distbis==np.nanmin(distbis))
        # Reduction of the NSIDC grid to an area around this approximate neareast neighbour
        halo = int(loose_grid_step/4)
        seaicefld = seaicefld.isel({'ygrid':slice(-halo+idyb[0]*loose_grid_step,(idyb[0]+1)*loose_grid_step+halo),'xgrid':slice(-halo+idxb[0]*loose_grid_step,(idxb[0]+1)*loose_grid_step+halo)})
        dist = dist[slice(-halo+idyb[0]*loose_grid_step,(idyb[0]+1)*loose_grid_step+halo),slice(-halo+idxb[0]*loose_grid_step,(idxb[0]+1)*loose_grid_step+halo)]
        # Location of the nearest neighbour to the ship
        (idy,idx) = np.where(dist==np.nanmin(dist))
        # Selection of the nearest neighbour and removal of the grid of the NSIDC file
        seaicefld = seaicefld.isel({'ygrid':idy[0],'xgrid':idx[0]}).reset_coords(['xgrid','ygrid']).drop_vars(['xgrid','ygrid'])
      # Turn the coordinates into variables (the coordinates will be the ship ones, the NSIDC are kept for check) and rename all variables to include the dataset id
      seaicefld = seaicefld.reset_coords(['time','latitude','longitude']).rename(newnames)
      # Concatenate the Xarray datasets extracted for each measurement time
      if 'seaice' in locals().keys():
        seaice = xr.concat([seaice,seaicefld],dim='time')
      else:
        seaice = seaicefld
      # Close the NSIDC netcdf file
      seaicefld.close()

    # Set the time axis to be the same as the campaign Xarray Dataset one
    seaice=seaice.assign_coords(time=lon.time)
    # Drop the general NSIDC attributes (but keep the ones for each variable) and keep a brief description
    # of the nsidc dataset
    seaice.attrs={'nsidc_'+dataset:comment} 

    return seaice
################################################################################
def distance(coord1, coord2):
    """
    This function computes the distance in kilometers between two points at the 
    surface of the Earth. It takes two arguments :
    - the coordinates of the first point (lon, lat) in degrees as a tuple
    - the coordinates of the second point (lon, lat) in degrees as a tuple
    The longitudes and latitudes can either be a single location for coord1
    and an array for coord2 or vice-versa or be two arrays with the same 
    dimensions.

    Author : virginie.guemas@meteo.fr - July 2022
    """

    (lon1, lat1) = coord1
    (lon2, lat2) = coord2

    if lon1.shape != lat1.shape or lon2.shape != lat2.shape :
      sys.exit("lon1 and lat1 should have the same dimensions on one hand as well as lon2 and lat2 on the other hand")

    if lon1.shape != lon2.shape :
      if lon1.shape!=() and lon2.shape!=():
        sys.exit("Either provide 2 coordinates array of same dimensions or one location and an array of any dimension.")

    R = 6371.0  # Earth radius

    # Latitude and longitude distance in radians
    dlon = np.radians(lon2) - np.radians(lon1)
    dlat = np.radians(lat2) - np.radians(lat1)

    #Haversine formula
    a = np.sin(dlat / 2)**2 + np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon / 2)**2
    distance = 2*R*np.arcsin(np.sqrt(a))

    return distance
################################################################################
def mosaic(site='tower', startdate='2019-10-09', enddate='2019-10-31'):
    """
    This function loads the MOSAiC data either from the tower, or from one 
    of the AFS station around. It takes three arguments :
    - site = 'tower' / 'asf30' / 'asf40' / 'asf50'. Default : 'tower'
    - startdate = YYYY-MM-DD. Example : '2019-10-05'
    - enddate = YYYY-MM-DD. Example : '2020-10-01'
    Please do not load the whole MOSAiC campaign at once. It takes too long.
    Data starts on 2019-10-05 and end on 2020-10-01.

    This function outputs an Xarray Dataset.

    Author : virginie.guemas@meteo.fr - October 2022
    """

    lstsites=('tower','asf40','asf50','asf30')
    if site not in lstsites: 
      sys.exit(('Argument site should be in ',lstsites))

    # Location of MOSAic files and basename
    directory1 = { 'tower':'Tower','asf30':'ASF/ASF30','asf40':'ASF/ASF40','asf50':'ASF/ASF50'}
    directory2 = { 'tower':'','asf30':'_asfs30','asf40':'_asfs40','asf50':'_asfs50'}
    path=rootpath+'MOSAiC/'+directory1[site]+'/ftpnoaa/3_level_archive'+directory2[site]+'/level3.4/'
    rootname = { 'tower':'metcity','asf30':'asfs30','asf40':'asfs40','asf50':'asfs50'}
    basename='mosseb.'+rootname[site]+'.level3.4.10min.'

    # Period covered by files
    start_date=datetime.date(int(startdate[0:4]),int(startdate[5:7]),int(startdate[8:10]))
    end_date=datetime.date(int(enddate[0:4]),int(enddate[5:7]),int(enddate[8:10]))
    print('MOSAiC data are loaded from ',start_date.strftime("%d %B %Y"),' to ',end_date.strftime("%d %B %Y"))

    # Loop on dates 
    date = start_date
    lstds=[]
    while date <= end_date:
      # Data are stored in daily files
      filename = path+basename+date.strftime("%Y%m%d")+'.000000.nc'
      # Whenever present (intermittency), files loaded as Xarray datasets
      if os.path.exists(filename):
        lstds.append(xr.open_dataset(filename))
      date += datetime.timedelta(days=1)    
     
    # Check whether the list of datasets contain any dataset
    if bool(lstds): 
      # Concatenation of datasets
      ds = xr.concat(lstds, dim = 'time')
    else:
      sys.exit('Not a single MOSAiC file over that period')

    if site == 'tower':
      # Creation of a height axis to reorganize data in (time, height) format
      height = [2, 6, 10]
      ds['height'] = xr.DataArray(height, dims=['height'], attrs={'units':'m'})

      # List of variables which are present as variables and variables_qc
      # lstvarsqc = ['temp','dew_point','rh','mixing_ratio','rhi','vapor_pressure','wspd_vec_mean','wdir_vec_mean','wspd_u_mean','wspd_v_mean','wspd_w_mean','temp_acoustic_mean','wspd_u_std','wspd_v_std','wspd_w_std','temp_acoustic_std','ustar'] # version level 2 october 2022
      lstvarsqc = ['temp','dew_point','rh','mixing_ratio','rhi','vapor_pressure','wspd_vec_mean','wdir_vec_mean','wspd_u_mean','wspd_v_mean','wspd_w_mean','wspd_u_std','wspd_v_std','wspd_w_std','turbulence']
      # List of variables which are present only as variables
      # lstvars = ['Hs','Cd','Tstar','zeta_level_n','WU_csp','WV_csp','UV_csp','WT_csp','UT_csp','VT_csp','phi_U','phi_V','phi_W','phi_T','phi_UT','epsilon_U','epsilon_V','epsilon_W','epsilon','Phi_epsilon','nSU','nSV','nSW','nST','NT','Phi_NT','Phix','DeltaU','DeltaV','DeltaT'] # version level 2 october 2022
      lstvars = ['Hs','Cd','ustar','Tstar','zeta_level_n','epsilon','sigU','sigV','sigW','wind_sector_qc_info']
      lstvars.extend(lstvarsqc)
      # List of variables which have dimensions (time, freq)
      #lstvarfreq = ['sUs','sVs','sWs','sTs','cWUs','cWVs','cUVs','cWTs','cUTs','cVTs'] # version level 2 october 2022
      lstvarfreq = [ ]
      lstvars.extend(lstvarfreq)
      # List of attributes that are present at each level but should be recomputed for the whole set, easy to compute afterwards so they are dropped
      lstattrs = ['min_val','max_val','avg_val','height_start','height_end','height_change_time','percent_missing']
      # Loop on variables to reshape
      for var in lstvars:
        if var in lstvarsqc:
          lstexts = ['', '_qc']
          if var in ['turbulence']:
            lstexts = ['_qc']
        else: 
          lstexts = ['']
        for ext in lstexts:
          # Dimensions differ whether we have (time) or (time, freq) initially
          if var in lstvarfreq:
            var0 = np.empty((len(ds['time']),60,len(height)))
            var0[:,:,:] = np.nan
          else:
            var0 = np.empty((len(ds['time']),len(height)))
            var0[:,:] = np.nan
          # Loop on measurement heights
          for hh in range(len(height)):
            # Initial names at each height
            name = var+'_'+str(height[hh])+'m'+ext
            # Storage of values in the output variables
            if var in lstvarfreq:
                var0[:,:,hh] = ds[name].values
            else:
              var0[:,hh] = ds[name].values
            # Keep the attributes of the last level
            attrs = ds[name].attrs
            # Remove the initial variables (at each height)
            ds=ds.drop(name)  

          # Store the output variable in a Xarray into the final dataset
          if var in lstvarfreq:
            ds[var+ext] = xr.DataArray(var0, dims=['time','freq','height'], attrs=attrs)
          else:
            if var == 'zeta_level_n':
              var = 'zeta' # Rename zeta variable 
            ds[var+ext] = xr.DataArray(var0, dims=['time','height'], attrs=attrs)
          # Drop the height from this atttribute
          ds[var+ext].attrs['location'] = 'met city tower' 

          # Modify the attributes that mention the height in the initial variables
          if var == 'zeta':
            ds[var].attrs['long_name'] = 'Monin-Obukhov stability parameter, z/L'
          if var == 'phi_U':
            ds[var].attrs['long_name'] = 'MO universal function for the standard deviations, U is the streamwise wind vector'
          if var == 'phi_V':
            ds[var].attrs['long_name'] = 'MO universal function for the standard deviations, V is the streamwise wind vector'
          if var == 'phi_W':
            ds[var].attrs['long_name'] = 'MO universal function for the standard deviations, W is the streamwise wind vector'
          if var == 'phi_T':
            ds[var].attrs['long_name'] = 'MO universal function for the standard deviations'
          if var == 'phi_UT':
            ds[var].attrs['long_name'] = 'MO universal function for the horizontal heat flux, U is the streamwise wind vector'
          if var == 'Phi_epsilon':
            ds[var].attrs['long_name'] = 'Monin-Obukhov universal function Phi_epsilon based on the median epsilon'
          if var == 'Phi_NT' :
            ds[var].attrs['long_name'] = 'Monin-Obukhov universal function Phi_NT'
          # Drop the uneless attributes
          for at in lstattrs:
            if (hasattr(ds[var+ext],at)) & (var!=lstvars[-1]):
              del ds[var+ext].attrs[at]

    # Include NSIDC sea ice concentrations
    if loadice:
      if site == 'tower':
        sic = nsidc(lat=ds.lat_tower,lon=ds.lon_tower,dataset='g02202v4')
      else:
        sic = nsidc(lat=ds.lat,lon=ds.lon,dataset='g02202v4')
      ds = xr.Dataset.merge(ds, sic)
      ds.attrs['nsidc_g02202v4'] = sic.nsidc_g02202v4

    return ds
################################################################################
def igp():
    """ 
    This function loads the IGP data.

    This function outputs an Xarray Dataset.

    Author : virginie.guemas@meteo.fr - October 2022
    """

    lstds=[]
    for jnumb in (292,):
      path=rootpath+'IGP/flight'+str(jnumb)+'/'
      filename = path+'core_masin_20180228_r001_flight'+str(jnumb)+'_igp-qc.nc'
      print(filename)
      lstds.append(xr.open_dataset(filename))
    ds = xr.concat(lstds, dim = 'time')

    return ds
