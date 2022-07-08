timecoord=[]
      if flight == 'FAAM': 
        for ii in range(len(values)):
          timecoord.append(datetime.datetime(2013,3,ds.calday[ii])+datetime.timedelta(seconds=ds.meantime.values[ii]))
          # Creation of time axis from the calday and meantime variables
