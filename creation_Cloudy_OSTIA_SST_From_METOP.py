
"""
Creation of groundtruthed cloudy SST dataset using OSTIA SST as groundtruth 
and METOP SST data as cloud masks.
Code associated with paper:
    Fablet et al., Data-driven Models for the Spatio-Temporal Interpolation of satellite-derived SST Fields,
    IEEE TCI, 2017
"""

import numpy as np
import netCDF4
from datetime import date, timedelta
from tqdm import tqdm
import matplotlib.pyplot as plt

# specify longitude, latitude of starting point to extract OSTIA SST grid
lon_start = 4.75  
lat_start = -44.5 

# specify size of extracted SST grid
num_lons = 600  
num_lats = 300

# specify extracting years
year_start = 2008
year_end = 2015

# resolution of OSTIA SST grid
res = 0.05

# converting coordinates of starting point to spatial 
pos_x = int((lon_start+180)/res)
pos_y = int((lat_start+90)/res)

# Exactring METOP missing data mask to create cloudy OSTIA
url_metop_mask = 'http://www.ifremer.fr/opendap/cerdap1/ghrsst/l3p/metop-a/avhrr/eur-l3p-glob_avhrr_metop_a/data/'
year_mask = 2009

delta = date(year_mask+1,1,1) - date(year_mask,1,1)

Metop_Cloud_Mask = np.zeros([delta.days,num_lats,num_lons])

for i in tqdm(range(delta.days)):
    date_tmp = date(year_mask,1,1) + timedelta(days=i)
    url_tmp = url_metop_mask + str(date_tmp.year)+'/'+str(format(i+1,'03d'))+'/'+str(date_tmp.year)+str(format(date_tmp.month,'02d'))+str(format(date_tmp.day,'02d'))+'-EUR-L3P_GHRSST-SSTsubskin-GLOB_AVHRR_METOP_A-eumetsat_sstglb_metop02_'+\
                     str(date_tmp.year)+str(format(date_tmp.month,'02d'))+str(format(date_tmp.day,'02d'))+'_000000-v01.7-fv01.0.nc?sea_surface_temperature[0:1:0]['+str(pos_y)+':1:'+str(pos_y+num_lats-1)+']['+str(pos_x)+':1:'+str(pos_x+num_lons-1)+']'
    dataset = netCDF4.Dataset(url_tmp,'r')
    sst_tmp = dataset.variables['sea_surface_temperature'][:]
    dataset.close()
    sst_metop = sst_tmp.data[0,:,:] - 273.15
    sst_metop[sst_tmp.mask[0,:,:]] = np.nan
    Metop_Cloud_Mask[i,:,:] = np.flipud(sst_metop)

# Extracting cloudy OSTIA
url_ostia = 'https://podaac-opendap.jpl.nasa.gov:443/opendap/allData/ghrsst/data/L4/GLOB/UKMO/OSTIA/'
delta = date(year_end+1,1,1) - date(year_start,1,1)
Ostia_Obs = np.zeros([delta.days,num_lats,num_lons])
count = 0
for year in range(year_start,year_end+1):
    print "extracting year " + str(year)
    d1 = date(year,1,1)
    d2 = date(year+1,1,1)
    delta = d2 - d1         # timedelta
    for i in tqdm(range(delta.days)):
        date_tmp = d1 + timedelta(days=i)
        url_tmp = url_ostia+str(date_tmp.year)+'/'+str(format(i+1,'03d'))+'/'+str(date_tmp.year)+str(format(date_tmp.month,'02d'))+str(format(date_tmp.day,'02d'))+'-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA.nc.bz2?time[0:1:0],analysed_sst[0:1:0]['+\
                         str(pos_y)+':1:'+str(pos_y+num_lats-1)+']['+str(pos_x)+':1:'+str(pos_x+num_lons-1)+']'
        dataset = netCDF4.Dataset(url_tmp,'r')
        sst_tmp = dataset.variables['analysed_sst'][:]
        dataset.close()
        sst_ostia = sst_tmp.data[0,:,:] - 273.15
        sst_ostia[sst_tmp.mask[0,:,:]] = np.nan
        sst_ostia = np.flipud(sst_ostia)
        sst_ostia[np.isnan(Metop_Cloud_Mask[i,:,:])] = np.nan
        Ostia_Obs[count,:,:] = sst_ostia
        count = count + 1 
 
### Write into netCDF file ######################################################################   

### Write metop mask
metop_mask_file = 'metop_mask.nc'
fid = netCDF4.Dataset(metop_mask_file,'w', format='NETCDF4')
# Define the dimensions
lat = fid.createDimension('lat', num_lats)       
lon = fid.createDimension('lon', num_lons) 
time = fid.createDimension('time', len(Metop_Cloud_Mask))

fid.title = "Metop Cloud Mask"
fid.summary = "Cloud mask created from Metop"
fid.acknowledgment = "The data from Ifremer"

nc_var = fid.createVariable('times', 'f8',('time'))
fid.variables['times'][:] = np.arange(len(Metop_Cloud_Mask))
fid.variables['times'].standard_name='day'
fid.variables['times'].units='day'
fid.variables['times'].comment='Day_Length'

nc_var = fid.createVariable('lats', 'f8',('lat'))
fid.variables['lats'][:] = np.arange(lat_start,lat_start+num_lats*res,res)
fid.variables['lats'].standard_name='latitude'
fid.variables['lats'].units='degrees_north'
fid.variables['lats'].comment='Latitude'

nc_var = fid.createVariable('lons', 'f8',('lon'))
fid.variables['lons'][:] = np.arange(lon_start,lon_start+num_lons*res,res)
fid.variables['lons'].standard_name='longitude'
fid.variables['lons'].units='degrees_east'
fid.variables['lats'].comment='Longitude'

nc_var = fid.createVariable('sst', 'i2',('time','lat', 'lon'),fill_value=-32768)
metop_tmp = np.ma.masked_array(Metop_Cloud_Mask,np.isnan(Metop_Cloud_Mask))
fid.variables['sst'][:] = metop_tmp*100
fid.variables['sst'].units='Celcius'
fid.variables['sst'].scale_factor= 0.01
fid.variables['sst'].valid_min = -300
fid.variables['sst'].valid_max = 4500
fid.variables['sst'].comment='analysed_sst'    
fid.close()

### Write cloudy ostia
cloudy_ostia_file = 'cloudy_ostia.nc'
fid = netCDF4.Dataset(cloudy_ostia_file,'w', format='NETCDF4')
# Define the dimensions
lat = fid.createDimension('lat', num_lats)        
lon = fid.createDimension('lon', num_lons)
time = fid.createDimension('time', len(Ostia_Obs)) 

fid.title = "Cloudy Ostia SST"
fid.summary = "Ostia cloud created by applying Metop mask"
fid.acknowledgment = "The data from CMEMS"

nc_var = fid.createVariable('times', 'f8',('time'))
fid.variables['times'][:] = np.arange(len(Ostia_Obs))
fid.variables['times'].standard_name='day'
fid.variables['times'].units='day'
fid.variables['times'].comment='Day_Length'

nc_var = fid.createVariable('lats', 'f8',('lat'))
fid.variables['lats'][:] = np.arange(lat_start,lat_start+num_lats*res,res)
fid.variables['lats'].standard_name='latitude'
fid.variables['lats'].units='degrees_north'
fid.variables['lats'].comment='Latitude'

nc_var = fid.createVariable('lons', 'f8',('lon'))
fid.variables['lons'][:] = np.arange(lon_start,lon_start+num_lons*res,res)
fid.variables['lons'].standard_name='longitude'
fid.variables['lons'].units='degrees_east'
fid.variables['lats'].comment='Longitude'

nc_var = fid.createVariable('sst', 'i2',('time','lat', 'lon'),fill_value=-32768)
ostia_tmp = np.ma.masked_array(Ostia_Obs,np.isnan(Ostia_Obs))
fid.variables['sst'][:] = ostia_tmp*100
fid.variables['sst'].units='Celcius'
fid.variables['sst'].scale_factor= 0.01
fid.variables['sst'].valid_min = -300
fid.variables['sst'].valid_max = 4500
fid.variables['sst'].comment='analysed_sst'    
fid.close()
                        
# display example
index = 2
plt.figure()
plt.imshow(Ostia_Obs[index,:,:],aspect='auto',cmap='jet')
plt.colorbar()
