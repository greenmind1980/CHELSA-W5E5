#! /usr/bin/env python

# This file is part of chelsa_w5e5.
#
# chelsa_w5e5 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# chelsa_w5e5 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with chelsa_w5e5.  If not, see <https://www.gnu.org/licenses/>.


# ***************************************
# import libraries
# ***************************************

import saga_api 
from functions.saga_functions import *
import sys
import os 
import argparse 
import datetime 
import os.path 
import cdsapi 
import psutil 
import shutil 
import gdal
import subprocess
import xarray as xr
import cfgrib


# ***************************************
# Global variables and settings
# ***************************************
saga_api.SG_Set_History_Depth(0)
process = psutil.Process(os.getpid())


# set times of the day and levels
times  = [
            '00:00','01:00','02:00',
            '03:00','04:00','05:00',
            '06:00','07:00','08:00',
            '09:00','10:00','11:00',
            '12:00','13:00','14:00',
            '15:00','16:00','17:00',
            '18:00','19:00','20:00',
            '21:00','22:00','23:00'
        ]
levels = [
            '1','2','3',
            '5','7','10',
            '20','30','50',
            '70','100','125',
            '150','175','200',
            '225','250','300',
            '350','400','450',
            '500','550','600',
            '650','700','750',
            '775','800','825',
            '850','875','900',
            '925','950','975',
            '1000'
        ]
plevels_temp = [
                '750',
                '775','800','825',
                '850','875','900',
                '925','950','975',
                '1000'
            ]


# *************************************************
# Get the command line arguments
# *************************************************
ap = argparse.ArgumentParser(
    description='''# This python code for CHELSA_V2.1_W5E5
the code is adapted to the ERA5 reanalysis and automatically downloads
data from the CDS via the CDS api. It aggregates hourly data to daily
resolution and then runs the CHELSA algorithm for min-, max-, and mean 
temperature, total downwelling solar radiation, and 
optional surface precipitation. The output directory needs the following 
subfolders: /rsds, /pr, /tasmin, /tasmax, /tas.
Dependencies for ubuntu_18.04:
libwxgtk3.0-dev libtiff5-dev libgdal-dev libproj-dev 
libexpat-dev wx-common libogdi3.2-dev unixodbc-dev
g++ libpcre3 libpcre3-dev wget swig-4.0.1 python2.7-dev 
software-properties-common gdal-bin python-gdal 
python2.7-gdal libnetcdf-dev libgdal-dev
python-pip cdsapi saga_gis-7.6.0
All dependencies are resolved in the chelsa_V2.1.cont singularity container
Tested with: singularity version 3.3.0-809.g78ec427cc
''',
    epilog='''author: Dirk N. Karger, dirk.karger@wsl.ch, Version 2.1'''
)

# collect the function arguments
ap.add_argument('-b','--year', type=int, help="year, integer")
ap.add_argument('-c','--month', type=int, help="month, integer")
ap.add_argument('-d','--day', type=int,  help="day, integer")
ap.add_argument('-t1','--array_task_id', type=int, help= 'slurm array task id to identify the respective band')
ap.add_argument('-i','--input', type=str, help="input directory, string")
ap.add_argument('-o','--output', type=str,  help="output directory, string")
ap.add_argument('-t','--temp', type=str, help="root for temporary directory, string")
ap.add_argument('-s','--srad', type=str, help="srad input directory, string")
ap.add_argument('-w','--w5e5', type=str, help="w5e5 input directory, string")
ap.add_argument('-cc', '--ccover',type=str, help='directory containing the cloudcover')
ap.add_argument('-l', '--lapse',type=str, help='directory containing the lapserates')
ap.add_argument('--calc_pr', action='store_true', help="calculate pr (default: do not))")
ap.add_argument('--calc_rsds', action='store_true', help="calculate rsds (default: do not))")
ap.add_argument('--calc_tas', action='store_true', help="calculate tas (default: do not))")
ap.add_argument('--calc_tasmax', action='store_true', help="calculate tasmax (default: do not))")
ap.add_argument('--calc_tasmin', action='store_true', help="calculate tasmin (default: do not))")

ap.add_argument('-pr', '--pr',type=str, help='file containing pr, format ncdf, string')
ap.add_argument('-hu', '--hurs',type=str, help='file containing hurs, format ncdf, string')
ap.add_argument('-ta', '--tas',type=str, help='file containing tas, format ncdf, string')
ap.add_argument('-tx', '--tasmax',type=str, help='file containing tasmax, format ncdf, string')
ap.add_argument('-tn', '--tasmin',type=str, help='file containing tasmin, format ncdf, string')
ap.add_argument('-sr', '--ssrd',type=str, help='file containing ssrd, format ncdf, string')

ap.add_argument('-u', '--uwind',type=str, help='era5 file containing 10m_u_component_of_wind for this year, format ncdf, string')
ap.add_argument('-v', '--vwind',type=str, help='era5 file containing 10m_v_component_of_wind for this year, format ncdf, string')
args = ap.parse_args()

# *************************************************
# Get times from arguments
# *************************************************
print(args)

year     = "%02d" % args.year
day      = "%02d" % args.day
month    = "%02d" % args.month
t        = "%0d"  % args.array_task_id
month1   = "%01d" % args.month

calc_pr = args.calc_pr
calc_rsds = args.calc_rsds
calc_tas = args.calc_tas
calc_tasmax = args.calc_tasmax
calc_tasmin = args.calc_tasmin
variable_long_name = {
    'pr':     'Precipitation',
    'rsds':   'Surface Downwelling Shortwave Radiation',
    'tas':    'Daily Mean Near-Surface Air Temperature',
    'tasmax': 'Daily Maximum Near-Surface Air Temperature',
    'tasmin': 'Daily Minimum Near-Surface Air Temperature',
}

# get day of the year
dt = datetime.datetime(args.year, args.month, args.day)
tt = dt.timetuple()
dayofyear = tt.tm_yday
dayofyear = "%0d" % dayofyear

# *************************************************
# Set the directories from arguments
# *************************************************
INPUT = args.input
OUTPUT= args.output
TEMP  = args.temp
SRAD  = args.srad
W5E5  = args.w5e5
LAPSE = args.lapse
CCOVER= args.ccover
tlapse_mean_path = LAPSE + 'CHELSA_tz_' + day + '_' + month + '_' + year + '_V.2.1.tif'

TEMP = TEMP + 'w5e5' + year + month + day + '/'
if os.path.isdir(TEMP):
    shutil.rmtree(TEMP)

if not os.path.exists(TEMP):
    os.makedirs(TEMP)


# **************************************************
# Download all necessary files from ERA5 using CDSAPI
# **************************************************
print(datetime.datetime.now())

subprocess.call(['gdal_translate', '-b', t, W5E5 + args.ssrd, TEMP + 'ssrd.tif'])
subprocess.call(['gdal_translate', '-b', t, W5E5 + args.pr, TEMP + 'pr.tif'])
subprocess.call(['gdal_translate', '-b', t, W5E5 + args.hurs, TEMP + 'hurs.tif'])
subprocess.call(['gdal_translate', '-b', t, W5E5 + args.tas, TEMP + 'tas.tif'])
subprocess.call(['gdal_translate', '-b', t, W5E5 + args.tasmax, TEMP + 'tasmax.tif'])
subprocess.call(['gdal_translate', '-b', t, W5E5 + args.tasmin, TEMP + 'tasmin.tif'])


use_mistral = False
if use_mistral:
    # get daily mean u and v wind
    for ipath, opath in [(args.uwind, TEMP + 'u.nc'), (args.vwind, TEMP + 'v.nc')]:
        # select all hours of the day
        subprocess.call(['cdo', '-f', 'nc', 'seldate,%s-%s-%s'%(year, month, day), ipath, opath])
        # use xarray to compute daily mean to also remove time dimension
        xr.open_dataset(opath).mean(dim='time').to_netcdf(opath)
else:
    c = cdsapi.Client()

    # u wind
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type':'reanalysis',
            'variable':'10m_u_component_of_wind',
            'year':year,
            'month':month,
            'day':day,
            'time':[
                '00:00','01:00','02:00',
                '03:00','04:00','05:00',
                '06:00','07:00','08:00',
                '09:00','10:00','11:00',
                '12:00','13:00','14:00',
                '15:00','16:00','17:00',
                '18:00','19:00','20:00',
                '21:00','22:00','23:00'
            ],
            'format':'netcdf'
        },
        TEMP+'u.nc')

    # v wind
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type':'reanalysis',
            'variable':'10m_v_component_of_wind',
            'year':year,
            'month':month,
            'day':day,
            'time':[
                '00:00','01:00','02:00',
                '03:00','04:00','05:00',
                '06:00','07:00','08:00',
                '09:00','10:00','11:00',
                '12:00','13:00','14:00',
                '15:00','16:00','17:00',
                '18:00','19:00','20:00',
                '21:00','22:00','23:00'
            ],
            'format':'netcdf'
        },
        TEMP+'v.nc')


# ************************************************
# Script
# ************************************************
if __name__ == '__main__':

    saga_api.SG_Get_Data_Manager().Delete_All()  #
    Load_Tool_Libraries(True)
    print(datetime.datetime.now())
    print('functions loaded')
    track_memory_usage('RUN STARTED')


    ### calculate windeffect *******************************************************************************************
    ## import the files
    hurs = import_gdal(TEMP + 'hurs.tif')
    uwind = import_ncdf(TEMP + 'u.nc')
    vwind = import_ncdf(TEMP + 'v.nc')
    track_memory_usage('windeffect start')
    ## get means pof winds
    if use_mistral:
        uwind = uwind.Get_Grid(0).asGrid()
        vwind = vwind.Get_Grid(0).asGrid()
    else:
        uwind = get_mean(uwind)
        vwind = get_mean(vwind)

    # get mean temperature
    tmean = import_gdal(TEMP + 'tas.tif')
    tmean = change_data_storage(tmean)
    cblev = grid_calculator(tmean,hurs,'(20+((a-273.15)/5))*(100-b)')
    cblev = closegaps(cblev)
    track_memory_usage('get means of winds')
    ## change longitudinal range
    uwind = change_latlong(uwind)
    vwind = change_latlong(vwind)
    #cblev = change_data_storage(cblev)
    cblev = calc_geopotential(cblev)
    track_memory_usage('change logitudinal range')
    ## set coordinate reference system
    set_2_latlong(uwind)
    set_2_latlong(vwind)
    set_2_latlong(cblev)
    ## change to shapefile for projection
    uwind_shp = gridvalues_to_points(uwind)
    vwind_shp = gridvalues_to_points(vwind)
    cblev_shp = gridvalues_to_points(cblev)
    track_memory_usage('change to shapefile for projection')
    ## reproject to mercator projection
    uwind_shp = reproject_shape(uwind_shp)
    vwind_shp = reproject_shape(vwind_shp)
    ## multilevel b spline
    demproj   = load_sagadata(INPUT + 'dem_merc3.sgrd')
    uwind_ras = multilevel_B_spline(uwind_shp,demproj)
    vwind_ras = multilevel_B_spline(vwind_shp,demproj)
    direction = polar_coords(uwind_ras,vwind_ras)
    track_memory_usage('multilevel b - spline')
    dem1      = load_sagadata(INPUT + 'dem_latlong3.sgrd')
    dem       = change_data_storage(dem1)
    dem_geo   = calc_geopotential(dem)
    cblev_ras = multilevel_B_spline(cblev_shp, dem)
    windef    = windeffect(direction,demproj)
    saga_api.SG_Get_Data_Manager().Delete(demproj)
    windef1   = proj_2_latlong(windef,dem)
    export_geotiff(windef1, 'we')
    track_memory_usage('windeffect calculated')

    # import dems
    dem_low  = load_sagadata(INPUT + 'orography.sgrd')
    dem_low  = change_data_storage(dem_low)
    dem_low  = change_latlong(dem_low)
    dem_high = load_sagadata(INPUT + 'dem_latlong.sgrd')
    dem_high = change_data_storage(dem_high)
    track_memory_usage('import dems')


    ### calculate surface precipitation ********************************************************************************
    if calc_pr:
        prec = import_gdal(TEMP + 'pr.tif')
        prec = grid_calculator_simple(prec,'a*86400')
        dummy= load_sagadata(INPUT + 'dummy_W5E5.sgrd')
        prec = resample(prec,dummy)
        prec = closegaps(prec)
        # correct wind effect by boundary layer height
        pblh           = grid_calculatorX(dem_low,cblev, 'a+(b/9.80665)')
        dist2bound     = calc_dist2bound(dem_high, pblh)
        maxdist2bound  = resample_up(dist2bound,pblh,7)
        maxdist2bound2 = invert_dist2bound(dist2bound,maxdist2bound)
        track_memory_usage('corrected windef by pblh')

        wind_cor    = grid_calculatorX(maxdist2bound2,windef1,'(b/(1-a/9000))')
        patchgrid   = load_sagadata(INPUT + 'patch.sgrd')

        exp_index   = load_sagadata(INPUT + 'expocor.sgrd')
        wind_cor    = grid_calculatorX(exp_index,windef1,'a*b')
        wind_cor    = patching(wind_cor, patchgrid)
        wind_coarse = resample_up(wind_cor,dummy,4)
        wind_coarse = closegaps(wind_coarse)
        track_memory_usage('windeffect corrected')

        # downscale precipitation and export
        wind_coarse25 = resample_up(wind_cor,pblh,4)
        prec2 = downscale_precip(wind_coarse25,wind_coarse,prec,'prec',1)

        precip = downscale_precip(wind_cor,wind_coarse25,prec2,'total surface precipitation',3)
        precip_ccINT = convert2uinteger10(precip, variable_long_name['pr'])
        export_geotiff(precip_ccINT, 'pr')
        track_memory_usage('precipitation downscaled')

        # clean up memory
        saga_api.SG_Get_Data_Manager().Delete(precip_ccINT)
        saga_api.SG_Get_Data_Manager().Delete(precip)
        saga_api.SG_Get_Data_Manager().Delete(wind_cor)
        saga_api.SG_Get_Data_Manager().Delete(maxdist2bound2)
        saga_api.SG_Get_Data_Manager().Delete(maxdist2bound)
        saga_api.SG_Get_Data_Manager().Delete(dist2bound)
        saga_api.SG_Get_Data_Manager().Delete(exp_index)
        saga_api.SG_Get_Data_Manager().Delete(windef1)
        track_memory_usage('cleaned up memory')


    ### downscale temperatures *****************************************************************************************
    variable_names = []
    if calc_tas: variable_names.append('tas')
    if calc_tasmax: variable_names.append('tasmax')
    if calc_tasmin: variable_names.append('tasmin')

    # get lapse rate and masks
    if variable_names:
        if os.path.isfile(tlapse_mean_path):
            tlapse_mean = import_gdal(tlapse_mean_path)
        else:
            tas = import_ncdf(TEMP + 't.nc')
            zg  = import_ncdf(TEMP + 'z.nc')
    
            ta950 = tas.Get_Grid(0).asGrid()
            ta850 = tas.Get_Grid(1).asGrid()
    
            zg950 = zg.Get_Grid(0).asGrid()
            zg850 = zg.Get_Grid(1).asGrid()
    
            tlapse_mean = tlapse(zg950, zg850, ta950,
                   ta850, '(d*9.80665-c*9.80665)/(b-a)')
    
            tlapse_mean = change_latlong(tlapse_mean)
            tlapse_mean = resample(tlapse_mean, dem_low)
    
            export_geotiff(tlapse_mean, 'tz')

        tlapse_mean = change_data_storage(tlapse_mean)
        oceans      = load_sagadata(W5E5 + 'oceans.sgrd')
        landseamask = load_sagadata(W5E5 + 'landseamask.sgrd')
        track_memory_usage('downscale temperature start')

    for variable_name in variable_names:
        # load and cast data
        tmean = import_gdal(TEMP + variable_name + '.tif')
        tmean = change_data_storage(tmean)
        track_memory_usage('loaded ' + variable_name + ' data')
    
        tmean_highres1 = lapse_rate_based_downscaling(dem_high, tlapse_mean, dem_low, tmean)
        tmean_highres = grid_calculator(tmean_highres1,landseamask,'a*b')
        tmean = grid_calculatorX(oceans,tmean,'a*b')
        tmean = closegaps(tmean)
        tmean_highres = patchingBspline(tmean_highres,tmean)
    
        tmean_srad_ccINT = convert2uinteger10(tmean_highres, variable_long_name[variable_name])
        export_geotiff(tmean_srad_ccINT, variable_name)
    
        # clean up memory
        saga_api.SG_Get_Data_Manager().Delete(tmean_srad_ccINT)
        saga_api.SG_Get_Data_Manager().Delete(tmean_highres)
        saga_api.SG_Get_Data_Manager().Delete(tmean_highres1)
        track_memory_usage(variable_name + ' cleaned up memory')

    # clean up memory
    saga_api.SG_Get_Data_Manager().Delete(dem_high)
    track_memory_usage('downscale temperature end')


    # Solar radiation **************************************************************************************************
    if calc_rsds:
        srad_tot = import_gdal(SRAD + 'CHELSA_stot_pj_' + dayofyear + '_V.2.1.tif')

        # import srad quotient and cloudcover
        ccbiascor = import_gdal(CCOVER + 'CHELSA_tcc_' + day + '_' + month + '_' +  year + '_V.2.1.tif')
        ccbiascor = grid_calculator_simple(ccbiascor, 'a/10000')
        srad_sur = surface_radiation(srad_tot, ccbiascor, variable_long_name['rsds'])
        track_memory_usage('srad loaded start')
    
        ## bias correct surface solar radiation downwards
        ssrd = import_gdal(TEMP + 'ssrd.tif')
    
        srad_resamp = resample(srad_sur, ssrd)
        # convert Jm^2 to Wm^2
        srad_resamp = grid_calculator_simple(srad_resamp, 'a*0.277778/24')
    
        srad_bias = grid_calculator(ssrd,srad_resamp,'(a+100)/(b+100)')
        srad_cor  = grid_calculatorX(srad_sur, srad_bias,'a*b')
        srad_cor = change_data_storage2(srad_cor, 3)
        track_memory_usage('srad loaded end')
    
        # export the bias corrected solar radiation
        export_geotiff(srad_cor, 'rsds')
        track_memory_usage('srad_biascor done')


    saga_api.SG_Get_Data_Manager().Delete_All()
    track_memory_usage('RUN FINISHED')

    # remove the temporary directory
    shutil.rmtree(TEMP)
