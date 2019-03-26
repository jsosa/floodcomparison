#!/usr/bin/env python

# inst: university of bristol
# auth: jeison sosa
# mail: sosa.jeison@gmail.com / j.sosa@bristol.ac.uk

import os
import sys
import numpy as np
import gdalutils as gu
from subprocess import call

def floodcomparison(obsf,modf,thresh,outfolder,buffer=0.01):

    """
    Compare two flood extent masks. Calc Hit rate,
    False alarm ratio, Critical success index and Bias.
    A threshold is applied to the `modf` (assumed to be a 
    water depth) file to generate a mask file where 1 is
    wet and 0 is dry (mask.tif). A buffer of 0.01 degrees 
    is applied to generate a buffer around the `obsf` file
    (buff.tif). Agreement between observed and modeled flood
    extents is calculated and saved as (aggr.tif)
    """

    # Create output folder
    try:
        os.mkdir(outfolder)
    except FileExistsError:
        pass

    # Reading input files
    obs = gu.get_data(obsf)
    geo1 = gu.get_geo(obsf)
    mod = gu.get_data(modf)
    geo2 = gu.get_geo(modf)

    # Grid extent in both files should have same dimensions
    if (geo1[4]==geo2[4]) & (geo1[5]==geo2[5]):

        # Calc fist intersection between obs and mod
        # Sometimes obs or mod is not available for some tribs in both
        # After intersectin is calculated a buffer can be applied
        tmp_mod = np.where(mod>thresh,1,0)
        tmp_obs = np.where(obs>0,1,0)
        tmp_buff_msk = np.where((tmp_mod+tmp_obs)==2,1,0)
        gu.write_raster(tmp_buff_msk,outfolder+'buffer_mask.tif',geo2,'Int16',0)

        # Calc proximity around `buffer_mask`
        call(['gdal_proximity.py','-distunits','GEO',
                                  '-co','COMPRESS=LZW',
                                  '-maxdist',str(buffer),
                                  '-nodata','-9999',
                                  outfolder+'buffer_mask.tif',
                                  outfolder+'buffer_dist.tif'])
        
        # Save previous as a mask (ones and zeros)
        buff = np.where(gu.get_data(outfolder+'buffer_dist.tif')>=0,1,0)
        gu.write_raster(buff,outfolder+'buff.tif',geo2,'Int16',0)
        
        # Clip `modf` with the buffer previously calculated
        nmod = np.where(buff==1,mod,0)
        nobs = np.where(buff==1,obs,0)

        # Create a mask of `nmod`
        mask = np.where(nmod>thresh,1,0)
        gu.write_raster(mask,outfolder+'mask.tif',geo2,'Int16',0)
        gu.write_raster(nobs,outfolder+'obs.tif',geo2,'Int16',0)

        # Calculate aggrement between observation and model within the buffer
        aggr = np.where((nobs+mask)==2,1,0)
        gu.write_raster(aggr,outfolder+'aggr.tif',geo2,'Int16',0)

        # Creating mask of wet and dry areas in both observation and
        # model within the buffer
        obs_wet = nobs
        mod_wet = mask
        obs_dry_tmp = np.where(obs_wet==1,0,1)
        mod_dry_tmp = np.where(mod_wet==1,0,1)
        obs_dry = np.where(buff==1,obs_dry_tmp,0)
        mod_dry = np.where(buff==1,mod_dry_tmp,0)

        # Calculating intersection
        TP = np.where((obs_wet+mod_wet)==2,1,0)
        TN = np.where((obs_wet+mod_dry)==2,1,0)
        FP = np.where((obs_dry+mod_wet)==2,1,0)
        FN = np.where((obs_dry+mod_dry)==2,1,0)

        # Calc parameters for scores
        a = TP.sum()
        b = FP.sum()
        c = TN.sum()
        d = FN.sum()

        # Calc scores
        # Bias : zero - best, positive - overprediction, negative - underprediction
        H = a/(a+c)
        F = b/(a+b)
        C = a/(a+b+c)
        B = (a+b)/(a+c)-1

        with open(outfolder + 'scores.txt','w') as f:
            f.write('Hit rate: ' + '%3.2f'%H + '\n')
            f.write('False alarm ratio: ' + '%3.2f'%F + '\n')
            f.write('Critical success index: ' + '%3.2f'%C + '\n')
            f.write('Bias: ' + '%3.2f'%B + '\n')

        return H,F,C,B
    
    else:
        sys.exit('Input files have different dimmensions')
