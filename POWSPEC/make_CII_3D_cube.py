import numpy as np
import pickle
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import z_at_value
import astropy.units as u
import scipy.constants as cst
from astropy.io import fits
import datetime

from IPython import embed

z_center = 6
Dc_center = cosmo.comoving_distance(z_center)
Nvox = 512 #number of pixel in each dimension

outfile = './CUBES_3D/cube_3D_z6_MJy_sr_'

cat = pickle.load(open('/Users/mbethermin/SIDES/PYSIDES/OUTPUTS/pySIDES_from_original.p', 'rb'))

ra_center = np.mean(cat['ra'])
dec_center = np.mean(cat['dec'])

delta_ra = np.max(cat['ra']) - np.min(cat['ra'])
delta_dec = np.max(cat['dec']) - np.min(cat['dec'])

#all the coordinates will be in comoving units
size_x = Dc_center.value * delta_ra * (np.pi/180) * np.cos(np.pi/180*dec_center)
size_y = Dc_center.value * delta_dec * (np.pi/180) 
box_size = np.min([size_x, size_y])

#compute the coordinates in the cube
ys = Dc_center.value * (cat['ra'] - ra_center) * (np.pi/180) * np.cos(np.pi/180*cat['dec'])
xs = Dc_center.value * (cat['dec'] - dec_center)* (np.pi/180)

bins = np.linspace(-box_size/2, box_size/2, num = Nvox + 1)

vox_size = bins[1]-bins[0]

z_bins = []
for Dc in Dc_center + bins * u.Mpc:
    z_bins.append(z_at_value(cosmo.comoving_distance, Dc))

for cii_mod in ['de_Looze', 'Lagache']:
 
    cube_MJy_per_sr, edges = np.histogramdd(sample = (cat['redshift'], ys, xs), bins = (z_bins, bins, bins), weights = cat['ICII_'+cii_mod])

    #convert to the proper unit (MJy/sr)
    for k in range(0,Nvox):
        area_voxel = (vox_size / Dc_center.value)**2
        dv_voxel = (z_bins[k+1]-z_bins[k]) / (1 + 0.5*(z_bins[k+1]+z_bins[k])) * cst.c * 1.e-3 #km/s
        cube_MJy_per_sr[k,:,:] *=  1.e-6 / area_voxel / dv_voxel
    
    #save the cube!
    f= fits.PrimaryHDU(cube_MJy_per_sr)
    hdu = fits.HDUList([f])
    hdr = hdu[0].header
    hdr["BUNIT"] = 'MJy/sr'
    for k in range(1,4):
        hdr["CDELT{}".format(k)] = vox_size
    hdr["DATE"] = (str(datetime.datetime.now()), "date of the creation")
    hdu.writeto(outfile+cii_mod+'.fits', overwrite=True)
    hdu.close()

embed()

