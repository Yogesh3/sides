import numpy as np
import astropy as ap
import pickle
from IPython import embed
from astropy.cosmology import Planck15 as cosmo

cat = pickle.load(open('../PYSIDES/OUTPUTS/pySIDES_from_original.p', 'rb'))

#Useful quantities#
fieldsize = 2. * (np.pi/180.)**2 #size of SIDES
nuCO = 115.27 #GHz
Jupperlist = [1,2,3,4,3,4,5,6,7,8,1] # the last bin is Riechers+18 from coldz
zcenlist = [0.27,1.40,2.56,3.73,0.49,0.95,1.43,1.91,2.39,2.87,2.4]
ref_2compare = ['ASPECS3mm','ASPECS3mm','ASPECS3mm','ASPECS3mm','ASPECS1mm','ASPECS1mm','ASPECS1mm','ASPECS1mm','ASPECS1mm','ASPECS1mm','ColdZ']
delta_z = 0.1

logLmin = 5.
logLmax = 12.

Nbin = 28

LF_array = np.zeros((Nbin, np.size(Jupperlist)))

Deltabin = (logLmax - logLmin) * 1. / Nbin #in dex

logLmean = logLmin + Deltabin * ( np.arange(Nbin) + 0.5 )

for k in range(0, np.size(Jupperlist)):

    inzbin = np.where(np.abs(cat['redshift'] - zcenlist[k]) < delta_z)
    ICO_inzbin = cat['ICO{:d}{:d}'.format(Jupperlist[k], Jupperlist[k]-1)][inzbin[0]]
    logLprim_inzbin = np.log10(3.25e7 * ICO_inzbin * cat['Dlum'][inzbin[0]]**2 / (1. + cat['redshift'][inzbin[0]])**3 / (nuCO * Jupperlist[k] / (1. + cat['redshift'][inzbin[0]]))**2 )
    
    histo = np.histogram(logLprim_inzbin, bins = Nbin, range = (logLmin, logLmax))



    Vol_zbin = 1./ 3. * (cosmo.comoving_distance(zcenlist[k]+delta_z)**3 - cosmo.comoving_distance(zcenlist[k]-delta_z)**3) * fieldsize # Mpc^3

    LF_array[:,k] = histo[0] / Deltabin / Vol_zbin

pickle.dump((LF_array, logLmean, Jupperlist, zcenlist, ref_2compare), open('LFarr_SIDES.p','wb'))


