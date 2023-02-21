import numpy as np
import astropy as ap
import pickle
from IPython import embed
from astropy.cosmology import Planck15 as cosmo

cat = pickle.load(open('../PYSIDES/OUTPUTS/pySIDES_from_original.p', 'rb'))

#Useful quantities#
fieldsize = 2. * (np.pi/180.)**2 #size of SIDES
zcenlist = [1.08,2.40]
transition = ['10', '21']
nu_rest = [492.16, 809.34]

delta_z = 0.2

logLmin = 5.
logLmax = 12.

Nbin = 28

LF_array = np.zeros((Nbin, np.size(zcenlist)))

Deltabin = (logLmax - logLmin) * 1. / Nbin #in dex

logLmean = logLmin + Deltabin * ( np.arange(Nbin) + 0.5 )

for k in range(0, np.size(zcenlist)):

    inzbin = np.where(np.abs(cat['redshift'] - zcenlist[k]) < delta_z)

    Vol_zbin = 1./ 3. * (cosmo.comoving_distance(zcenlist[k]+delta_z)**3 - cosmo.comoving_distance(zcenlist[k]-delta_z)**3) * fieldsize # Mpc^3
                    
    logLCI_inzbin = np.log10(cat['ICI'+transition[k]][inzbin[0]] *  (1.04e-3 * cat['Dlum'][inzbin[0]]**2 * nu_rest[k] / (1 + cat['redshift'][inzbin[0]])))

    
    histo = np.histogram(logLCI_inzbin, bins = Nbin, range = (logLmin, logLmax))
    
    LF_array[:,k] = histo[0] / Deltabin / Vol_zbin

pickle.dump((LF_array, logLmean, zcenlist, transition, nu_rest), open('LFCIarr_SIDES.p','wb'))
