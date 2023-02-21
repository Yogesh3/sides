import numpy as np
import astropy as ap
import pickle
from IPython import embed
from astropy.cosmology import Planck15 as cosmo

cat = pickle.load(open('../PYSIDES/OUTPUTS/pySIDES_from_original.p', 'rb'))

#Useful quantities#
fieldsize = 2. * (np.pi/180.)**2 #size of SIDES
zcenlist = [4.5,5.5, 6.0]
delta_z = [0.1,0.4,0.2]

method = ['Lagache', 'de_Looze']

logLmin = 5.
logLmax = 12.

Nbin = 28

nu_rest = 1900.54

LF_array = np.zeros((Nbin, np.size(zcenlist), len(method)))

Deltabin = (logLmax - logLmin) * 1. / Nbin #in dex

logLmean = logLmin + Deltabin * ( np.arange(Nbin) + 0.5 )

for k in range(0, np.size(zcenlist)):

    inzbin = np.where(np.abs(cat['redshift'] - zcenlist[k]) < delta_z[k])

    Vol_zbin = 1./ 3. * (cosmo.comoving_distance(zcenlist[k]+delta_z[k])**3 - cosmo.comoving_distance(zcenlist[k]-delta_z[k])**3) * fieldsize # Mpc^3
                    
    for i in range(0, np.size(method)):
        logLCII_inzbin = np.log10(cat['ICII_'+method[i]][inzbin[0]] *  (1.04e-3 * cat['Dlum'][inzbin[0]]**2 * nu_rest / (1 + cat['redshift'][inzbin[0]]))) #np.log10(cat['LCII_'+method[i]][inzbin[0]])
    
        histo = np.histogram(logLCII_inzbin, bins = Nbin, range = (logLmin, logLmax))

        LF_array[:,k,i] = histo[0] / Deltabin / Vol_zbin

pickle.dump((LF_array, logLmean, zcenlist, method), open('LFCIIarr_SIDES.p','wb'))
