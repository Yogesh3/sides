import numpy as np
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
import pickle

from IPython import embed


for cii_mod in ['de_Looze', 'Lagache']:

    hdu = fits.open('./CUBES_3D/cube_3D_z6_MJy_sr_'+cii_mod+'.fits')
    hdr = hdu[0].header
    cube = hdu[0].data
    
    pow_sqr = np.absolute(np.fft.fftn(cube)**2 * hdr['CDELT1'] * hdr['CDELT2'] *hdr['CDELT3'] / (hdr['NAXIS1'] * hdr['NAXIS2'] * hdr['NAXIS3']) )
    
    u_freq = np.fft.fftfreq(hdr['NAXIS1'], d=hdr['CDELT1'])
    v_freq = np.fft.fftfreq(hdr['NAXIS2'], d=hdr['CDELT2'])
    w_freq = np.fft.fftfreq(hdr['NAXIS3'], d=hdr['CDELT3'])
    
    kfreq = np.sqrt(u_freq[:,np.newaxis,np.newaxis]**2 + v_freq[np.newaxis,:,np.newaxis]**2 + w_freq[np.newaxis,np.newaxis,:]**2)
    
    kbins = np.linspace(np.min(kfreq), np.max(kfreq), num = int(hdr['NAXIS2'] / 2))
    
    norm, bin_edges = np.histogram(kfreq, bins=kbins, range=range)
    hist, bin_edges = np.histogram(kfreq, bins=kbins, weights=pow_sqr)
    
    Pk = hist / norm * (u.MJy / u.sr)**2 * u.Mpc**3
    kbins = kbins / u.Mpc
    
    kmean = 0.5*(kbins[0:-1] + kbins[1:])
    
    dict = {'kmean': kmean[1:], 'Pk': Pk[1:].to('Jy2 Mpc3 / sr2')}
    pickle.dump(dict, open('Pk_3D_z6_'+cii_mod+'.p', 'wb'))
    
    #plt.plot(kmean, 1.e12*Pk)
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.show()

    hdu.close()


embed()
