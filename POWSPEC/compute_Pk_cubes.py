import numpy as np
from powspec.powspec import power_spectral_density
from astropy.io import fits
import matplotlib.pyplot as plt
import astropy.units as u
import pickle

from IPython import embed

todo = ['continuum','all_lines_Lagache', 'all_lines_de_Looze', 'CO_all', 'CO10', 'CO21', 'CO32', 'CO43', 'CO54', 'CO65', 'CO76', 'CO87', 'CI21',  'CI10',  'CII_Lagache',  'CII_de_Looze'] 

for i in [0, 1]:

    if i == 0:
        cube_path = '/data/SIDES/PYSIDES_ORIGINAL_OUTPUTS/CUBES/pySIDES_from_original_CONCERTO_'
    if i == 1:
        cube_path = '/data/SIDES/PYSIDES_ORIGINAL_OUTPUTS/CUBES/MASKED/pySIDES_from_original_CONCERTO_masked_'

    cube_ext = '_nobeam_MJy_sr.fits'

    Pk_dict = {}

    Pk_dict['cube_types'] = todo

    for cube_type in todo:
        
        print('Compute P(k) of the '+cube_type+' cube...')
        
        path = cube_path + cube_type + cube_ext
        
        Pk_arr = []
        
        with fits.open(path) as cube:
            
            res = cube[0].header['CDELT1'] * u.deg
            
            Nslice = np.shape(cube[0])[0]
            
            for k in range(0, Nslice):
                pk , kbins = power_spectral_density(cube[0].data[k,:,:] * u.MJy / u.sr , res = res.to('arcmin'))
                pk = pk.to('Jy^2/sr')
                
                Pk_arr.append(pk)
                
                if not 'kbins' in Pk_dict.keys():
                    Pk_dict['kbins'] = kbins
                else:
                    if np.sum(np.abs(Pk_dict['kbins']- kbins)).value > 1.e-9:
                        print('Error!!! The bins of each cubes do not have the same size. All the cubes should have the same resolution.')
                        exit

        unit_Pk = pk.unit
        Pk_arr = np.array(Pk_arr) * unit_Pk
        Pk_dict['Pk_arr_'+cube_type] = Pk_arr
    

    Pk_dict['kmean'] = kmean = 0.5 * (Pk_dict['kbins'][0:-1] + Pk_dict['kbins'][1:])

    if i ==0:
        pickle_file_name = 'Pk_results.p'
    if i == 1:
        pickle_file_name = 'Pk_results_masked.p'

    print('Save the results in ', pickle_file_name, '...')

    pickle.dump(Pk_dict,  open(pickle_file_name, 'wb'))

embed()
