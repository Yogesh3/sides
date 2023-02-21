from pysides.make_maps import *
from pysides.load_params import *
from astropy.table import Table

params_sides = load_params('PAR_FILES/SIDES_from_original.par')

params_maps = load_params("PAR_FILES/Herschel_maps.par")

params_maps["sides_cat_path"] = '/Users/mbethermin/SIDES/PYSIDES/OUTPUTS/pySIDES_from_original.fits'
params_maps["output_path"] = '/Users/mbethermin/SIDES/PYSIDES/OUTPUTS/MAPS/'
params_maps["run_name"] = 'pySIDES_Herschel_test_astrometry'

params_maps["filter_list"]  = ['SPIRE250']
params_maps["beam_fwhm_list"] = [18.2]              
params_maps["pixel_size"]  = [6.] 

cat = Table.read(params_maps["sides_cat_path"])
cat = cat.to_pandas()

cat['SSPIRE250'][0] = 1.e3
cat['ra'][0] = 0.1
cat['dec'][0] = 0.2
cat['redshift'][0] = 0.

cat['SSPIRE250'][1] = 2.e3
cat['ra'][1] = 1.2
cat['dec'][1] = 0.3
cat['redshift'][1] = 0.

cat['SSPIRE250'][2] = 3.e3
cat['ra'][2] = 0.2 #+3/3600 (to offset by half a pixel to test half pixel problems => was also run and still OK!)
cat['dec'][2] = 1.30 # +3/3600
cat['redshift'][2] = 0.

make_maps(cat, params_maps, params_sides)
