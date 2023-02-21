from pysides.make_cube import * 
from pysides.load_params import *
from pysides.load_sides_csv import *
from pysides.gen_sfr_props import *
from pysides.gen_magnification import *
from pysides.gen_fluxes import * 
from pysides.gen_lines import *

from IPython import embed
import numpy as np

params_sides = load_params('PAR_FILES/SIDES_from_original.par')

params_cube = load_params("PAR_FILES/CONCERTO_cubes.par")

params_cube['run_name'] = 'Test_0p09deg2_astrometry'

params_cube['save_each_transition']  = False    #save the cubes for each tranistion (CO10, CO21, [CI]10, [CI]21...)
params_cube['save_each_line']        = False   #save the cubes for each line (CO, [CI], [CII])
params_cube['save_continuum_only']   = False    #save the cube(s) with only continuum
params_cube['save_all_lines']        = False    #save the cube(s) with all the lines
params_cube['save_full']             = True    #save the full cube(s) with everything
params_cube['gen_cube_CII_Lagache']  = False

params_cube["sides_cat_path"] = '/Users/mbethermin/SIDES/PYSIDES/OUTPUTS/pySIDES_from_original.fits'
params_cube["output_path"] = '/Users/mbethermin/SIDES/PYSIDES/OUTPUTS/CUBES/'

cat = Table.read(params_cube["sides_cat_path"]) 
cat = cat.to_pandas()

cat = cat.loc[(np.abs(cat['ra']) < 0.3) & (np.abs(cat['dec']) < 0.3)]
cat.reset_index(inplace = True)

#Put a supidely bright sources to test calibration
cat['ICO21'][0] = 1.e6
cat['ra'][0] = 0.1 + 5./3600.
cat['dec'][0] = 0.1 + 5./3600.
cat['redshift'][0] = 0.

cat['ICO21'][1] = 2.e6
cat['ra'][1] = 0.1
cat['dec'][1] = 0.2
cat['redshift'][1] = 0.

cat['ICO21'][2] = 3.e6
cat['ra'][2] = 0.2
cat['dec'][2] = 0.1
cat['redshift'][2] = 0.

make_cube(cat ,params_sides, params_cube)
