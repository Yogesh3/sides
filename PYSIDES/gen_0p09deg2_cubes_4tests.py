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

params_cube['run_name'] = 'Test_0p09deg2'

cat = Table.read(params_cube["sides_cat_path"]) 
cat = cat.to_pandas()

cat = cat.loc[(np.abs(cat['ra']) < 0.3) & (np.abs(cat['dec']) < 0.3)]
cat.reset_index(inplace = True)

#Put a supidely bright source to test calibration
#cat['ICO21'][0] = 100.

make_cube(cat ,params_sides, params_cube)
