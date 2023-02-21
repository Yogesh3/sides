from pysides.make_cube import * 
from pysides.load_params import *
from pysides.load_sides_csv import *
from pysides.gen_sfr_props import *
from pysides.gen_magnification import *
from pysides.gen_fluxes import * 
from pysides.gen_lines import *


params_sides = load_params('PAR_FILES/SIDES_from_original.par')

params_cube = load_params("PAR_FILES/CONCERTO_cubes.par")

cat = Table.read(params_cube["sides_cat_path"]) 
cat = cat.to_pandas()

make_cube(cat ,params_sides, params_cube)
