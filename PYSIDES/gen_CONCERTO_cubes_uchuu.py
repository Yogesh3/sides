from pysides.make_cube import * 
from pysides.load_params import *
from pysides.load_sides_csv import *
from pysides.gen_sfr_props import *
from pysides.gen_magnification import *
from pysides.gen_fluxes import * 
from pysides.gen_lines import *
import os


params_sides = load_params('PAR_FILES/SIDES_from_uchuu.par')
params_cube = load_params("PAR_FILES/CONCERTO_cubes_uchuu.par")

dirpath = "/data/SIDES/PYSIDES_UCHUU_OUTPUTS/vpeak/"
with os.scandir(dirpath) as it:
    for entry in it:
        if entry.is_file():
            tile0=(entry.name.split('.')[0]).split('_')[-2]
            tile1=(entry.name.split('.')[0]).split('_')[-1]
            #params_cube['sides_cat_path'] = dirpath + "pySIDES_from_uchuu_tile_"+tile0+"_"+tile1+".fits"
            params_cube['sides_cat_path'] = dirpath + entry.name
            print('Creating cubes for catalog ' + entry.name)
            params_cube['run_name'] = "pySIDES_from_uchuu_"+tile0+"_"+tile1+"_CONCERTO"

            cat = Table.read(params_cube["sides_cat_path"]) 
            cat = cat.to_pandas()
            #cat = cat.loc[(cat.LCII_de_Looze < 1e+8)]
            #cat.reset_index(inplace = True)

            make_cube(cat ,params_sides, params_cube)
            print('Cubes under the name ' + params_cube['run_name'] + ' were created')
