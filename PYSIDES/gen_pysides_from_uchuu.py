from pysides.load_params import *
from pysides.load_uchuu_cats import *
from pysides.gen_sfr_props import *
from pysides.gen_magnification import *
from pysides.gen_fluxes import *
from pysides.gen_fluxes_filter import *
from pysides.gen_lines import *
from pysides.gen_outputs import *
import os

#loop that creates the cubes for all the UCHUU catalogs
uchuu_cats_folder_path='/data/SIDES/UCHUU/UCHUUcatalogs/vpeak/'
with os.scandir(uchuu_cats_folder_path) as it:
    for entry in it:

        print('Creation of the cubes for the catalog {}'.format(entry.name))
        
        uchuu_cat_path = uchuu_cats_folder_path + entry.name

        params = load_params('PAR_FILES/SIDES_from_uchuu.par')
        
        tile0=(entry.name.split('.')[0]).split('_')[-2]
        tile1=(entry.name.split('.')[0]).split('_')[-1]

        params['run_name']='pySIDES_from_uchuu_tile_'+tile0+'_'+tile1

        cat = load_uchuu_cats(uchuu_cat_path)

        cat = gen_sfr_props(cat, params) #defines the galaxy type and generates the SFR

        cat = gen_magnification(cat, params, magnify = True) #picks the mu value according to the galaxy's z and the corresponding probability distribution (taken from Hilbert et al.2007)

        cat = gen_fluxes(cat, params) #calculates the fluxes

        #cat = gen_fluxes_filter(cat, params) #generates the final template's flux (lambda in filter bands)

        cat = gen_CO(cat, params) # inserting CO lines in the continuum of the galaxies

        cat = gen_CII(cat, params) # inserting CII lines in the continuum of the galaxies

        cat = gen_CI(cat, params) # inserting CI lines in the continuum of the galaxies

        gen_outputs(cat, params) #gives back all the results

