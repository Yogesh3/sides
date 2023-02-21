from pysides.make_cube import * 
from pysides.load_params import *
from pysides.load_sides_csv import *
from pysides.gen_sfr_props import *
from pysides.gen_magnification import *
from pysides.gen_fluxes import * 
from pysides.gen_lines import *

from IPython import embed

params_sides = load_params('PAR_FILES/SIDES_from_original.par')

params_cube = load_params("PAR_FILES/CONCERTO_cubes.par")

#change some parameters compared to the normal run
params_cube['output_path'] += 'MASKED/'
params_cube['save_cube_nobeam_Jy_pix'] = False
params_cube['gen_cube_smoothed_Jy_beam'] = False
params_cube['gen_cube_smoothed_MJy_sr'] = False
params_cube['run_name'] += '_masked' 

cat = Table.read(params_cube["sides_cat_path"]) 
cat = cat.to_pandas()

#Keep only in the catalog the unmasked (faint) sources. The mass cut is based on the mass limit of the Laigle et al. catalog (full sample, deep). 
zbins = [0., 0.35, 0.65, 0.95, 1.3, 1.75, 2.25, 2.75, 3.5, 4.0, 4.8]
logMlim = [8.1, 8.7, 9.1, 9.3, 9.7, 9.9, 10.0, 10.1, 10.1, 10.1]

keep = np.ones(len(cat), dtype = bool)

for k in range(0, len(logMlim)):
    sel = np.where( (cat['redshift'] >= zbins[k]) & (cat['redshift'] < zbins[k+1]) )
    keep[sel[0]] = (cat['Mstar'][sel[0]] < 10.**logMlim[k])

cat = cat.loc[keep == True]
cat.reset_index(inplace = True)

make_cube(cat ,params_sides, params_cube)
