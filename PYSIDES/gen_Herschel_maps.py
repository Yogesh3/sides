from pysides.make_maps import *
from pysides.load_params import *
from astropy.table import Table

params_sides = load_params('PAR_FILES/SIDES_from_original.par')

params_maps = load_params("PAR_FILES/Herschel_maps.par")

cat = Table.read(params_maps["sides_cat_path"])
cat = cat.to_pandas()

make_maps(cat, params_maps, params_sides)
