from pysides.load_params import *
from pysides.load_sides_csv import *
from pysides.gen_sfr_props import *
from pysides.gen_magnification import *
from pysides.gen_fluxes import *
from pysides.gen_fluxes_filter import *
from pysides.gen_lines import *
from pysides.gen_outputs import *

csv_idl_path = '/data/SIDES/CESAM/SIDES_Bethermin2017_corr.csv'

params = load_params('PAR_FILES/SIDES_from_original.par')

cat = load_sides_csv(csv_idl_path)

cat = gen_sfr_props(cat, params)

cat = gen_magnification(cat, params)

cat = gen_fluxes(cat, params)

cat = gen_fluxes_filter(cat, params)

cat = gen_CO(cat, params)

cat = gen_CII(cat, params)

cat = gen_CI(cat, params)

gen_outputs(cat, params)

