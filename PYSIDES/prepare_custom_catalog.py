import pickle
from IPython import embed
from pysides.gen_fluxes import *
from pysides.load_params import *
from pysides.gen_fluxes_filter import *
from pysides.gen_outputs import *

cat = pickle.load(open('OUTPUTS/pySIDES_from_original.p', 'rb'))
params = load_params('PAR_FILES/SIDES_from_original.par', force_pysides_path = '/Users/mbethermin/SIDES/PYSIDES/')

#IMPORTANT MODIFY THE RUN NAME IF YOU DO NOT WANT TO ERASE THE ORIGINAL FILE#
params['run_name'] = 'pySIDES_custom'

#Add new filters
params['filter_list'] = ['LABOCA870'] #modify the lits of filters compared to the default parameters
cat = gen_fluxes_filter(cat, params)

#Add new monochromatic fluxes
new_lambda = [1130] #list of new monochramatic fluxes in microns
cat = add_fluxes(cat, params, new_lambda)


#Cut the catalog (useful if people are working on detected galaxy samples to reduce file size)
cat = cat.loc[cat['S1130'] > 2.e-5] #you can also cut before if the quantity used for the selections are already generated (more efficient)

#Save the catalog
cat = cat[['ra', 'dec', 'redshift', 'S1130']] #if you want to limit the number of columns to export
params['gen_pickle'] = False #no need to create a pickle if you want to send it to a colleague later, else comment
gen_outputs(cat, params)

embed()
