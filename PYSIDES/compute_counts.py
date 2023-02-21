import pickle
import numpy as np
from IPython import embed
from pysides.gen_fluxes import *

Omega_sides = 2. # deg^2

def compute_counts(wl):

    #lambda in microns

    cat = pickle.load(open("pySIDES_from_original.p", "rb"))

    #for debugging
    cat = cat.loc[cat['SFR']>0]
    cat.reset_index(inplace = True)
    
    if 'S{:0.0f}'.format(wl) in cat.columns:
        print('wavelength is already standard!')
    else:
        cat = gen_fluxes(cat, lambda_filters = [wl])

    S_sides = cat['S{:0.0f}'.format(wl)]

    Scuts_sides = 10.**np.arange(-8.,1.,0.05)#compute SIDES number counts

    NsupS_sides = []
    eNsupS_sides = []

    print('compute the number counts')

    #for Scut in Scuts_sides:
    #    N = np.sum(S_sides >= Scut)
    #    NsupS_sides.append(N)
    #    eNsupS_sides.append(np.sqrt(N))
    #NsupS_sides = np.array(NsupS_sides) / Omega_sides #to go to deg^-2
    #eNsupS_sides = np.array(eNsupS_sides) / Omega_sides #to go to deg^-2


    Sdiff_sides = []
    for k in range(0,len(Scuts_sides)-1):
        Sdiff_sides.append(np.mean([Scuts_sides[k], Scuts_sides[k+1]]))
    Sdiff_sides = np.array(Sdiff_sides)
               
    diff_sides = []
    ediff_sides = []
    for i in range(0,len(Scuts_sides)-1):
        N = np.sum( (S_sides >= Scuts_sides[i]) &  (S_sides < Scuts_sides[i+1]))
        diff_sides.append(N / (Scuts_sides[i+1] - Scuts_sides[i]))
        ediff_sides.append(np.sqrt(N) / (Scuts_sides[i+1] - Scuts_sides[i]))
        
    diff_sides = np.array(diff_sides) * (Sdiff_sides)**2.5 * (180. / np.pi)**2 / Omega_sides
    ediff_sides = np.array(ediff_sides) * (Sdiff_sides)**2.5 * (180. / np.pi)**2 / Omega_sides

    str = '#SIDES number counts at {:0.0f} microns\n'.format(wl)
    str += '# column 1: Flux density in Jy\n'
    str += '# column 2: Euclidian-normalized differential counts dN/dS*S^2.5 in gal*Jy^1.5/sr\n'
    str += '# column 3: uncertainty\n'

    for k in range(0, len(Sdiff_sides)):
        str += '{} {} {}\n'.format(Sdiff_sides[k], diff_sides[k],  ediff_sides[k])

    txt = open('COUNTS_TABLES/SIDES_counts_{:0.0f}.txt'.format(wl), 'w')

    txt.write(str)

    txt.close()


    embed()

    return 0

compute_counts(1200)
