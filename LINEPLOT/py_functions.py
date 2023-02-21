import numpy as np
import matplotlib.pyplot as plt
import vaex as vx
from astropy.cosmology import Planck18_arXiv_v2 as cosmo
import astropy.units as u
import os
import sys


def gen_lf(Lline, field_size, z1, z2, bins = 28, min_logLbin_value = 5, max_logLbin_value = 12):
    
    '''
    Function that constructs the luminosity function given:
    - L_line: the line (e.g., CO(J - J-1), [CII], [CI])
    - field_size: the size of the field (survey size)
    - z1, z2: the redshift slice
    - bins: the number of luminosity bins (default = 28)
    - min_bin_value, max_bin_value: the range of the luminosities (default = [1e+5, 1e+12]*units)
    OUTPUT: the luminosity bins (L) and the Φ(L) value for each L
    '''
    
    counts = Lline.count(binby = np.log10(Lline), limits = [min_logLbin_value, max_logLbin_value], shape = bins)
    omega = (np.pi / 180) ** 2 * field_size
    Vz = omega / 3 * ((cosmo.luminosity_distance(z2) / (1+z2)) ** 3 - (cosmo.luminosity_distance(z1) / (1+z1)) ** 3)
    Lbins = 10**np.linspace(min_logLbin_value, max_logLbin_value, bins+1)
    Lbins_center = 0.5 * (Lbins[1:] + Lbins[:-1])
    dex = np.diff(np.log10(Lbins))[0]
    phi = counts / Vz / dex
    
    return Lbins_center, phi.value

def load_catalogs(SIDES_cats_path, field_size, z1, z2, tile):

    '''
    Given the surey size (field_size) and the redshift slice loads the corresponding part of the total Uchuu field
    ATTENTION selecting a field_size above 1 sq.deg will force the loading of the total dataframe which will take around 10 min to run
    '''
    
    names=[]
    with os.scandir(SIDES_cats_path) as it:
        for entry in it:
            if entry.is_file():
                names.append(SIDES_cats_path + entry.name)
    filenames = sorted(names)
    
    omega = (np.pi / 180) ** 2 * field_size
    dV = (omega / 3 * ((cosmo.luminosity_distance(z2) / (1 + z2)) ** 3 - (cosmo.luminosity_distance(z1) / (1 + z1)) ** 3)).value

    if field_size <= 1:
        bdf = vx.open(filenames[tile])
        ramin, decmin = bdf.ra.min(), bdf.dec.min()
        radec_filter = (bdf.ra < ramin + np.sqrt(field_size)) & (bdf.dec < decmin + np.sqrt(field_size))
        z_filter = (bdf.redshift > z1) & (bdf.redshift <= z2)
        df = bdf[z_filter & radec_filter]
        df.extract()
    else:
        print('You have selected a field size > 1!!! Loading the whole data frame of 117 sq.deg will take some time (~ 10 mins) ...')
        bdf = vx.open_many(filenames)
        z_filter = (bdf.redshift > z1) & (bdf.redshift <= z2)
        df = bdf[z_filter]
        df.extract()

    return df, dV

def line_computation(df, line, J):
    
    '''
    Given the line and the transition computes the luminosity of all the galaxies in the catalog
    Possible options: CO(1-0), CO(2-1), CO(3-2), CO(4-3), CO(5-4), CO(6-5), CO(7-6), CO(8-7), [CII], [CI](1-0), [CI](2-1)
    '''
    
    if line == 'CO':
        nu_rest = np.array([115.2712, 230.5380, 345.7960, 461.0408, 576.2679, 691.4731, 806.652, 921.7997])
        df['LprimCO{}{}'.format(J,J-1)] = 3.25e7 * df['ICO{}{}'.format(J,J-1)] * df['Dlum']**2 / (1. + df['redshift'])**3 * ((1. + df['redshift']) / (nu_rest[J-1]))**2
        Lline = df['LprimCO{}{}'.format(J,J-1)]
        units = 'K \, km \, m^{-1}pc^2'
        transition = '{}{}'.format(J, J-1)
        
    elif line == 'CII':
        Lline = df['LCII_de_Looze']
        units = 'L_{\odot}'
        transition = ' '
        
    else:
        nu_rest = np.array([492.16, 809.344])
        df['LCI{}{}'.format(J,J-1)] = 3.25e7 * df['ICI{}{}'.format(J,J-1)] * df['Dlum']**2 / (1. + df['redshift'])**3 * ((1. + df['redshift']) / (nu_rest[J-1]))**2
        Lline = df['LCI{}{}'.format(J,J-1)]
        units = 'L_{\odot}'
        transition = '{}{}'.format(J, J-1)
        
    return Lline, units, transition



def LF_computation(SIDES_cats_path, line, J, field_size, z1, z2, Lbins, min_logLbin_value, max_logLbin_value, tile, variance = False, plot = True):
    
    '''
    Computes the luminosity function of the requested line as well as the field-to-field variance if selected
    -----
    INPUT:
    SIDES_cats_path: path where the SIDES-Uchuu catalogs are stored
    line, J: the line and the transition (if CO or [CI]) you aim to build the LF of
    field_size: survey size (Ω)
    z1, z2: redshift slice (minimum and maximum redshift)
    Lbins, min_logLbin_value, max_logLbin_value: number of luminosity bins, minimum and maximum luminosity value in log scale
    tile: ID of the sub-field that will be used to construct the LF
    variance: If True it will compute the field-to_field variance at each luminosity bin, False is the default
    plot: if True it will plot the resulting LF and/or the field-to-field variance
    '''
    
    if variance:
        
        if field_size > 1:
            print('Computing the variance for field size > 1 is not offered!!! Please select a smaller field size.')
            sys.exit()
        else:
            phis = np.empty((117, Lbins))
            for tile in range(117):
                df, dV = load_catalogs(SIDES_cats_path, field_size, z1, z2, tile)
                Lline, units, transition = line_computation(df, line, J)
                Ls, phi = gen_lf(Lline, field_size, z1, z2, Lbins, min_logLbin_value, max_logLbin_value)
                phis[tile, :] = phi
            mean, median, std = np.mean(phis, axis = 0), np.median(phis, axis = 0), np.std(phis, axis = 0)
            p5, p16, p84, p95 = np.percentile(phis, [5, 16, 84, 95], axis = 0)
            rPearson = np.corrcoef(phis, rowvar = False)
        
            if plot:
                plt.figure(figsize = (6,3), dpi = 150)
                plt.loglog(Ls, mean, c = 'k', lw = 2, label = 'mean')
                plt.loglog(Ls, median, c = 'k', ls = '--', lw = 2, label = 'median')
                plt.fill_between(Ls, p5, p95, color='k', alpha = 0.3, label = '2 sigma')
                plt.fill_between(Ls, p16, p84, color='k', alpha = 0.6, label = '1 sigma')
                plt.annotate('{} @ {} < z < {}'.format(line + transition, z1, z2), xy=(0.5,0.9), xycoords='axes fraction')
                plt.xlabel(r'$\rm L_{} \, [{}]$'.format('{' + line + transition + '}', units))
                plt.ylabel(r'$\rm \Phi (L) \, [Mpc^{-3} dex^{-1}]$')
                plt.legend(loc = 3)
                plt.show()
        
    else:
        
        df, dV = load_catalogs(SIDES_cats_path, field_size, z1, z2, tile)
        Lline, units, transition = line_computation(df, line, J)
        Ls, phi = gen_lf(Lline, field_size, z1, z2, Lbins, min_logLbin_value, max_logLbin_value)
        mean, median, std, p5, p16, p84, p95, rPearson = [0 for _ in range(8)]
        
        if plot:
            plt.figure(figsize = (6,3), dpi = 150)
            plt.loglog(Ls, phi, c = 'k', lw = 2)
            plt.annotate('{} @ {} < z < {}'.format(line + transition, z1, z2), xy = (0.5,0.9), xycoords = 'axes fraction')
            plt.xlabel(r'$\rm L_{} \, [{}]$'.format('{' + line + transition + '}', units))
            plt.ylabel(r'$\rm \Phi (L) \, [Mpc^{-3} dex^{-1}]$')
            plt.show()
            
    return Ls, phi, mean, median, std, p5, p16, p84, p95, rPearson


def rho_mol_computation(SIDES_cats_path, J, z1, z2, field_size, aco = 3.6, variance = True):

    '''
    INPUT: CO rotational transition (J), redshift (z) and survey size (Ω)
    OUTPUT: ρ_mol(z,Ω) and std[ρ_mol(z,Ω)]
    If variance is False then it only computes the rho_mol value for the selected J, z, and Ω (fast computation)
    '''
    
    if variance:
    
        L_integral = np.empty(117)
        for tile in range(117):
            df, dV = load_catalogs(SIDES_cats_path, field_size, z1, z2, tile)
            L_integral[tile] = df.sum(df.LprimCO10) / dV

        mean_rho_mol, median_rho_mol, std_rho_mol = aco * np.mean(L_integral), aco * np.median(L_integral), aco * np.std(L_integral)
        p5_rho_mol, p16_rho_mol, p84_rho_mol, p95_rho_mol = aco * np.percentile(L_integral,[5,16,84,95])
        rho_mol = 0
        
    else:
        mean_rho_mol, median_rho_mol, std_rho_mol, p5_rho_mol, p16_rho_mol, p84_rho_mol, p95_rho_mol = [0 for _ in range(7)]
        df = load_catalogs(SIDES_cats_path, field_size, z1, z2, tile)
        rho_mol = aco * df.sum(df.LprimCO10) / dV

    return rho_mol, mean_rho_mol, median_rho_mol, std_rho_mol, p5_rho_mol, p16_rho_mol, p84_rho_mol, p95_rho_mol
