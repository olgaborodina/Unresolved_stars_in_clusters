# N-количесетво строк в функции светимости, K-количество строк в изохроне 
# vis_mag_char, fi_char,mass_char,abs_mag_char - вид зв величина звезд скопления, количество звезд (значение функции светимости), масса, абсолютная зв величина - для считывания из файла
# константы: EBV - погл. r - расстояние до скопления ALPHA - доля двойных в скоплении, step-шаг с которым строится гистограмма для подсчета массы
# vis_mag: 0- видимая зв величина, начиная с которой рассматривается гиста, fin-предельная зв величина, 1,2 - вид зв величины для интерполяции, чтобы найти значение функции светимости на середине отрезка vis_mag0 и vis_mag0+step
# fi: 1,2 значение функции в точках vis_mag1 и vis_mag2, fi- интерполированное значение функции в середине отрезка vis_mag0 и vis_mag0+step
# abs_mag: 0 - абсолютная зв величина из vis_mag0, 1,2 - абс зв величины для интерполяции, чтобы найти значение массы для abs_mag0
# mass: 1,2 - масса в точках abs_mag1 и abs_mag2, mass- интерполированное значение массы в точке abs_mag0
# MASS_DOUBLE: _i- масса ТОЛЬКО ДВОЙНЫХ звезд из интервала vis_mag0 и vis_mag0+step, MASS_DOUBLE - масса ВСЕХ ДВОЙНЫХ звезд от начальной до конечной звездной величины
# MASS_SINGLE: _i- масса ТОЛЬКО ОДИНОЧНЫХ звезд из интервала vis_mag0 и vis_mag0+step, MASS_SINGLE - масса ВСЕХ ОДИНОЧНЫХ звезд от начальной до конечной звездной величины
# MASS_CLUSTER_NEW: масса скопления с приближением, что есть двойные
# MASS_CLUSTER_OLD: масса скопления с приближением, что все одиночные
# q-отношение массы главного компонента, к вторичному

import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
import math as m
from scipy.interpolate import interp1d
from random import random
import scipy.optimize
import csv
 
def lummass (M):
    return pow(10, (-0.705)*np.log10(M)*np.log10(M) + 4.655*np.log10(M) - 0.025)

def get_parameters(mode): 
    if (mode == 'IC2714'):
        return 10.484, 0.340
    elif (mode == 'NGC1912'):
        return 10.294, 0.253
    elif (mode == 'NGC2099'):
        return 10.74, 0.301
    elif (mode == 'NGC6834'):
        return 11.588, 0.706
    elif (mode == 'NGC7142'):
        return 10.2, 0.25
    else:
        raise ValueError 
      
def to_lf (a, b, c, step, alpha):
    func = interp1d(a, b)
    return func(c) * step * alpha

def to_mf (a, b, c, d):
    func = interp1d(a, b)
    return func(c)*d

def neyman(data_q):
    iter = True
    while iter == True:
        q = random()
        Fr = random()
        F = interp1d(data_q['q'], data_q['F'])
        if Fr <= F(q):
            iter = False
            return q
        else:
            iter = True
def func(x,q,lf):
    return np.log10(lummass(x))*np.log(10) + np.log(1 + 
                   pow(m.e, np.log(10)*(-0.705*np.log10(q)*np.log10(q) +
                4.655*np.log10(q) - 1.41*np.log10(x)*np.log10(q))))-np.log(lf)    
        
in_file_lf = input("Input the Luminosity Function")
in_file_is = input("Input the isochrone")
in_file_q = "func_gauss.txt"

dist_mod, EBV = get_parameters (in_file_lf.split('_')[0])

data_lf = pd.read_csv(in_file_lf, delimiter=' ', header=None, names = ['J, mag','F','F_low', 'F_high'])
data_is = pd.read_csv(in_file_is, delimiter=' ', header=None, names = ['M, SM','J, mag'])
data_q = pd.read_csv(in_file_q, delimiter=' ', header=None, names = ['q','F'])

data_lf['Jabs, mag'] = data_lf['J, mag']- 2.43 * 0.37 * EBV - dist_mod

step = (data_lf['Jabs, mag'].iloc[-1] - data_lf['Jabs, mag'].iloc[0]) / 30
main_table = pd.DataFrame(np.linspace(data_lf['Jabs, mag'].iloc[0], data_lf['Jabs, mag'].iloc[-1], num = 31, dtype=None)+step/2)
main_table.rename(columns={0 : 'Jabs'}, inplace = True)
main_table = main_table.drop(main_table.index[30])
f = open('output.txt', 'w+')
results = pd.DataFrame(data = {'ALPHA': [0], 'MASS CLUSTER OLD': [0], 'MASS CLUSTER NEW': [0], 'RATIO': [0] })
for ALPHA in np.arange (0., 1., 0.1):
    ratio = 0
    MASS_BINARIES = 0
    
    #getting the mass of the cluster without any binary
    main_table['Num_all_ss'] = np.round(to_lf(data_lf['Jabs, mag'], data_lf['F'], main_table['Jabs'], step, 1))
    main_table['mass_all_ss'] = to_mf(data_is['J, mag'], data_is['M, SM'], main_table['Jabs'], main_table['Num_all_ss'])
    MASS_CLUSTER_OLD = main_table['mass_all_ss'].sum()
    
    #getting the mass of binaries in the cluster
    main_table['Num_bin'] = np.round(to_lf(data_lf['Jabs, mag'], data_lf['F'], main_table['Jabs'], step, ALPHA))
    main_table['Lum0_bin'] = lummass(to_mf(data_is['J, mag'], data_is['M, SM'], main_table['Jabs'],1))
    
    data_binaries = main_table.loc[main_table.index.repeat(main_table['Num_bin'].astype(int))]
    data_binaries.index = pd.RangeIndex(len(data_binaries.index))
    
    for star in range (0, len(data_binaries)):
        q=neyman(data_q)
        MASS_BINARIES += scipy.optimize.brentq(func, 0.001, 10, args = (q,data_binaries['Lum0_bin'].iloc[star])) * (1 + q)
    
    main_table['Num_part_ss'] = np.round(to_lf(data_lf['Jabs, mag'], data_lf['F'], main_table['Jabs'], step, 1)) - np.round(to_lf(data_lf['Jabs, mag'], data_lf['F'], main_table['Jabs'], step, ALPHA))
    main_table['mass_part_ss'] = to_mf(data_is['J, mag'], data_is['M, SM'], main_table['Jabs'], main_table['Num_part_ss'])
    MASS_SINGLES = main_table['mass_part_ss'].sum()
   
    MASS_CLUSTER_NEW = MASS_BINARIES + MASS_SINGLES;

    results = results.append(pd.Series([ALPHA, MASS_CLUSTER_OLD, MASS_CLUSTER_NEW, MASS_CLUSTER_NEW/MASS_CLUSTER_OLD], index =['ALPHA', 'MASS CLUSTER OLD', 'MASS CLUSTER NEW', 'RATIO']), ignore_index=True)
   #ratio += MASS_CLUSTER_NEW/MASS_CLUSTER_OLD
results = results.drop(0)
results.round({'ALPHA': 1, 'MASS CLUSTER OLD': 2, 'MASS CLUSTER NEW': 2, 'RATIO': 4 })
results.to_csv(f"alpha_{in_file_q.split('.txt')[0]}.txt"  , sep='\t', index=False)
    # finding error despertion
    