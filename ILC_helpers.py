import healpy as hp
import numpy as np
import os
import matplotlib.pyplot as plt

def read_map(band):

    path = os.path.dirname(__file__)+'/wmap_band_imap_r9_nineyear_v5'
    os.chdir(path)

    maps = []
    for b in band:

        map_fits =  'wmap_band_imap_r9_9yr_' + b + '_v5.fits'
        map = hp.read_map(map_fits)
        maps.append(map)

    return np.array(maps)

def map2alm(maps, lmax):

    alms = []
    for map in maps:

        alm = hp.map2alm(map, lmax = lmax, iter = 10)
        alms.append(alm)

    return np.array(alms)

def read_band_beam(band_beam, lmax):

    path = os.path.dirname(__file__)+'/wmap_ampl_bl_9yr_v5p1'
    os.chdir(path)

    beam_transfer_functions = []
    for bb in band_beam:

        beam_txt = 'wmap_ampl_bl_' + bb + '_9yr_v5p1.txt'
        with open(beam_txt, 'r') as f:
            lines = f.readlines()

        beam_transfer_function = np.zeros(lmax+1)
        i = 0
        l = 0
        while l < lmax+1:

            if lines[i][0] != '#':
                beam_transfer_function[l] = float(lines[i].split()[1])
                l += 1

            i += 1

        beam_transfer_functions.append(beam_transfer_function)

    beam_transfer_functions = np.array(beam_transfer_functions)

    K_ave = beam_transfer_functions[0]
    Ka_ave = beam_transfer_functions[1]
    Q_ave = sum(beam_transfer_functions[2:4])/2
    V_ave = sum(beam_transfer_functions[4:6])/2
    W_ave = sum(beam_transfer_functions[6:10])/4

    return np.array([K_ave, Ka_ave, Q_ave, V_ave, W_ave])

class ILC:

    def __init__(self, alms, beam_transfer_functions, lmax):

        self.alms = alms
        self.btfs = beam_transfer_functions
        self.lmax = lmax

    def Cl_cal(self):

        Cl = np.zeros((self.lmax+1, 5, 5))
        for i in range(5):
            for j in range(5):
                pre_Cl = hp.alm2cl(self.alms[i], self.alms[j], lmax = self.lmax)

                for l in range(self.lmax+1):
                    Cl[l][i][j] = pre_Cl[l]

        return Cl

    def weight_cal(self, Cl):

        e = np.ones((5, 1))
        weights = []
        for l in range(self.lmax+1):

            A = np.linalg.solve(Cl[l], e)
            weight = A/np.dot(e.T, A)
            weights.append(weight)

        #weightの配列を直す
        weights_upd = weights[0]

        for l in range(1, self.lmax+1):
            weights_upd = np.concatenate([weights_upd, weights[l]], axis = 1)

        return weights_upd


    def Do_ILC(self, weights):


        alm_clean = np.zeros(len(self.alms[0]), dtype = complex)
        for l in range(self.lmax+1):
            for m in range(l+1):
                j = hp.sphtfunc.Alm.getidx(self.lmax, l, m)

                for i in range(5):
                    alm_clean[j] += weights[i][l] * self.alms[i][j]

        return alm_clean

def alm2map(alm_clean, Nside, fwhm = 0.0):

    map_clean = hp.alm2map(alm_clean, nside = Nside, fwhm = fwhm*np.pi/180)

    return map_clean

def Dl(alm_clean, lmax):

    l = np.arange(lmax+1)
    Cl_clean = hp.alm2cl(alm_clean, lmax = lmax)
    Dl = Cl_clean*l*(l+1)/(2*np.pi)*10**6

    plt.plot(l, Dl, label = 'Dl')
    plt.xscale('log')

    plt.xlabel('Multipole l')
    plt.ylabel('Dl[μK^2]')

    plt.xlim([2,751])
    plt.ylim([0, 6000])
    plt.show()
