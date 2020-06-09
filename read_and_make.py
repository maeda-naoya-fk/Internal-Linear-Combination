import healpy as hp
import numpy as np
import pysm
import pysm.units as u
import os
import warnings
warnings.simplefilter('ignore')


class Read_WMAP:

    def __init__(self, band, lmax, nside):

        self.band = band
        self.lmax = lmax
        self.nside = nside

    def read_imap(self, imap_dir):

        path = os.path.dirname(__file__) + '/' + imap_dir
        os.chdir(path)

        imaps = []
        for b in self.band:

            imap_fit = 'wmap_band_imap_r9_9yr_' + b + '_v5.fits'
            imap = hp.read_map(imap_fit, verbose = False)
            imaps.append(imap)

        return np.array(imaps)

    def map_to_alm(self, imaps):

        alms = []
        for imap in imaps:

            alm = hp.map2alm(imap, lmax = self.lmax, iter = 10)
            alms.append(alm)

        return np.array(alms)

    def read_beam(self, band_beam, beam_dir):

        path = os.path.dirname(__file__) + '/' + beam_dir
        os.chdir(path)

        beam_transfer_functions = []
        for bb in band_beam:

            beam_txt = 'wmap_ampl_bl_' + bb + '_9yr_v5p1.txt'
            with open(beam_txt, 'r') as f:
                lines = f.readlines()

            beam_transfer_function = np.zeros(self.lmax+1)
            i = 0
            l = 0
            while l < self.lmax+1:

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


class Make_Plank(Read_WMAP):

    def __init__(self, band, lmax, nside, fwhm, model):
        super().__init__(band, lmax, nside)
        self.fwhm = fwhm
        self.model = model

    def make_imap(self):

        sky = pysm.Sky(nside = self.nside, preset_strings = self.model)

        imaps = []
        for b, f in zip(self.band, self.fwhm):

            map = sky.get_emission(b * u.GHz)
            imap = hp.smoothing(map[0], fwhm = f/60*np.pi/180, verbose = False)
            imaps.append(imap)

        return np.array(imaps)
