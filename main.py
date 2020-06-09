from read_and_make import *
from ILC_helpers import *
import configparser
import json
import os

def goWMAP():

    band = eval(WMAP.get('band'))
    band_beam = eval(WMAP.get('band_beam'))
    lmax = int(WMAP.get('lmax'))
    nside = int(WMAP.get('nside'))
    imap_dir = WMAP.get('imap_directory')
    beam_dir = WMAP.get('beam_directory')

    print('>Read imap')
    RW = Read_WMAP(band, lmax, nside)
    imaps = RW.read_imap(imap_dir)
    print('>Change imap to alm')
    alms = RW.map_to_alm(imaps)
    print('>Read beam')
    beam_transfer_functions = RW.read_beam(band_beam, beam_dir)
    print('>Calculate　Cl')
    ILCW = ILC_WMAP(alms, lmax, beam_transfer_functions)
    Cl = ILCW.Cl_cal()
    print('>Calculate weight')
    weights = ILCW.weight_cal(Cl)
    print('>Do ILC')
    alm_clean = ILCW.do_ILC(weights)
    print('>Show clean cmb map')
    Show(alm_clean, nside, lmax, fwhm = 1.0).moview()

def goPlank():

    band = json.loads(Plank.get('band'))
    fwhm = np.array(json.loads(Plank.get('fwhm')))
    model = eval(Plank.get('model'))
    nside = int(Plank.get('nside'))
    lmax = 3 * nside - 1

    print('>Create imap')
    MP = Make_Plank(band, lmax, nside, fwhm, model)
    imaps = MP.make_imap()
    print('>Change imap to alm')
    alms = MP.map_to_alm(imaps)
    print('>Calculate　Cl')
    ILCP = ILC_Plank(alms, lmax)
    Cl = ILCP.Cl_cal()
    print('>Calculate weight')
    weights = ILCP.weight_cal(Cl)
    print('>Do ILC')
    alm_clean = ILCP.do_ILC(weights)
    print('>Show clean cmb map')
    Show(alm_clean, nside, lmax, fwhm = 1.0).moview()


if __name__ == '__main__':

    # setup(read config.ini)
    config = configparser.ConfigParser()
    config.read('config.ini')
    WMAP = config['WMAP']
    Plank = config['Plank']

    if WMAP.get('WMAP') == 'True':
        print('WMAP')
        print('--------------------------------')
        print('Do ILC, using WMAP data')
        goWMAP()

    elif Plank.get('Plank') == 'True':
        print('Plank')
        print('--------------------------------')
        print('Do ILC, using Plank data')
        goPlank()

    else:
        print('ERROR: you have to choose WMAP or Plank or LiteBIRD = True in config.ini')
