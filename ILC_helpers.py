import healpy as hp
import numpy as np
import matplotlib.pyplot as plt

class ILC:

    def __init__(self, alms, lmax):

        self.alms = alms
        self.lmax = lmax
        self.len = len(self.alms)

    def Cl_cal(self):

        Cl = np.zeros((self.lmax+1, self.len, self.len))
        for i in range(self.len):
            for j in range(self.len):
                pre_Cl = hp.alm2cl(self.alms[i], self.alms[j], lmax = self.lmax)

                for l in range(self.lmax+1):
                    Cl[l][i][j] = pre_Cl[l]

        return Cl

    def weight_cal(self, Cl):

        e = np.ones((self.len, 1))
        weights = []
        for l in range(self.lmax+1):

            A = np.linalg.solve(Cl[l], e)
            weight = A/np.dot(e.T, A)
            weights.append(weight)

        weights_upd = weights[0]

        for l in range(1, self.lmax+1):
            weights_upd = np.concatenate([weights_upd, weights[l]], axis = 1)

        return weights_upd


    def do_ILC(self, weights):

        alm_clean = np.zeros(len(self.alms[0]), dtype = complex)
        for l in range(self.lmax+1):
            for m in range(l+1):
                j = hp.sphtfunc.Alm.getidx(self.lmax, l, m)

                for i in range(self.len):
                    alm_clean[j] += weights[i][l] * self.alms[i][j]

        return alm_clean



class ILC_WMAP(ILC):

    def __init__(self, alms, lmax, beam_transfer_functions):
        super().__init__(alms, lmax)
        self.btfs = beam_transfer_functions

class ILC_Plank(ILC):
    pass



class Show:

    def __init__(self, alm_clean, nside, lmax, fwhm = 0.0):

        self.alm_clean = alm_clean
        self.nside = nside
        self.lmax = lmax
        self.fwhm = fwhm

    def moview(self):

        map_clean = hp.alm2map(self.alm_clean, nside = self.nside, fwhm = self.fwhm*np.pi/180, verbose=False)
        hp.mollview(map_clean, title = 'CMB map', cmap = 'jet', max = 0.2, min = -0.2)

        plt.show()

    def Dl(self):

        l = np.arange(self.lmax+1)
        Cl_clean = hp.alm2cl(self.alm_clean, lmax = self.lmax)
        Dl = Cl_clean*l*(l+1)/(2*np.pi)*10**6

        plt.plot(l, Dl, label = 'Dl')
        plt.xscale('log')

        plt.xlabel('Multipole l')
        plt.ylabel('Dl[μK^2]')

        plt.xlim([2, self.lmax])
        plt.ylim([0, 6000])
        plt.show()
