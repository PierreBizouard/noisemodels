


import numpy as np


def GaussianShotnoise( field,lam, tvec ):
    c = 2.99792458e8
    h = 6.62606959e-34
    dt = tvec[1]-tvec[0]
    Ephot = h * c / lam
    fieldampl = np.sqrt(Ephot)
    #number of photons per time bin
    Np = np.abs(field)**2 *dt / Ephot
    noisefield = np.zeros(np.shape(field),dtype=complex)
    newfield = np.zeros(np.shape(field),dtype=complex)
    lenfield = len(field)
    for i in range(lenfield):
#        if i%5000==0:
#            print "GaussianNoise %.2f percent done."%(i*100.0/lenfield)
        std = np.round( np.sqrt(Np[i]) )
        if std > 0.0:
            number_noise_photons =  int( np.abs( np.random.normal( 0, std, 1) ))
            phasen = np.random.rand( number_noise_photons ) * 2 * np.pi
            noisefield[i] = fieldampl * np.sum( np.exp( 1.0j * phasen ))
         #   print "Photonen %d    Rauschphotonen  %d    |noisefield|**2 %.f "%(Np[i],number_noise_photons, np.abs(noisefield[i])**2)
    newfield = field + noisefield
    return newfield




def GetOppmField( omvec, tvec):
    """def GetOppmField( omvec, tvec):
    returns the temporal field for one-photon-per-mode- noise"""
    tvec= np.ndarray.flatten( tvec)
    omvec= np.ndarray.flatten( omvec)
    c = 2.99792458e8
    h = 6.62606957e-34
    phasen=np.random.random(len(tvec))*2*np.pi
    dom = omvec[1]-omvec[0]
    dt = tvec[1]-tvec[0]
    photonenenergien = h * omvec / np.pi / 2
    spekfeld = np.multiply( np.sqrt(photonenenergien/dom), np.exp(1.0j * phasen ))
    oppfeld = np.fft.ifft( (np.fft.fftshift(spekfeld))[::-1]) * np.sqrt( len(omvec)  * dom / dt)
#    oppfeld = np.fft.ifft( (np.fft.fftshift(spekfeld))) * np.sqrt( len(omvec)  * dom / dt) <- this yields a spectrally mirrored field. not good
 #   print "spektrale Energie = ",np.sum( np.abs(spekfeld)**2 *dom)
 #   print "zeitl. Energie    = ",np.sum( np.abs(oppfeld)**2 * dt)
    return oppfeld



def getrandomphase3( Bandwidth, FreqJitterFactor, nufwhm, tvec): #Sept,2,2013
    dt = tvec[1] - tvec[0]
    sigmasq = FreqJitterFactor * nufwhm / 2 / np.pi * Bandwidth
    sigma = np.sqrt(sigmasq) #random.normal expects the standard deviation, not the variance!
    nur = np.random.normal( 0 ,sigma, np.shape(tvec) )
    phase = np.cumsum( nur ) * dt * 2 * np.pi
    return phase

def PdModellGauss3( field, FreqJitterFactor, GaussSpecWidthFactor, nufwhm , relomvec,  tvec  ):
    """
def PdModellGauss3( feld, FreqJitterFactor, GaussSpecWidthFactor, nufwhm , relomvec,  tvec  ):
     """

    dom = relomvec[2]-relomvec[1]
    N = len(tvec)
    SimBandwidth = (np.max(relomvec) - np.min( relomvec))/2/np.pi
    #
    # phase diffusion modell
    #
    ppPD = getrandomphase3(SimBandwidth, FreqJitterFactor, nufwhm, tvec)
    fPD = np.multiply(field,np.exp( 1.0j * ppPD))
    ftPD0 = np.fft.fftshift(np.fft.fft(fPD))
    oldenergy = np.sum(np.abs(ftPD0)**2) * dom
    #
    # shape spectrum to gauss
    #
    faktor1 = 2.4295 / (nufwhm * GaussSpecWidthFactor)
    faktor2 = np.exp(  -0.5 * (relomvec/np.pi/2)**2 / (GaussSpecWidthFactor * nufwhm/1.6651)**2)
    faktor3 = np.sqrt( (relomvec/np.pi/2)**2 + (GaussSpecWidthFactor * nufwhm/2)**2)
    ftPD1 = np.multiply( ftPD0, faktor1*faktor2*faktor3)
    #
    # rescale to preserve energy
    #
    newenergy = np.sum(np.abs(ftPD1)**2) * dom
##    print "oldenergy ",oldenergy
###    print "newenergy ",newenergy
    ftPD1 = ftPD1 * np.sqrt( oldenergy / newenergy)   
###    print "newenergy2",np.sum(np.abs(ftPD1)**2)*dom
    fPD1 = np.fft.ifft(np.fft.fftshift(ftPD1))
    return fPD1
