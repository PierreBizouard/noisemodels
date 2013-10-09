#noisemodels.py

rev 0 oct 09 2013

... model some laser noise

### GaussianShotnoise( field,lambda, tvec )
* models shot noise based on poisson distribution for a given **field** at the wavlength **lambda**

### GetOppmField( omvec, tvec)
* one-photon-per-mode

### getrandomphase3( Bandwidth, FreqJitterFactor, nufwhm, tvec)
* classical phase diffusion phase

### PdModellGauss3 field, FreqJitterFactor, GaussSpecWidthFactor, nufwhm , relomvec,  tvec  )
* returns a Phase diffusion noise with gaussian spectrum
* see [1]M. H. Frosz, “Validation of input-noise model for simulations of supercontinuum generation and rogue waves,” Opt. Express, vol. 18, no. 14, pp. 14778–14787, Jul. 2010.
