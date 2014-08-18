# this file is part of the XXX project
# Copyright 2014 Melissa Ness

import scipy
from scipy import interpolate 
from scipy import ndimage 
# import numpy as np

def get_data():

  T_est,g_est,feh_est = loadtxt("starsin_new.txt", usecols = (5,7,9), unpack =1) 
  thismeta = array([T_est, feh_est, g_est])
  a = open('starsin_new.txt','r')
  al = a.readlines() 
  bl = []
  for each in al:
    bl.append(each.split()[0]) 
  for jj,each in enumerate(bl):
    each = each.strip('\n')
    a = pyfits.open(each) 
    b = pyfits.getheader(each) 
    start_wl =  a[0].header['CRVAL1']
    diff_wl = a[0].header['CDELT1']
    #T_est = a[0].header['RVTEFF']
    #Feh_est =a[0].header['RVFEH']
    #g_est  = a[0].header['RVLOGG']
    #thisdata = atleast_2d(a[1].data)[0] 
    print atleast_2d(a[1].data).shape
    if jj == 0:
      nmeta = len(thismeta)
      nlam = len(a[1].data[0])
    val = diff_wl*(len_data) + start_wl 
    wl_full_log = arange(start_wl,val, diff_wl) 
    ydata = (atleast_2d(a[1].data))[0] 
    assert len(ydata) == nlam
    wl_full = [10**aval for aval in wl_full_log]
    testdata = scipy.ndimage.gaussian_filter(ydata, 20 ) 
    xdata= wl_full
    xdata =array(xdata)
    ydata =array(ydata)
    a1 = xdata
    b1 = ydata
    b2a = b1
    b2 = scipy.ndimage.gaussian_filter(b2a,1) 
    b3 = scipy.ndimage.gaussian_filter(b2a,20) 
    b4 = scipy.ndimage.gaussian_filter(b2a,100) 
    a1 = array(a1) 
    b4 = array(b4) 
    xdata2 = array(list(a1[600:800])+list(a1[900:2700])+ list(a1[4000:5800])+ list(a1[6800:8150]))
    ydata2 = array(list(b3[600:800])+list(b3[900:2700])+ list(b3[4000:5800])+ list(b3[6800:8150]))
    #fit1 = pylab.polyfit(array(xdata2), array(ydata2), 25) 
    fit1 = pylab.polyfit(array(xdata2), array(ydata2), 5) 
    xgrid1 = a1[365:-450]
    #y1 = polyval(array(fit1), xgrid1) 
    y1 = polyval(array(fit1), xgrid1) 
    y2new = b2[365:-450]#/y1
    starname2 = each.split('.fits')[0]+'.txt'
    # this is a hack and must be replaced 
    ivarnew = ones_like(y2new)
    if jj == 0:
      npix = len(xgrid1) 
      dataall = zeros((npix, len(bl), 3))
      metaall = zeros((len(b1), nmeta))
    if jj > 0:
      assert xgrid1[0] == dataall[0, 0, 0]
    dataall[:, jj, 0] = xgrid1
    dataall[:, jj, 1] = y2new
    dataall[:, jj, 2] = ivarnew
    metaall[jj, :] = thismeta
  return dataall, metaall

def do_one_regression(data, meta):
  """
  ## inputs
  - data [nobjs, 3] wavelengths, fluxes, invvars
  - meta [nobjs, nmeta] Teff, Feh, etc, etc

  ## outputs
  - coefficients of the fit
  """
  nobj, nmeta = meta.shape
  assert data.shape == (nobj, 3)
  # least square fit
  Cinv = data[:, 2] # invvar slice of data
  M = hstack((ones(nobj, 1), meta))
  MTCinvM = dot(M.T, Cinv[:, None] * M) # craziness b/c Cinv isnt a matrix
  x = data[:, 1] # intensity slice of data
  MTCinvx = dot(M.T, Cinv * x)
  print 1 
  return linalg.solve(MTCinvM, MTCinvx)

def do_regressions(dataall, metaall):
  return map(do_one_regression, dataall, metaall)

dataall, metaall = get_data()
coeffs = do_regressions(dataall, metaall)

