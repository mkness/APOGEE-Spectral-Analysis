# this file is part of the XXX project
# Copyright 2014 Melissa Ness

import scipy
from scipy import interpolate 
from scipy import ndimage 
# import numpy as np

def get_data():

  T_est,g_est,feh_est = loadtxt("starsin_new.txt", usecols = (5,7,9), unpack =1) 
  thismeta = array([T_est, feh_est, g_est])
  thismeta = [T_est, feh_est, g_est]
  #thismeta = zip(T_est, feh_est, g_est)
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
    print atleast_2d(a[1].data).shape
    if jj == 0:
      nmeta = len(thismeta)
      nlam = len(a[1].data[0])
    #val = diff_wl*(len_data) + start_wl 
    val = diff_wl*(nlam) + start_wl 
    wl_full_log = arange(start_wl,val, diff_wl) 
    ydata = (atleast_2d(a[1].data))[0] 
    ydata_err = (atleast_2d(a[2].data))[0] 
    ydata_flag = (atleast_2d(a[3].data))[0] 
    assert len(ydata) == nlam
    wl_full = [10**aval for aval in wl_full_log]
    testdata = scipy.ndimage.gaussian_filter(ydata, 20 ) 
    xdata= wl_full
    xdata =array(xdata)
    ydata =array(ydata)
    ydata_err =array(ydata_err)
    ydata_flag =array(ydata_err)
    #takevalues = logical_and(ydata_err > 10, logical_and(ydata_flag !=1, ydata_flag <= 4500)  
    takevalues = logical_and(ydata  >20, logical_and(ydata_flag >1, ydata_flag <= 3000) )  
    takevalues = logical_and(ydata > 20, logical_and(ydata_flag >1, ydata_flag <= 3000) ) 
    a1 = xdata#[takevalues]
    b1 = ydata#[takevalues]
    b2 = scipy.ndimage.gaussian_filter(b1,1) 
    b3 = scipy.ndimage.gaussian_filter(b1[takevalues],100) 
    a1 = array(a1) 
    b1 = array(b1) 
    fit1 = pylab.polyfit(array(a1[takevalues]), array(b3), 3) 
    y1 = polyval(array(fit1), a1) 
    ynew = b2/y1
    y2new = b2/y1
    xgrid1 = a1 
    starname2 = each.split('.fits')[0]+'.txt'
    ivarnew = ones_like(y2new)
    if jj == 0:
      npix = len(xgrid1) 
      dataall = zeros((npix, len(bl), 3))
      #metaall = zeros((npix, len(bl), 3))
      #metaall = zeros( (npix, len(bl), 3) )
      #metaall = zeros((len(b1), nmeta))
      #metaall = ones((len(bl), nmeta+1))
      metaall = ones((len(bl), nmeta))
    if jj > 0:
      assert xgrid1[0] == dataall[0, 0, 0]
    dataall[:, jj, 0] = xgrid1
    dataall[:, jj, 1] = y2new
    dataall[:, jj, 2] = ivarnew
    metaall[:, 0] = T_est
    metaall[:, 1] = g_est
    metaall[:, 2] = feh_est
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
  M = hstack((ones([nobj, 1]), meta))
  MTCinvM = dot(M.T, Cinv[:, None] * M) # craziness b/c Cinv isnt a matrix
  x = data[:, 1] # intensity slice of data
  MTCinvx = dot(M.T, Cinv * x)
  print 1 
  return linalg.solve(MTCinvM, MTCinvx)

def do_regressions(dataall, metaall):
  return map(do_one_regression, dataall, metaall)

dataall, metaall = get_data()
coeffs = do_regressions(dataall, metaall)

