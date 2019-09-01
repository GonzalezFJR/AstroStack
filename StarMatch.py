'''
  @ autor: Xuan
  @ mail: gonzalezrodrigo@uniovi.es

'''

from math import sqrt, atan2, sin, cos, pi
import numpy as np
import os, sys
from photutils import datasets
from photutils import DAOStarFinder
import matplotlib.pyplot as plt
import astropy
from matplotlib.pyplot import figure
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from photutils import aperture_photometry, CircularAperture
from astropy.visualization import *

#=========== Read fits and extract info
#========================================================================

def LoadFit(path):
  if not os.path.isfile(path):
    print 'ERROR: not found file ' + path
    return[0]
  hdu = astropy.io.fits.open(path)
  return hdu


def LoadFitImage(path, area = []):
  ''' Opens a .fit image and returns the float 2D array and the filter (read from the header) '''
  if not os.path.isfile(path): 
    print 'ERROR: not found file ' + path
    return[0]
  hdu = astropy.io.fits.open(path)[0]  #datasets.load_star_image()
  #hdu.verify('fix')
  filter = ''
  if 'FILTER' in hdu.header: filter = hdu.header['FILTER']
  image = hdu.data.astype(float)
  if not len(area) > 0:
    return image, filter
  x0, x1, y0, y1 = area
  recImg = hdu.data[x0:x1, y0:y1].astype(float)
  return image, filter, recImg

def GetBkgMedianAndSigma(image):
  ''' Obtains the background as median of the image and the sigma '''
  from astropy.stats import mad_std
  bkg = np.median(image)
  bkg_sigma = mad_std(image)
  return [bkg, bkg_sigma]

def GetSources(image, fwhm = 5., thr = 3, dataformat = '%.2g'):
  ''' Find stars with a given fwhm and threshold '''
  bkg,bkgsigma = GetBkgMedianAndSigma(image)
  daofind = DAOStarFinder(fwhm = fwhm, threshold = thr*bkgsigma)
  sources = daofind(image)
  for col in sources.colnames: sources[col].info.format = dataformat
  return sources

def GetBrightestSource(a):
  ''' Find the brightest star in an image or a group of stars '''
  if   isinstance(a, astropy.table.table.Table): # Table of stars
    maxflux = max(a['flux'])
    minmag  = min(a['mag'])
    bs = a[0]
    for s in a:
      #if s['flux'] == maxflux: bs = s
      if s['mag']  == minmag: bs = s
    return bs
  elif isinstance(a, np.ndarray): # Image
    return GetBrightestSources(GetSources(a))
  else:
    print '[GetBrightestSources] ERROR: unknown type ', type(a)

def GetListOfBightestSources(sources, nSources = 10, area = []):
  x = sources['xcentroid']; y = sources['ycentroid']; coor = np.transpose([x,y])
  magb = sources['mag']
  c = [z for _, z in sorted(zip(magb,coor))]
  xb,yb = np.transpose(c[:])#[:nSources])
  magb = sorted(magb)[:]#[:nSources]
  x = []; y = []; mag = []
  if len(area) >= 4:
    x0, x1, y0, y1 = area
    for i in range(len(xb)):
      xs = xb[i]; ys = yb[i];
      if xs < x1 and xs > x0 and ys < y1 and ys > y0:
        x.append(xs); y.append(ys); mag.append(magb[i])
    x = x[:nSources]; y = y[:nSources]; mag = mag[:nSources]
  else:
    x = xb[:nSources]; y = yb[:nSources]; mag = magb[:nSources]
  return [x, y, mag]

def GetSourceInPos(sources, px, py, thr = 5):
  ''' Looks for a source in [px, py] from all stars in the given list '''
  for s in sources:
    x = s['xcentroid']
    y = s['ycentroid']
    if abs(x-px) < thr and abs(y-py) < thr: return s
  print '[GetSourceInPos] WARNING: not found source in position [x, y] = [%i, %i]'%(px, py)
  return 0

def GetStarInPos(stars, px, py, thr = 5):
  ''' Looks for a star in [px, py] from all stars in the given list '''
  mindist = 1000000
  for s in stars:
    dist = s.GetDistToPoint(px, py)
    if dist < mindist:
      mindist = dist
      closestStar = s
    if mindist < thr: return s
    else: print '[GetStarInPos] WARNING: not found star in position [x, y] = [%i, %i]'%(px, py)
  return 0


################################################################################################
## Match stars

def GetBrightestStars(fname, n = 2, fwhm = 5, thr = 3):
  if not fname.endswith('.fit'): fname += '.fit'
  a = LoadFitImage(fname)
  img = a[0]
  sources = GetSources(img, fwhm=fwhm, thr=thr)
  return GetListOfBightestSources(sources, n)

def GetDistance(x1, y1 = '', x2 = '', y2 = ''):
  if y1 == '':
    x1, y1, x2, y2 = x1
  d = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)
  return sqrt(d)

def GetMeanStarsDistance(stars):
  v = []
  for s in stars: v.append(GetDistance(s))
  return np.average(v)

def GetClosest(x, y, slist):
  X, Y = slist
  n = len(X)
  mind = 100000; index = -1
  for i in range(n):
    d = GetDistance(x, y, X[i], Y[i])
    if d < mind:
      mind = d
      index = i
  return X[index], Y[index]

def GetMatchedList(mainList, secondList):
  mX = []; mY = []
  for x, y in zip(mainList[0], mainList[1]):
    xx, yy = GetClosest(x, y, secondList)
    mX.append(xx)
    mY.append(yy)
  return mX, mY

def GetGoodMatchedStars(file1, file2, dmin = 0.1, dmax = 500, n = 100, fwhm = 5, verbose = False, trying = 0):
  X1, Y1, m1 = GetBrightestStars(file1, n, fwhm=fwhm)
  X2, Y2, m2 = GetBrightestStars(file2, n, fwhm=fwhm)
  X2, Y2 = GetMatchedList([X1, Y1], [X2, Y2])
  outcoor = []
  outdist  = []
  for i in range(n):
    d = GetDistance(X1[i], Y1[i], X2[i], Y2[i])
    if d < dmin or d > dmax: continue
    coor = [X1[i], Y1[i], X2[i], Y2[i]]
    outcoor.append(coor); outdist.append(d)
    if verbose: print '[%1.1f, %1.1f]  [%1.1f, %1.1f], d = %1.2f'%(X1[i], Y1[i], X2[i], Y2[i], d)
  stars = [x for _,x in sorted(zip(outdist,outcoor))]
  dist  = [GetDistance(s) for s in stars]
  median = dist[len(dist)/2]
  mean   = sum(dist)/len(dist)
  cut = 1.2*median
  goodstars = []
  for s in stars: 
    if GetDistance(s) < cut: goodstars.append(s)
  if   len(goodstars) < (n/5):
    if trying < 3:
      print 'ERROR: matched %i good stars out of %i!! Repeating with %i stars...'%(len(goodstars), n, n*1.5)
      return GetGoodMatchedStars(file1, file2, dmin, dmax, int(n*1.5), fwhm=fwhm+trying, verbose=verbose, trying=trying+1)
    else:
      print 'ERROR: IMPOSIBLE TO MATCH STARS!!!!! We are very sorry... here you have out best attempt:'
      for i in range(len(stars)):
        x1, y1, x2, y2 = stars[i]
        print '  --- [%i] [%1.2f, %1.2f] - [%1.2f, %1.2f] :::  d = %1.2f'%(i, x1, y1, x2, y2, dist[i])
    return []
  elif len(goodstars) < (n/2):
    print 'WARNING: matched %i good stars out of %i!! '%(len(goodstars), n)
  return goodstars




#========== Class Star
#===============================================================

class Star:
  ''' Stores the magnitudes for each filter for a given star. The calibration must be done before! '''

  def SetDistance(self, d):
    ''' Sets the distance to the star, in parsec '''
    self.distance = d

  def SetMag(self, mag, filter = 'L'):
    ''' Sets the magnitude (has to be already calibrated!) for a filter '''
    if   filter == 'L': self.L = mag
    elif filter == 'R': self.R = mag
    elif filter == 'G': self.V = mag
    elif filter == 'V': self.V = mag
    elif filter == 'B': self.B = mag
    else: print '[Star.SetMag] WARNING: wrong filter (%s)'%filter

  def SetName(self, name):
    ''' Sets the name of the star... only for named stars or reference stars '''
    self.name = name

  def GetDistance(self):
    ''' Returns the distance of the star in parsec '''
    return self.distance

  def GetDistToPoint(self, x, y):
    ''' Returns the distance, in pixels, of the star to a given point in the image '''
    d2 = (x - self.GetX())*(x - self.GetX()) + (y - self.GetY())*(y - self.GetY())
    return np.sqrt(d2)

  def GetMag(self, filter = ''):
    ''' Returns the magnitude for a given filter '''
    if   filter == 'L': return self.L
    elif filter == 'R': return self.R
    elif filter == 'G': return self.V
    elif filter == 'V': return self.V
    elif filter == 'B': return self.B
    else: return self.V
    #mag = self.L if self.L != 0 else sefl.V
    #  return mag

  def GetCoor(self):
    ''' Returns the coordenates (in pixels) in the image '''
    return self.coor

  def GetName(self):
    ''' Returns the name of the star '''
    return self.name

  def GetT(self):
    ''' Returns the temperature from the B-V color index '''
    return GetTforBmV(self.GetBmV())

  def GetAbsMag(self, distance = -1):
    ''' Returns the absolute magnitude from the rel mag and the distance '''
    if distance != -1: self.SetDistance(distance)
    return GetAbsMagnitudeFromMagAndDistance(self.GetMag(), self.GetDistance())

  def GetBolometricMag(self):
    ''' Returns the bolometric magnitude from the abs magnitude and the temperature '''
    return GetBolometricMagnitude(self.GetAbsMag(), self.GetT())

  def GetLumi(self, distance = -1):
    ''' Returns the luminosity of the star from the bolometric magnitude and the distance ''' 
    if distance != -1: self.SetDistance(distance)
    return GetLumiFromMbol(self.GetBolometricMag())

  def GetRadius(self):
    ''' Returns the radius of the star using Stefan-Boltzmann law'''
    Lum = self.GetLumi()
    T = self.GetT()
    return GetStarRadius(Lum, T)

  def GetX(self):
    ''' Returns the X coordinate in the image (in pixels) '''
    return self.GetCoor()[0]

  def GetY(self):
    ''' Returns the Y coordinate in the image (in pixels) '''
    return self.GetCoor()[1]

  def GetBmV(self):
    ''' Returns the B - V color index '''
    return self.GetMag('B') - self.GetMag('V')

  def __init__(self, coor, R = -99, V = -99, B = -99, L = -99, distance = -1, name = ''):
    self.R = R
    self.V = V
    self.B = B
    self.L = L
    self.coor = coor
    self.distance = distance
    self.name = name
  
