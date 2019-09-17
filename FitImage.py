from StarAlign import *
from math import sqrt, atan2, sin, cos, pi
import numpy as np
from StarAlign import GetImgAlignment, ShiftAndRotate
from scipy.ndimage.interpolation import shift, rotate
from astropy.io import fits

#path = 'fireworks/fireworks_002.fit'
def fitw(w, n=20): 
  while len(w) < n: w += ' '
  return w

class FitImage:
  def SetName(self, n):
    self.name = n

  def SetType(self, im):
    self.type = im

  def SetImg(self, im):
    self.img = im

  def IsBiasCorr(self, s = True):
    self.kIsBiasCorr = s

  def IsDarkCorr(self, s = True):
    self.kIsBiasCorr = s

  def IsFlatCorr(self, s = True):
    self.kIsBiasCorr = s

  def AddImgToList(self, img):
    self.imglist.append(img)

  def Load(self):
    if self.img != None: self.report('WARNING: Reloading the image...', 1)
    if not os.path.isfile(self.name):
      print 'ERROR: not found file ' + self.name
      return 
    self.report('Loading image... '+self.name)
    self.fits = ffits = astropy.io.fits.open(self.name)
    self.hdu = ffits[0]
    self.img = self.hdu.data.astype(float)
    if self.img.max() <= 32767 or self.img.min() < 0:
      self.img += 32768
    self.SetHeaderFromImgtype()

  def report(self, t, verbose = -1):
    vb = verbose if verbose != -1 else self.verbose
    self.reportlist.append(t+'\n')
    if vb > 0: print t

  def Save(self, path, outname, overwrite = False):
    if not os.path.isdir(path): os.mkdir(path)
    self.fits[0].data = self.img
    if not outname.endswith('.fit') and not outname.endswith('.fits'): outname += '.fit'
    str(self.fits[0].header)
    fname = path + '/' + outname
    if os.path.isfile(fname):
      lname = fname.split('.')
      extension = lname[-1]
      newname = ''; 
      for l in lname[:-1]: newname+=l
      newname += '_old' + extension
      if fname.startswith('.'): newname = '.'+newname
      os.rename(fname, newname)
    self.fits.writeto(path + '/' + outname, overwrite=overwrite)

  #########################################################
  ## Header

  def GetHeaderInfo(self):
    keys = self.hdu.header.keys()
    t = ''
    if 'OBSERVER' in keys : t += '  ...%s : %s\n'%(fitw('OBSERVER'), self.hdu.header['OBSERVER'])
    if 'OBJECT'   in keys : t += '  ...%s : %s\n'%(fitw('OBJECT'),   self.hdu.header['OBJECT'])
    if 'DATE-OBS' in keys : t += '  ...%s : %s\n'%(fitw('DATE-OBS'), self.hdu.header['DATE-OBS'])
    if 'EXPTIME'  in keys : t += '  ...%s : %s\n'%(fitw('EXPTIME'),  self.hdu.header['EXPTIME'])
    if 'FILTER'   in keys : t += '  ...%s : %s\n'%(fitw('FILTER'),   self.hdu.header['FILTER'])
    if 'IMAGETYP' in keys : t += '  ...%s : %s\n'%(fitw('IMAGETYP'), self.hdu.header['IMAGETYP'])
    if 'INSTRUME' in keys : t += '  ...%s : %s\n'%(fitw('INSTRUME'), self.hdu.header['INSTRUME'])
    if 'SWCREATE' in keys : t += '  ...%s : %s\n'%(fitw('SWCREATE'), self.hdu.header['SWCREATE'])
    if 'BITPIX'   in keys : t += '  ...%s : %s\n'%(fitw('BITPIX'),   self.hdu.header['BITPIX'])
    if 'NAXIS1'   in keys : t += '  ...%s : %s\n'%(fitw('NAXIS1'),   self.hdu.header['NAXIS1'])
    if 'NAXIS2'   in keys : t += '  ...%s : %s\n'%(fitw('NAXIS2'),   self.hdu.header['NAXIS2'])
    if 'XPIXSZ'   in keys : t += '  ...%s : %s\n'%(fitw('XPIXSZ'),   self.hdu.header['XPIXSZ'])
    if 'YPIXSZ'   in keys : t += '  ...%s : %s\n'%(fitw('YPIXSZ'),   self.hdu.header['YPIXSZ'])
    if 'XBINNING' in keys : t += '  ...%s : %s\n'%(fitw('XBINNING'), self.hdu.header['XBINNING'])
    if 'YBINNING' in keys : t += '  ...%s : %s\n'%(fitw('YBINNING'), self.hdu.header['YBINNING'])
    return t

  def SetHeaderInfo(self, dic):
    for k in dic.keys():
      self.SetHeaderParam(k, dic[k])

  def SetHeaderParam(self, param, val):
    self.hdu.header[param] = val

  def GetExpTime(self):
    return float(self.hdu.header['EXPTIME'])

  def GetFilterHeader(self):
    return self.hdu.header['FILTER'] if 'FILTER' in self.hdu.header.keys() else ''

  def GetImgTypeHeader(self):
    return self.hdu.header['IMAGETYP'] if 'IMAGETYP' in self.hdu.header.keys() else ''

  def ReportHeader(self):
    t  = '\n#################################################\n'
    t += '### Header info\n'
    t += self.GetHeaderInfo()
    self.report(t)

  def SetHeaderFromImgtype(self):
    headimtype = self.GetImgTypeHeader()
    filt       = self.GetFilterHeader()
    if   self.type == 'light': 
      if headimtype != '' and headimtype != 'LIGHT': print 'WARNING: possible type conflict "%s" and "%s"... setting "%s"'%(self.type, headimtype, 'LIGHT')
      self.SetImgTypeHeader('LIGHT')
      if   filt in ['R', 'G', 'B']:
        print 'WARNING: setting light image "%s" to filter "%s" (based on header info)'%(self.name, filt)
        self.SetFilterHeader(filt)
      else: 
        self.SetFilterHeader('L')
    elif self.type == 'dark' : 
      if headimtype != '' and headimtype != 'DARK': print 'WARNING: possible type conflict "%s" and "%s"... setting "%s"'%(self.type, headimtype, 'DARK')
      self.SetImgTypeHeader('DARK')
    elif self.type == 'bias' : 
      if headimtype != '' and headimtype != 'BIAS': print 'WARNING: possible type conflict "%s" and "%s"... setting "%s"'%(self.type, headimtype, 'BIAS')
      self.SetImgTypeHeader('BIAS')
    elif self.type == 'flat' : 
      if headimtype != '' and headimtype != 'FLAT': print 'WARNING: possible type conflict "%s" and "%s"... setting "%s"'%(self.type, headimtype, 'FLAR')
      self.SetImgTypeHeader('FLAT')
    elif self.type == 'blue' : 
      if filt != '' and filt != 'B': print 'WARNING: possible color conflict "%s" and "%s"... setting "%s"'%(self.type, filt, 'B')
      self.SetImgTypeHeader('LIGHT')
      self.SetFilterHeader('B')
    elif self.type == 'red'  : 
      if filt != '' and filt != 'R': print 'WARNING: possible color conflict "%s" and "%s"... setting "%s"'%(self.type, filt, 'R')
      self.SetImgTypeHeader('LIGHT')
      self.SetFilterHeader('R')
    elif self.type == 'green': 
      if filt != '' and filt != 'G': print 'WARNING: possible color conflict "%s" and "%s"... setting "%s"'%(self.type, filt, 'G')
      self.SetImgTypeHeader('LIGHT')
      self.SetFilterHeader('G')

  def SetFilterHeader(self, val):
    self.SetHeaderParam('FILTER', val)

  def SetImgTypeHeader(self, val):
    self.SetHeaderParam('IMAGETYP', val)

  def SetExpTimeHeader(self, val):
    self.SetHeaderParam('EXPTIME', val)

  #########################################################
  ## More data

  def resetData(self):
    self.fwhm = -1
    self.bkg  = 0
    self.bkg_sigma = 0
    self.maxval = -1
    self.minval = -1
    self.exptime = 0

  def SetFWHM(self, f):
    self.fwhm = f

  def CalcBkg(self):
    if not isinstance(self.img, np.ndarray): return
    from astropy.stats import mad_std
    bkg   = np.median(self.img)
    sigma = mad_std(self.img)
    self.bkg       = bkg
    self.bkg_sigma = sigma

  def CalcMinMaxVals(self):
    if not isinstance(self.img, np.ndarray): return
    self.maxval = self.img.max()
    self.minval = self.img.min()

  def ReportData(self):
    self.report('\n#################################################')
    self.report('### Image data')
    self.report(' >> Background      : %1.2f'%self.bkg,       self.verbose)
    self.report(' >> Background sigma: %1.2f'%self.bkg_sigma, self.verbose)
    self.report(' >> Minimum value   : %1.2f'%self.minval,    self.verbose)
    self.report(' >> Maximum value   : %1.2f'%self.maxval,    self.verbose)

  #########################################################
  ## Aligns
  def GetImgAlignment(self, refImg, maxstars = 100, testAngle = 0.0001, nSteps = 10000):
    self.alg_maxstars = maxstars
    self.alg_testAngle = testAngle
    self.alg_nSteps = nSteps
    self.report('## Start alignment for images:')
    self.report('## Reference: %s'%refImg.name)
    self.report('## Target   : %s'%self.name)
    imname, center, shift, angle, dist, text = GetImgAlignment(refImg.img, self.img, self.alg_maxstars, self.alg_testAngle, self.alg_nSteps, self.verbose)
    self.SetAlignParams(center, shift, angle, dist)
    self.alg_text = text

  def SetAlignParams(self, center = [], shift = [], angle = 0, dist = -1):
    self.alg_center = center
    self.alg_shift  = shift
    self.alg_angle  = angle
    self.alg_dist   = dist

  def Align(self, center = [], shift = [], angle = 0, dist = -1):
    if center != [] : self.alg_center = center
    if shift  != [] : self.alg_shift  = shift
    if angle  !=  0 : self.alg_angle  = angle
    if dist   != -1 : self.alg_dist   = dist
    self.img = ShiftAndRotateImage(self.img, self.alg_center, self.alg_shift, self.alg_angle)
    self.alignStatus = True

  def ReportAlign(self):
    self.report('\n#################################################\n')
    self.report('### Align data for img %s'%self.name)
    self.report(' >> Ref image: %s'%self.imgrefname)
    self.report(self.alg_text)
    self.report(' >> Center   : [%1.1f, %1.1f]'%(self.alg_center[0], self.alg_center[1]))
    self.report(' >> Shift    : [%1.1f, %1.1f]'%(self.alg_shift[0], self.alg_shift[1]))
    self.report(' >> Angle   : %1.2f'%self.alg_angle)
    self.report(' >> Mean dist: %1.2f'%self.alg_dist)

  #########################################################
  ## Stack

  def Stack(self, imglist, mode = 'sum', scaleTo = -1):
    ''' Stack modes: sum, median, average '''
    self.stackmode = mode
    self.stackScaleTo = scaleTo
    self.imglist = [i.name for i in imglist]
    gimg = []
    lx   = len(self.img)
    ly   = len(self.img[0])
    for img in imglist:
      lenx = len(img.img)
      leny = len(img.img[0])
      if not lenx == lx or not leny == ly:
        print 'ERROR: stacking images with different dimenstions!!! %ix%i and %ix%i'%(lx,ly,lenx,leny)
        continue
      gimg.append(img)
    fact = 1
    if mode == 'sum':
      for img in gimg: 
        if scaleTo != -1:
          fact = scaleTo/img.GetExpTime()
        self.img += (img.img*fact)
    elif mode == 'average' or mode == 'mean':
      n = len(gimg)
      for img in gimg: 
        if scaleTo != -1: fact = scaleTo/img.GetExpTime()
        self.img += (img.img*fact)
      self.img *= 1./float(n+1)
    elif mode == 'median':
      for i in range(lx):
        for j in range(ly):
          p = [self.img[i][j]]
          for img in gimg: p.append(img.img[i][j])
          self.img[i][j] = np.median(p)
    else: self.report("ERROR: unknown stacking mode %s!!!"%mode, 10)

  def ReportStack(self):
    t  = '\n#################################################\n'
    t += '### Stack info\n'
    i = 0
    t += ' ++ Stack mode "%s"\n'%self.stackmode
    if self.stackScaleTo != -1: t += ' ++ Scaling to "%1.1f seconds..."\n'%self.stackScaleTo
    t += ' ++ Total: %i images stacked! \n'%(len(self.imglist)+1)
    for iname in self.imglist:
      t += ' >> Stacked image [%i]: %s\n'%(i,iname)
      i += 1
    self.report(t)
 
  #########################################################
  ## Calibrate

  def CorrectBias(self, biasImg):
    self.biasImg = biasImg.name
    self.img -= biasImg.img
    self.IsBiasCorr(True)

  def CorrectDark(self, darkImg):
    self.darkImg = darkImg.name
    timecorr = self.GetExpTime()/darkImg.GetExpTime()
    self.img -= (darkImg.img*timecorr)
    lx   = len(self.img)
    ly   = len(self.img[0])
    for i in range(lx):
      for j in range(ly):
        if self.img[i][j] < 0: self.img[i][j] = 0
    self.IsDarkCorr(True)

  def CorrectFlat(self, flatImg):
    self.flatImg = flatImg.name
    flatMean = flatImg.img.mean()
    self.img = self.img/flatImg.img*flatMean
    self.IsFlatCorr(True)

  def LookForHotPixels(self):
    self.hotPixels = []
    self.deathPixels = []
    for i in range(len(self.img)):
      raw = self.img[i]
      for j in range(len(raw)):
        v = self.img[i][j]
        if v == self.maxval:
          self.hotPixels.append([i,j])
        elif v == 0:
          self.deathPixels.append([i,j])

  def InterpolatePixel(self, x, y):
    vals = []; gvals = []
    maxx = len(self.img)
    maxy = len(self.img[0])
    if y+1 < maxy: vals.append(self.img[x  ][y+1])
    if x+1 < maxx: vals.append(self.img[x+1][y  ])
    if y   > 0   : vals.append(self.img[x  ][y-1])
    if x   > 0   : vals.append(self.img[x-1][y  ])
    for v in vals:
      if v != 0 and v != self.maxval: gvals.append(v)
    mean = float(sum(gvals))/len(gvals)
    return int(mean)

  def CorrectHotPixels(self):
    if not hasattr(self, 'hotPixels'): self.LookForHotPixels()
    for pix in self.hotPixels:
      x, y = pix
      val = self.InterpolatePixel(x, y)
      self.img[x][y] = val
    for pix in self.deathPixels:
      x, y = pix
      val = self.InterpolatePixel(x, y)
      self.img[x][y] = val
    self.report('  >> Pixel corrections applied using pixel interpolation! [%i hot pixels] [%i death pixels]'%(len(self.hotPixels),len(self.deathPixels)))

  def ReportCalibration(self):
    t  = '\n#################################################\n'
    t += '### Calibration info\n'
    t += ' >> Bias correction: %s\n'%self.biasImg if self.kIsBiasCorr else ''
    t += ' >> Dark correction: %s\n'%self.darkImg if self.kIsDarkCorr else ''
    t += ' >> Flat correction: %s\n'%self.flatImg if self.kIsFlatCorr else ''
    if hasattr(self, 'deathPixels'):
      t += ' >> Death pixels   : \n%s'%str(self.deathPixels)
    if hasattr(self, 'hotPixels'):
      t += ' >> Hot   pixels   : \n%s'%str(self.hotPixels)
    self.report(t)

  #########################################################
  ## Init
 
  def __init__(self, name, imgtype = 'light', img = None, verbose = 0):
    self.SetName(name)
    self.SetImg(img)
    self.SetType(imgtype)
    self.IsBiasCorr(False)
    self.IsDarkCorr(False)
    self.IsFlatCorr(False)
    self.imglist = []
    self.AddImgToList(name)
    self.SetAlignParams(); 
    self.reportlist = []
    self.verbose = verbose
    self.alignStatus = False

'''
f = FitFileManager('fireworks', pathToDarks = 'darks', verbose = True)
imgs = f.lights
out = GetImgAlignment(f.lights[0], f.lights, maxstars = 100, testAngle = 0.001, nSteps = 10000, verbose = 1)
out = ShiftAndRotateGroup(out, dthr = 1)

for name in out.keys():
  center, shift, angle, d = out[name]
  ShiftAndRotate(name, center, shift, angle)


def ShiftAndRotateGroup(dic, dthr = 1):
  names = dic.keys()
  out = {}; repeat = []; done = []
  for name in names:
    center, shift, angle, d = dic[name]
    if d < dthr:
      img = ShiftAndRotate(name, center, shift, angle)
      rotinfo = [center, shift, angle, d]

'''

