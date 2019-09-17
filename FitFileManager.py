import os, random
from math import sqrt, atan2, sin, cos, pi
import numpy as np
from FitImage import FitImage

class FitFileManager:
  def SetPath(self, path):
    self.path = path

  def SetImgName(self, img):
    self.name = img

  def SetFileTag(self, nam = ''):
    self.tag = nam

  def SetFileEndName(self, end='.fit'):
    self.endname = end

  def SetTagDark(self, t = 'D'):
    self.tagDark = t

  def SetTagBias(self, t = 'I'):
    self.tagBias = t

  def SetTagFlat(self, t = 'F'):
    self.tagFlat = t

  def SetTagBlue(self, t = 'B'):
    self.tagBlue = t

  def SetTagRed(self, t = 'R'):
    self.tagRed = t

  def SetTagGreen(self, t = 'B'):
    self.tagGreen = t

  def SetTagLight(self, t = ''):
    self.tagLight= t

  def SetPathDark(self, t = ''):
    self.pathDark = t

  def SetPathBias(self, t = ''):
    self.pathBias = t

  def SetPathFlat(self, t = ''):
    self.pathFlat = t

  def SetPathBlue(self, t = ''):
    self.pathBlue = t

  def SetPathRed(self, t = ''):
    self.pathRed = t

  def SetPathGreen(self, t = ''):
    self.pathGreen = t

  def Init(self, path = '', ftag = '', pdarks = '', pflats = '', pbias = '', pblue = '', pred = '', pgreen = ''):
    self.SetPath(path)
    self.SetFileTag(ftag)
    self.SetFileEndName()
    self.SetTagDark()
    self.SetTagFlat()
    self.SetTagBias()
    self.SetTagBlue()
    self.SetTagRed()
    self.SetTagGreen()
    self.SetTagLight()
    self.SetPathDark(pdarks)
    self.SetPathFlat(pflats)
    self.SetPathBias(pbias)
    self.SetPathBlue(pblue)
    self.SetPathRed(pred)
    self.SetPathGreen(pgreen)

  def IsGoodFile(self, f):
    if not f.endswith(self.endname): return False
    if self.tag != '' and not self.tag in f: return False 
    return True

  def IsDark(self, f):
    if not self.IsGoodFile(f): return False
    if ('_'+self.tagDark+'_') in f or ('_'+self.tagDark+'.') in f or 'dark' in f or 'Dark' in f or 'DARK' in f: return True
    return False

  def IsFlat(self, f):
    if not self.IsGoodFile(f): return False
    if ('_'+self.tagFlat+'_') in f or ('_'+self.tagFlat+'.') in f or 'flat' in f or 'Flat' in f or 'FLAT' in f: return True
    return False

  def IsBias(self, f):
    if not self.IsGoodFile(f): return False
    if ('_'+self.tagBias+'_') in f or ('_'+self.tagBias+'.') in f or 'bias' in f or 'Bias' in f or 'BIAS' in f: return True
    return False

  def IsBlue(self, f):
    if not self.IsGoodFile(f): return False
    if ('_'+self.tagBlue+'_') in f or ('_'+self.tagBlue+'.') in f or 'blue' in f or 'Blue' in f or 'BLUE' in f: return True
    return False

  def IsGreen(self, f):
    if not self.IsGoodFile(f): return False
    if ('_'+self.tagGreen+'_') in f or ('_'+self.tagGreen+'.') in f or 'green' in f or 'Green' in f or 'GREEN' in f: return True
    return False

  def IsRed(self, f):
    if not self.IsGoodFile(f): return False
    if ('_'+self.tagRed+'_') in f or ('_'+self.tagRed+'.') in f or 'red' in f or 'Red' in f or 'RED' in f: return True
    return False

  def IsLight(self, f):
    if not self.IsGoodFile(f): return False
    if self.tagLight != '' and not self.tagLight in f: return False 
    if self.IsDark(f) or self.IsBias(f) or self.IsFlat(f) or self.IsBlue(f) or self.IsGreen(f) or self.IsRed(f): return False
    return True

  def Read(self, path, tag = ''):
    for f in os.listdir(path):
      fname = '%s/%s'%(path, f)
      if os.path.isfile(fname):
         if not self.IsGoodFile(fname): continue
         if   self.IsDark( fname): self.darks .append(fname)
         elif self.IsFlat( fname): self.flats .append(fname)
         elif self.IsBias( fname): self.bias  .append(fname)
         elif self.IsRed(  fname): self.red   .append(fname)
         elif self.IsGreen(fname): self.green .append(fname)
         elif self.IsBlue( fname): self.blue  .append(fname)
         elif self.IsLight(fname): self.lights.append(fname)
      elif os.path.isdir(fname):
        self.Read(fname)

  def splitColors(self):
    return len(self.blue) > 0 or len(self.green) > 0 or len(self.red) > 0

  def __init__(self, path, pathToDarks = '', pathToBias = '', pathToFlats = '', pathToRed = '', pathToGreen = '', pathToBlue = '', tag = '', verbose = False):
    self.Init(path, tag, pathToDarks, pathToFlats, pathToBias, pathToBlue, pathToRed, pathToGreen)
    self.lights = [] 
    self.darks = []
    self.bias = []
    self.flats = []
    self.blue = []
    self.red = []
    self.green = []
    self.Read(self.path)
    if pathToDarks != '': self.Read(pathToDarks, 'dark')
    if pathToFlats != '': self.Read(pathToFlats, 'flat')
    if pathToBias  != '': self.Read(pathToBias,  'bias')
    if pathToRed   != '': self.Read(pathToRed,   'red')
    if pathToBlue  != '': self.Read(pathToBlue,  'blue')
    if pathToGreen != '': self.Read(pathToGreen, 'green')
    self.verbose = verbose
    if verbose: print self

  def __str__(self):
    c = '[FitFileManager]\n'
    if len(self.lights) > 0:
      c += 'LIGHTS (%i)\n'%len(self.lights)
      for l in self.lights: c += ' >> %s\n'%l
    if len(self.darks) > 0:
      c += 'DARKS (%i)\n'%len(self.darks)
      for l in self.darks: c += ' >> %s\n'%l
    if len(self.flats) > 0:
      c += 'FLATS (%i)\n'%len(self.flats)
      for l in self.flats : c += ' >> %s\n'%l
    if len(self.bias) > 0:
      c += 'BIAS (%i)\n'%len(self.bias)
      for l in self.bias: c += ' >> %s\n'%l
    if len(self.red) > 0:
      c += 'RED (%i)\n'%len(self.red)
      for l in self.red: c += ' >> %s\n'%l
    if len(self.green) > 0:
      c += 'GREEN (%i)\n'%len(self.green)
      for l in self.green: c += ' >> %s\n'%l
    if len(self.blue) > 0:
      c += 'BLUE (%i)\n'%len(self.blue)
      for l in self.blue: c += ' >> %s\n'%l
    return c


########################################################################################
########################################################################################
########################################################################################
class ImgManager:

  def report(self, t, verbose = -1):
    vb = verbose if verbose != -1 else self.verbose
    self.reportlist.append(t+'\n')
    if vb > 0: print t

  def SetOutFolder(self,s):
    self.outfolder = s

  def SetImageName(self, n):
    self.imgname = n

  def SetHeader(self, h):
    self.header = h

  def GetLightImg(self, imgname):
    for img in (self.lights + self.blue + self.red + self.green):
      if img.name == imgname: return img
    return None

  def SetRefImage(self, imgname = ''):
    if self.GetLightImg(imgname) != None: self.imgref = imgname
    else:
      if   len(self.lights) > 0 : self.imgref = self.lights[0].name
      elif len(self.red)   > 0 : self.imgref = self.red[0]   .name
      elif len(self.green) > 0 : self.imgref = self.green[0] .name
      elif len(self.blue)  > 0 : self.imgref = self.blue[0]  .name
      else: self.imgref = ''

  def Load(self, imgname, imgtype = 'light'):
    if isinstance(imgname, str) and ',' in imgname:
      imgnamelist = imgname.replace(' ', '').split(',')
      self.Load(imgnamelist, imgtype)
      return
    if isinstance(imgname, list):
      for iname in imgname: self.Load(iname, imgtype)
      return
    self.report('Loading image %s as %s'%(imgname, imgtype), self.verbose)
    img = FitImage(imgname, imgtype, verbose = 1); 
    img.Load()
    img.SetHeaderInfo(self.header)
    img.SetFWHM(self.FWHM)
    img.CalcBkg()
    img.CalcMinMaxVals()
    img.ReportData()
    img.ReportHeader()
    filt   = img.GetFilterHeader()
    imtype = img.GetImgTypeHeader()
    if   imtype == 'LIGHT': 
      if   filt == 'L': self.lights.append(img)
      elif filt == 'B': self.blue.append(img)
      elif filt == 'R': self.red.append(img)
      elif filt == 'G': self.green.append(img)
      else:
        print 'WARNING: unknown filter %s'%filt
    elif imtype == 'DARK': self.darks.append(img)
    elif imtype == 'FLAT': self.flats.append(img)
    elif imtype == 'BIAS': self.bias.append(img)
    else:
      print 'WARNING: unknown image type %s!! Adding to lights...'%imtype
      self.lights.append(img)

  ################################################################
  ### Calibration

  def Calibrate(self):
    # Get master bias
    self.biasMaster = self.GetBiasMaster()
    # Apply bias correction to darks and flats
    if self.biasMaster != None: 
      for img in (self.darks + self.flats): img.CorrectBias(self.biasMaster)
    # Look for hot and death pixels in darks and flats:
    for img in (self.darks + self.flats):
      img.LookForHotPixels()
      img.CorrectHotPixels()
    # Get master dark and flat
    self.darkMaster = self.GetDarkMaster()
    self.flatMaster = self.GetFlatMaster()
    # Correct light images
    for img in (self.lights + self.blue + self.red + self.green):
      if self.correctHotPixels:
        img.LookForHotPixels()
        img.CorrectHotPixels()
      if self.biasMaster != None: img.CorrectBias(self.biasMaster)
      if self.darkMaster != None: img.CorrectDark(self.darkMaster)
      if self.flatMaster != None: img.CorrectFlat(self.flatMaster)
    
  def GetBiasMaster(self):
    if len(self.bias) == 0: return None
    self.biasMaster = self.bias[0]
    if len(self.bias) > 1:
      moreBias = self.bias[1:]
      self.biasMaster.Stack(moreBias, 'mean')
    self.darkMaster.Save(self.outfolder, 'bias_master')
    return self.biasMaster

  def GetDarkMaster(self):
    if len(self.darks) == 0: return None
    times = [x.GetExpTime() for x in self.darks]
    mediantime = np.median(times)
    self.darkMaster = self.darks[0]
    if len(self.darks) > 1:
      moreDarks = self.darks[1:]
      self.darkMaster.Stack(moreDarks, 'mean', scaleTo = mediantime)
    self.darkMaster.Save(self.outfolder, 'dark_master')
    return self.darkMaster

  def GetFlatMaster(self):
    if len(self.flats) == 0: return None
    self.flatMaster = self.flats[0]
    if len(self.flats) > 1:
      moreFlats = self.flats[1:]
      self.flatMaster.Stack(moreFlats, 'mean')
    self.darkMaster.Save(self.outfolder, 'flat_master')
    return self.flatMaster
    
  ##############################################################################
  ### Align and stack
  def SetAlgParams(self, alg_distthr, alg_maxStars, alg_testAngle, alg_nSteps):
    if alg_distthr   > 0: self.alg_distthr   = alg_distthr
    if alg_maxStars  > 0: self.alg_maxStars  = alg_maxStars
    if alg_testAngle > 0: self.alg_testAngle = alg_testAngle
    if alg_nSteps    > 0: self.alg_nSteps    = alg_nSteps

  def Align(self, refImageName = '', distthr = -1, maxstars = -1, testAngle = -1, nSteps = -1):
    self.SetAlgParams(distthr, maxstars, testAngle, nSteps)
    self.SetRefImage(refImageName)
    refimg = self.GetLightImg(self.imgref)
    refimg.alignStatus = True
    repeatAlg = []; usedRef = [refimg.name]
    for img in (self.lights + self.blue + self.red + self.green):
      if img.name == refimg.name: continue
      img.GetImgAlignment(refimg, self.alg_maxStars, self.alg_testAngle, self.alg_nSteps)
      if img.alg_dist > self.alg_distthr: repeatAlg.append(img.name)
      else: img.Align()
    if len(repeatAlg) != 0:
      self.report('WARNING!! Problems with alignment for the following images: \n' + str(repeatAlg),1)
      tries = 1
      while tries <= 5 and len(repeatAlg) > 0:
        tries += 1
        self.report(' >> Repeating alignment... try %i'%tries,1)
        # Get new ref image
        listOfStars = self.lights + self.blue + self.red + self.green; random.shuffle(listOfStars)
        for rimg in listOfStars:
          if rimg.name not in (repeatAlg + usedRef):
            usedRef.append(rimg.name)  
            break
        # Try with a new image
        testImages = (self.GetLightImg(x) for x in repeatAlg)
        repeatAlg = []
        for img in testImages:
          img.GetImgAlignment(rimg, self.alg_maxStars, self.alg_testAngle, self.alg_nSteps)
          if img.alg_dist > self.alg_distthr: repeatAlg.append(img.name)
          else: img.Align()
    if len(repeatAlg) != 0:
      self.report('WARNING: could not align the following images: \n'+ str(repeatAlg),1)
              
  def Stack(self, mode = 'sum'):
    self.report("### Begin stack with method: %s"%mode)
    glights = filter(lambda x : x.alignStatus, self.lights)
    gred    = filter(lambda x : x.alignStatus, self.red)
    ggreen  = filter(lambda x : x.alignStatus, self.green)
    gblue   = filter(lambda x : x.alignStatus, self.blue)
    if len(glights) > 0:
      self.report('Stacking %i light images...'%len(glights))
      self.lightStack = glights[0] 
      self.lightStack.Stack(glights, mode)
      self.lightStack.Save(self.outfolder, 'light_%s'%mode)
    if len(gred) > 0:
      self.report('Stacking %i red images...'%len(gred))
      self.redStack = gred[0] 
      self.redStack.Stack(gred, mode)
      self.redStack.Save(self.outfolder, 'red_%s'%mode)
    if len(ggreen) > 0:
      self.report('Stacking %i green images...'%len(ggreen))
      self.greenStack = ggreen[0] 
      self.greenStack.Stack(ggreen, mode)
      self.greenStack.Save(self.outfolder, 'green_%s'%mode)
    if len(gblue) > 0:
      self.report('Stacking %i blue images...'%len(gblue))
      self.blueStack = gblue[0] 
      self.blueStack.Stack(gblue, mode)
      self.blueStack.Save(self.outfolder, 'blue_%s'%mode)
              
  ##############################################################################
  ### Init
    
  def __init__(self, fileManager = None, outfolder = './project', imgname = '', header = {}, FWHM = 5, refImg = '', lights = [], darks = [], flats = [], bias = [], red = [], green = [], blue = [], saveAllSteps = True, verbose = 0):
    self.lights = [] 
    self.darks = []
    self.bias = []
    self.flats = []
    self.blue = []
    self.red = []
    self.green = []
    self.verbose = verbose
    self.SetAlgParams(1, 100, 0.0001, 10000)
    self.FWHM = FWHM
    self.correctHotPixels = True
    self.reportlist = []
    self.SetOutFolder(outfolder)
    self.SetImageName(imgname)
    self.SetHeader(header)
    self.saveAllSteps = saveAllSteps
    self.Load(lights, 'light')
    self.Load(darks,  'dark')
    self.Load(flats,  'flat')
    self.Load(bias,   'bias')
    self.Load(red,    'red')
    self.Load(green,  'green')
    self.Load(blue,   'blue')
    if fileManager != None:
      self.Load(fileManager.lights, 'light')
      self.Load(fileManager.darks,  'dark')
      self.Load(fileManager.flats,  'flat')
      self.Load(fileManager.bias,   'bias')
      self.Load(fileManager.red,    'red')
      self.Load(fileManager.green,  'green')
      self.Load(fileManager.blue,   'blue')

#f = FitFileManager('fireworks', pathToDarks = 'darks')
#print f
#imgs = f.lights
