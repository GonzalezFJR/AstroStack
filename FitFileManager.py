import os
from math import sqrt, atan2, sin, cos, pi
import numpy as np

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

  def __init__(self, path, pathToDarks = '', pathToBias = '', pathToFlats = '', pathToRed = '', pathToGreen = '', pathToBlue = '', tag = ''):
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

#f = FitFileManager('fireworks', pathToDarks = 'darks')
#print f
#imgs = f.lights
