from math import sqrt, atan2, sin, cos, pi
import numpy as np
from StarMatch import *

def GetCentroid(Ax1, Ay1 = '', Ax2 = '', Ay2 = '', Bx1 = '', By1 = '', Bx2 = '', By2 = ''):
  if   Ay1 == '':
    return GetCentroid(Ax1[0], Ax1[1])
  elif Ax2 == '':
    star1 = Ax1;
    star2 = Ay1
    Ax1, Ay1, Ax2, Ay2 = star1
    Bx1, By1, Bx2, By2 = star2
  Axm = (Ax1 + Ax2)/2; Aym = (Ay1 + Ay2)/2
  Bxm = (Bx1 + Bx2)/2; Bym = (By1 + By2)/2
  tanA = (Ay2 - Ay1)/(Ax2 - Ax1); A = atan2(Ay2 - Ay1, Ax2 - Ax1)
  tanB = (By2 - By1)/(Bx2 - Bx1); B = atan2(By2 - By1, Bx2 - Bx1)
  sinA = sin(A); sinB = sin(B)
  cosA = cos(A); cosB = cos(B)
  Px = ( (Bym - Aym)*sinA*sinB + Bxm*sinA*cosB - Axm*cosA*sinB ) / sin(A - B)
  Py = ( (Axm - Bxm)*cosA*cosB + Aym*sinA*cosB - Bym*cosA*sinB ) / sin(A - B)
  return Px, Py

def RotatePoint(x, y, x0, y0, a):
  cosa = cos(a); sina = sin(a)
  xp = x0 + (x - x0)*cosa - (y - y0)*sina
  yp = y0 + (x - x0)*sina + (y - y0)*cosa
  return xp, yp

def RotateStarsAndGetMeanDist(matchedStars, Cx, Cy, angle, n = 10, verbose = False):
  distances = []
  istar = 0
  for star in matchedStars:
    if istar >= n: break
    x, y, xf, yf = star
    xp, yp = RotatePoint(x, y, Cx, Cy, angle)
    inid = GetDistance(x,  y,  xf, yf)
    d    = GetDistance(xp, yp, xf, yf)
    if verbose: print 'Rotation angle %1.3f: [%1.2f, %1.2f] --> [%1.2f, %1.2f] | target: [%1.2f, %1.2f]... init d = %1.2f,  current d = %1.2f'%(angle*180/pi, x, y, xp, yp, xf, yf, inid, d)
    distances.append(d)
    istar += 1
  return np.average(distances)

def GetOptimalRotationAngle(stars, centroid, nStars = 10, testAngle = 0.01, nSteps = 50000):
  angmin = 0; dmin = 10000
  Cx, Cy = centroid
  for ang in np.linspace(-testAngle, testAngle, nSteps):
    d = RotateStarsAndGetMeanDist(stars, Cx, Cy, ang, nStars)
    if d < dmin: 
      dmin   = d
      angmin = ang
  return angmin, dmin

'''
def GuessOptimalCentroidsAndRotAngle(stars, n = 8, ang = 0.001, steps = 50000):
  angmin = 0; dmin = 1000
  Cxmin = 0; Cymin = 0
  for i in range(n):
    for j in range(n):
      if j >= i: continue
      Cx, Cy = GetCentroid(stars[j], stars[i])
      ang, d = GetOptimalRotationAngle(stars, [Cx, Cy], n, 0.01, 50000)
      if d < dmin:
        dmin = d
        angmin = ang
        Cxmin = Cx; Cymin = Cy
      #print '[%1.0f, %1.0f], [angle, dmin] = [%1.8f, %1.4f]'%(Cx, Cy, ang, d)
  return dmin, angmin, Cxmin, Cymin
'''

def GetCenters(stars):
  CMx1 = 0; CMy1 = 0; CMx2 = 0; CMy2 = 0;
  n = len(stars)
  for star in stars:
    x1, y1, x2, y2 = star
    CMx1 += x1; CMy1 += y1
    CMx2 += x2; CMy2 += y2
  CMx1 = CMx1/n; CMy1 = CMy1/n
  CMx2 = CMx2/n; CMy2 = CMy2/n
  #print 'CM before: [%1.2f, %1.2f]'%(CMx1, CMy1)
  #print 'CM after : [%1.2f, %1.2f]'%(CMx2, CMy2)
  return CMx1, CMy1, CMx2, CMy2

def TranslateStars(stars, centers = ''):
  if   centers == '':
    CMx1, CMy1, CMx2, CMy2 = GetCenters(stars)
    OX = CMx1 - CMx2; OY = CMy1 - CMy2
  elif len(centers) == 2:
    OX, OY = centers
  elif len(centers) == 4:
    CMx1, CMy1, CMx2, CMy2 = centers
    OX = CMx1 - CMx2; OY = CMy1 - CMy2
  shiftedStars = []
  for s in stars:
    x1, y1, x2, y2 = s
    shiftedStars.append([x1, y1, x2+OX, y2+OY])
  return shiftedStars

def GetImgAlignment(refImg, img, maxstars = 100, testAngle = 0.0001, nSteps = 10000, verbose = 1):
  if isinstance(img, list):
    if refImg in img:
      index = img.index(refImg)
      img.pop(index)
    nImgs = len(img)
    outputs = {}
    i = 0
    for im in img:
      i += 1
      if verbose >= 0: print '[IMG %i / %i]'%(i, nImgs)
      ofname, ocentr, oshift, oang, odmin  = GetImgAlignment(refImg, im, maxstars, testAngle, nSteps, verbose)
      outputs[ofname] = [ocentr, oshift, oang, odmin]
    return outputs
  if verbose: print ' ## Matching stars in files "%s" and "%s"...'%(refImg, img)
  matchedStars = GetGoodMatchedStars(refImg, img, n = maxstars)
  if len(matchedStars) == 0:
    print '0 stars -- Something went worng... skipping image %s'%img
    return [img, 0, 0, 0, 0]
  if verbose: print '    >> Using %i stars!'%len(matchedStars)
  dOrig  = GetMeanStarsDistance(matchedStars)
  if verbose: print '    >> Img orig distance: %1.2f'%dOrig
  if verbose: print ' ## Obtaining centroids...'
  CMx1, CMy1, CMx2, CMy2 = GetCenters(matchedStars)
  OX = CMx1 - CMx2; OY = CMy1 - CMy2
  if verbose: print '    >> [%1.2f, %1.2f] and [%1.2f, %1.2f]... shifts: [%1.2f, %1.2f]'%(CMx1, CMy1, CMx2, CMy2, OX, OY)
  shiftedStars = TranslateStars(matchedStars)
  dTrans = GetMeanStarsDistance(shiftedStars)
  if verbose: print '    >> Distance after translation: %1.2f '%dTrans
  if (dOrig-dTrans)/dOrig < 0.5: 
    print 'WARNING: probably a bad translation... consider increasing the number of stars!'
  if verbose: print ' ## Opimazing rotation angle...'
  angmin, dmin = GetOptimalRotationAngle(shiftedStars, [CMx1, CMy1], nStars = len(shiftedStars), testAngle=testAngle, nSteps=nSteps)
  if verbose: print '    >> Rotating an angle of %1.5f (%1.5f dgr)'%(angmin, 180*angmin/pi)
  if verbose: print '    >> Distance after rotation: %1.2f '%dmin
  if verbose > 1: RotateStarsAndGetMeanDist(shiftedStars, CMx1, CMy1, angmin, len(shiftedStars), 1)
  return [img, [CMx1, CMy1], [OX, OY], angmin, dmin]

def ShiftAndRotate(fname, center, shift, angle, outname = ''):
  loadimg = LoadFitImage(path)
  image = loadimg[0]
  shifted = shift(image, shift)
  padX = [shifted.shape[1] - center[0], center[0]]
  padY = [shifted.shape[0] - center[1], center[1]]
  imgdesp  = np.pad(shifted, [padY, padX], 'constant')
  imgrot   = rotate(imgdesp, angle, reshape=False)
  imgfinal = imgrot[padY[0] : -padY[1], padX[0] : -padX[1]]
  if outname == '':
    extension = ''
    if   fname.endswith('.fit' ): 
      outname = fname[:-4]
      extension = '.fit'
    elif fname.endswith('.fits'): 
      outname = fname[:-5]
      extension = '.fits'
    outname += '_align'+extension
  if os.path.isfile(outname):
    os.system('mv %s %s.bak'%outname)
  f = LoadFit(path)
  f[0].data = imgfinal
  f.writeto(outname,overwrite=True)


