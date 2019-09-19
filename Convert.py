import numpy as np
from PIL import Image # pillow package
from astropy.io import fits
import matplotlib.pyplot as plt
import rawpy
#from astropy.visualization import astropy_mpl_style
#plt.style.use(astropy_mpl_style)

def NewFitsFromData(data, filt = '', imgtyp = '', name = '', doSave = False, outpath = './'):
  # filt: 'R', 'G', 'B' (for imgtype LIGHT)
  # imgtyp: 'LIGHT', 'BIAS', 'FLAT', 'DARK'
  fit = fits.PrimaryHDU(data)
  if filt in ['B', 'R', 'G', 'L'] and imgtyp == '': imgtyp = 'LIGHT'
  if filt   != '': fit.header['FILTER'] = filt
  if imgtyp != '': fit.header['IMAGETYP'] = imgtyp

  if doSave: 
    if   filt == '' and imgtyp == '': fname = '%s.fit'%name
    elif filt == '' and imgtyp != '': fname = '%s_%s.fit'%(name,img)
    elif filt != '' and imgtyp == '': fname = '%s_%s.fit'%(name,filt)
    elif filt != '' and imgtyp in ['LIGHT', 'light', 'Light', 'L']: fname = '%s_%s.fit'%(name,filt)
    else:                             fname = '%s_%s_%s.fit'%(name,filt,img)
    fit.writeto(outpath + fname)
  return fit

def GetRGBfromRaw(imgname):
  raw = rawpy.imread(imgname)
  rgb = raw.postprocess(gamma=(1,1), no_auto_bright=True, output_bps=16)
  return rgb

def ConvertColorImage(imgname, imgtype = '', doSave = False, outpath = './'):
  image = Image.open(imgname)
  name = imgname.split('.')[0]
  ext = imgname.split('.')[-1]
  xsize, ysize = image.size
  print '>>>>>> Converting image %s (%s) %s to data...'%(name, ext, 'as "%s"'%imgtype if imgtype != '' else '')
  #plt.imshow(image)

  ##############################################################################
  # Split the three channels (RGB) and get the data as Numpy arrays.
  r, g, b = image.split()
  r_data = np.array(r.getdata()) # data is now an array of length ysize*xsize
  g_data = np.array(g.getdata())
  b_data = np.array(b.getdata())
  #print(r_data.shape)
  # Reshape the image arrays to be 2-dimensional:
  r_data = r_data.reshape(ysize, xsize)
  g_data = g_data.reshape(ysize, xsize)
  b_data = b_data.reshape(ysize, xsize)
  if   imgtype in ['dark', 'DARK', 'Dark', 'D', 'd']:
    fg = NewFitsFromData(g_data, filt='', imgtype = 'DARK', name=name, doSave=doSave, outpath=outpath)
  elif imgtype in ['flat', 'FLAT', 'Flat', 'F', 'f']:
    fg = NewFitsFromData(g_data, filt='', imgtype = 'FLAT', name=name, doSave=doSave, outpath=outpath)
  elif imgtype in ['bias', 'BIAS', 'Bias'          ]:
    fg = NewFitsFromData(g_data, filt='', imgtype = 'BIAS', name=name, doSave=doSave, outpath=outpath)
  else:
    fr = NewFitsFromData(r_data, filt='R', name=name, doSave=doSave, outpath=outpath)
    fg = NewFitsFromData(g_data, filt='G', name=name, doSave=doSave, outpath=outpath)
    fb = NewFitsFromData(b_data, filt='B', name=name, doSave=doSave, outpath=outpath)
    return [fr, fg, fb]
  return fg

def GetMetaData(image, printAll = False):
  ''' Get metadata dic from a JPG image... the image is open by PIL.Image... '''
  from PIL.ExifTags import TAGS, GPSTAGS
  if isinstance(image, str): image = Image.open(image)
  _TAGS_r = dict(((v, k) for k, v in TAGS.items()))
  _GPSTAGS_r = dict(((v, k) for k, v in GPSTAGS.items()))
  exifd = image._getexif() if hasattr(image, '_getexif') else image.tag
  keys = list(exifd.keys())
  if _TAGS_r["MakerNote"]   in keys: keys.remove(_TAGS_r["MakerNote"])
  if _TAGS_r["UserComment"] in keys: keys.remove(_TAGS_r["UserComment"])
  keys = [k for k in keys if k in TAGS]
  if printAll:
    print '################ Medatada:'
    print("\n".join([str((TAGS[k], exifd[k]))for k in keys]))
  iget = lambda x : exifd[_TAGS_r[x]]
  isk  = lambda x : _TAGS_r[x] in keys
  outdic = {}
  if isk('ExposureTime'):
    n,d = iget("ExposureTime")
    exptime = float(n)/d
    outdic['EXPTIME']  = exptime
  if   isk('DateTimeOriginal'): outdic['DATE-OBS'] = str(iget('DateTimeOriginal'))
  elif isk('DateTime'):         outdic['DATE-OBS'] = str(iget('DateTime'))
  if   isk('Make') and isk('Model'): outdic['INSTRUME'] = '%s - %s'%(str(iget('Make')), str(iget('Model')))
  return outdic 

def GetMetaDataRaw(image, printAll = False):
  ''' Get metadata from a raw image '''
  import exifread
  f = open(image)
  tags = exifread.process_file(f)
  outdic = {}
  outdic['EXPTIME']  = float(str(tags['EXIF ExposureTime']))
  outdic['DATE-OBS'] = str(tags['EXIF DateTimeOriginal'])
  outdic['INSTRUME'] = '%s - %s'%( str(tags['Image Model']), str(tags['Image Make']) ) 
  return outdic

#imgname = 'yopimage.jpg'
#imgname = '_MG_4817.JPG'
#imgname = '_MG_0296.CR2'
#outdic = GetMetaData(imgname,1)
#print outdic

