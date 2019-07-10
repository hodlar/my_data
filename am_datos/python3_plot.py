import numpy as np
import matplotlib.pyplot as plt
import math
import numpy
import string
import sys

# Utilities for reading column-oriented files
#
# H. Ferguson - revised 10/20/03 to allow use of Numeric or numarray
# H. Ferguson - revised 2/10/08 to use numpy
# H. Ferguson - revised to add some options and remove Numeric and numarray

"""Routines for reading general whitespace-delimited, column-oriented files.
   Returned values are either numpy one-dimensional arrays.
   The read routines parse the input looking for
   decimal points or non-numeric characters to decide on the format of
   the output.  Reading is therefore a bit slow, but the program interface
   is extremely simple. For example, if the file 'foo' has three columns,
   read them in as follows:
	a,b,c = fgetcols('foo')

   A few other options:
       a,b,c,d = fgetcols('foo',1,3,5,7)  # Read some selected columns
       a = fgetcols('foo')     # Read all the columns (a is then a tuple of arrays)
       a,b,c = fgetcols('foo',fs=',') # Change the field separator to a comma 
       a,b,c = fgetcols('foo',cmt='!') # Change the comment character to '!'

   The module also provides an object-oriented interface to save re-reading
   the file if multiple getcol calls are desired:
       f = readcol('foo')
       a,b = f.getcols(1,2)
       c,d = f.getcols(3,4)
       f.close()

   Ignores comment lines.
   Ignores blank lines. 
   Optionally changes INDEF to a desired value (e.g. -99.99).

   As of version 5.0, only numpy is offered (Numeric and numarray used to be 
   options).
"""

__version__ = '5.0' # Numpy is now the default
__author = 'Henry C. Ferguson, STScI'

def remove_comments(l,cmt='#'):
    comments = []
    for i in range(len(l)):
        l[i] = l[i].strip()
        if l[i] == '\n' or l[i][0] == cmt:
            comments = comments + [i]
    ngone = 0
    for i in comments:
        comment = l.pop(i-ngone)
        ngone = ngone+1
    return l

def replace_indef(l,indef):
    for i in range(len(l)):
        while string.find(l[i],"INDEF") > -1:
            idx = string.find(l[i],"INDEF")
            l[i] = l[i][:idx]+indef+l[i][idx+5:]
    return l

class readcol:
    """Column-oriented file methods."""
    def __init__(self,cfile,arraytype=numpy,indef=""):
        """Open file, read in all the lines, and return numpy arrays.
          
           Arguments:
           cfile -- file to read
           arraytype -- numpy (used to allow Numeric or numarray)
           indef -- string replacement for INDEF (e.g. NaN)
        """
        f = open(cfile,'r')
        self.l = f.readlines()
        self.l = remove_comments(self.l)
        if indef:
            self.l = replace_indef(self.l,indef)
        f.close()
        self.N = arraytype
    def getcol(self,col,fs=None):
        """Read in a single column (columns start at 1)."""
        return getcol(col,self.l,self.N,fs=fs)
    def getcols(self,*args,**kwargs):
        """Read in a multiple columns (columns start at 1)."""
        if 'fs' in list(keywords.keys()):
            fs = keywords['fs']
        else:
            fs = None
        ret = []
        for i in range(len(args)):
           ret = ret + [getcol(args[i],self.l,self.N,fs=fs)]
        return ret
    def close(self):
        """Release the memory associated with the lines read by __init__"""
        del(self.l)
        

def getcol(col,lines,N,fs=None):
  """Read in a single column from a list of strings. Parse each column to
     determine the type of variable (integer, float, string) and return 
     either an array of that type (int64, float64) or a character array.

     Arguments:
     col -- desired column (starting at 1)	
     lines -- list of strings (one per line) read from input file
     N -- numpy
  """
  i = col-1
  nlines = len(lines)
  if fs != None: # If delimiter is not whitespace, remove the whitespace
      oldlines = lines
      lines = []
      for ol in oldlines:
          lines += [string.join(ol.split())] 
  a = lines[0].split(fs) # Determine the type from the first line
#  if string.find(a[i],'.') < 0:
  if 1 < 0:
    try:
      x = int(a[i]) 
    except:
      values = list(range(nlines))
      getstrings(col,lines,values,fs=fs)
      values = N.array(values)
    else:
      values = N.zeros((nlines),N.int64)
      if type(getints(col,lines,values,fs=fs)) == type(1):
        values = N.zeros((nlines),N.float64)
        getfloats(col,lines,values,fs=fs)
  else:
    try:
      x = float(a[i]) 
    except:
      values = list(range(nlines))
      getstrings(col,lines,values,fs=fs)
      values = N.array(values)
    else:
      values = N.zeros((nlines),N.float64)
      getfloats(col,lines,values,fs=fs)
  return values

def getstrings(col,lines,values,fs=None):
  n = 0
  for l in lines:
    a = l.split(fs)
    values[n] = a[col-1]
    n = n+1

def getints(col,lines,values,fs=None):
  n = 0
  for l in lines:
    a = l.split(fs)
    if string.find(a[col-1],'.') > 0:
      return -1
    else:
      values[n] = int(a[col-1])
    n = n+1
  return values    


def getfloats(col,lines,values,fs=None):
  n = 0
  for l in lines:
    a = l.split(fs)
    values[n] = float(a[col-1])
    n = n+1


def fgetcol(cfile,col,arraytype="numpy",cmt='#',indef="-99.99"):
    """Read in a single column from a file. Parse the column to
       determine the type of variable (integer, float, string) and return 
       either an array of that type (int64, float64) or a character array.

       Arguments:
       cfile -- file to be read
       col -- desired column (starting at 1)	
       arraytype -- numpy
       indef="-99.99" (INDEF replacement string)
    """
    f = open(cfile,'r')
    l = f.readlines()
    f.close()
    l = remove_comments(l,cmt=cmt)
    if indef:
        l = replace_indef(l,indef)
    if arraytype == "numpy":
        N = numpy
    return getcol(col,l,N)

def fgetcols(cfile,*args,**keywords):
    """Read multiple columns from a file. Parse each column to
       determine the type of variable (integer, float, string) and return 
       either one-dimensional arrays of the appropriate type (int64, float64) 
       or a character array.

       Arguments:
       cfile -- file to be read
       *args -- desired columns (starting at 1)	
       **keywords -- indef="-99.99" (INDEF replacement string)
                  -- cmt="#" (comment character)
                  -- fs=None (field separator; defaults to whitespace)

       Examples:
         If the file 'foo' has three columns, read them in as follows:
	     a,b,c = fgetcols('foo')

         A few other examples:
             a,b,c,d = fgetcols('foo',1,3,5,7) # read selected columns 
             a = fgetcols('foo')               # read all columns 
             a,b,c = fgetcols('foo',fs=',')    # Change the field separator
             a,b,c = fgetcols('foo',cmt='!')   # Change the comment character to '!'

    """
    f = open(cfile,'r')
    l = f.readlines()
    f.close()
    if 'cmt' in list(keywords.keys()):
        cmt = keywords['cmt']
    else:
        cmt = '#'
    l = remove_comments(l,cmt=cmt)
    if 'indef' in list(keywords.keys()):
        indef = keywords['indef']
        l = replace_indef(l,indef)
    N = numpy
    if 'arraytype' in list(keywords.keys()):
        arraytype = keywords['arraytype']
        if arraytype != "numpy":
            print("readcol: As of v5.0, only numpy arrays are returned")
    if 'fs' in list(keywords.keys()):
        fs = keywords['fs']
    else:
        fs = None
    ret = []
    ncols = len(args)
    colnumbers = args
    if ncols == 0:       # If no columns are listed, read them all
        ncols = len(l[0].split(fs))
        colnumbers = N.array(list(range(ncols)))+1
    for i in range(ncols):
        ret = ret + [getcol(colnumbers[i],l,N,fs=fs)]
    return ret

m_base1 = '25'
m_base2 = '183'
m_base3 = '570'

##Reading data
atom = sys.argv[1]
file0 = '0-0/' + atom + '.dat'
file1 = m_base1 + '/' + atom +'.dat'
file2 = m_base2 + '/' + atom +'.dat'
file3 = m_base3 + '/' + atom +'.dat'
field1 = '25G'
field2 = '205G'
field3 = '570G'

plot_limit = 103
#Reading each column
z,h0 = fgetcols(file0)
z,h1 = fgetcols(file1)
z,h2 = fgetcols(file2)
z,h3 = fgetcols(file3)

##limiting values to significant variants
zdens = z[1:plot_limit]
h0dens = h0[1:plot_limit]
h1dens = h1[1:plot_limit]
h2dens = h2[1:plot_limit]
h3dens = h3[1:plot_limit]

##generating density plot
figure = plt.figure(num='perfil de densidades')
ax = figure.add_subplot(1, 1, 1)

plt.xlabel('z [km]')
plt.ylabel('H [cm-3]')
plt.plot(zdens, h0dens, label='0G')
plt.plot(zdens, h1dens, label = m_base1)
plt.plot(zdens, h2dens, label = m_base2)
plt.plot(zdens, h3dens, label = m_base3)
ax.set_yscale('log')
plt.legend()

##Generating absolute diff from plots
plt.figure(num = 'diferencias absolutas')
zabs = z[1:plot_limit]
abs1 = h1[1:plot_limit] - h0[1:plot_limit]
abs2 = h2[1:plot_limit] - h0[1:plot_limit]
abs3 = h3[1:plot_limit] - h0[1:plot_limit]
plt.xlabel('z [km]')
plt.ylabel('Difference in H [cm-3]')
plt.plot(zabs, abs1, label = m_base1 + ' - 0G', color='C1')
plt.plot(zabs, abs2, label = m_base2 + ' - 0G', color='g')
plt.plot(zabs, abs3, label = m_base3 + ' - 0G', color='r')
ax.set_yscale('log')
plt.legend()


##Generating relative plot
plt.figure(num = 'diferencias relativas')
zrel = z[1:plot_limit]
rel1 = h1[1:plot_limit] - h0[1:plot_limit]
rel2 = h2[1:plot_limit] - h1[1:plot_limit]
rel3 = h3[1:plot_limit] - h2[1:plot_limit]
plt.xlabel('z [km]')
plt.ylabel('Difference in H [cm-3]')
plt.plot(zrel, rel1, label = m_base1 + ' - 0G', color='C1')
plt.plot(zrel, rel2, label = m_base2 + ' - ' + m_base1, color='g')
plt.plot(zrel, rel3, label = m_base3 + ' - ' + m_base2, color='r')
ax.set_yscale('log')
plt.legend()


##Generating pressure plots
file0 = '0-0/P.dat'
file1 = m_base1 + '/P.dat'
file2 = m_base2 + '/P.dat'
file3 = m_base3 + '/P.dat'

zpres,pres = fgetcols(file0)
zpres,pres1 = fgetcols(file1)
zpres,pres2 = fgetcols(file2)
zpres,pres3 = fgetcols(file3)
zpres = zpres[1:plot_limit]
pres = pres[1:plot_limit]
pres1 = pres1[1:plot_limit]
pres2 = pres2[1:plot_limit]
pres3 = pres3[1:plot_limit]

figure = plt.figure(num='Presion')
ax = figure.add_subplot(1, 1, 1)

plt.xlabel('z [km]')
plt.ylabel('P [Pa]')
plt.plot(zpres, pres, label='0G')
plt.plot(zpres, pres1, label = m_base1)
plt.plot(zpres, pres2, label = m_base2)
plt.plot(zpres, pres3, label = m_base3)
ax.set_yscale('log')
plt.legend()


##Generating temperature plot
file0 = '0-0/T.dat'

z,T = fgetcols(file0)
z = z[1:plot_limit]
T = T[1:plot_limit]

figure = plt.figure(num='Temperatura')
ax = figure.add_subplot(1, 1, 1)
plt.xlabel('z [km]')
plt.ylabel('T [k]')
plt.plot(z, T)
ax.set_yscale('log')


##Generating magnetic plots
file0 = '0-0/magnetic.dat'
file1 = m_base1 + '/magnetic.dat'
file2 = m_base2 + '/magnetic.dat'
file3 = m_base3 + '/magnetic.dat'

zm,magnetic0 = fgetcols(file0)
zm,magnetic1 = fgetcols(file1)
zm,magnetic2 = fgetcols(file2)
zm,magnetic3 = fgetcols(file3)
zm = z
magnetic0 = magnetic0[1:plot_limit]
magnetic1 = magnetic1[1:plot_limit]
magnetic2 = magnetic2[1:plot_limit]
magnetic3 = magnetic3[1:plot_limit]

figure = plt.figure(num='Campo Magnetico')
ax = figure.add_subplot(1, 1, 1)

plt.xlabel('z [km]')
plt.ylabel('P [Pa]')
#plt.plot(zm, magnetic0, label='0G')
plt.plot(zm, magnetic1, label = m_base1, color='C1')
plt.plot(zm, magnetic2, label = m_base2, color='g')
plt.plot(zm, magnetic3, label = m_base3, color='r')
#ax.set_yscale('log')
plt.legend()



##plot everything
plt.show()
