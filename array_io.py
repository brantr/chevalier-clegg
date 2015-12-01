import numpy as np

def write_array(fname, xa):
  fp = open(fname,"w")
  for x in xa:
    s = '%lf\n' % x 
    fp.write(s)
  fp.close()

def read_array(fname):
  fp = open(fname,"r")
  fl = fp.readlines()
  x = np.fromfile(fname, dtype=np.float32, sep='\n')
  fp.close()
  return x

def write_two_arrays(fname, xa, ya):
  fp = open(fname,"w")
  for i in range(len(xa)):
    s = '% 10e\t% 10e\n' % (xa[i],ya[i])
    fp.write(s)
  fp.close()

def read_two_arrays(fname):
  fp = open(fname,"r")
  fl = fp.readlines()
  n = len(fl)
  x = np.zeros(n)
  y = np.zeros(n)
  for i in range(n):
    x[i] = float(fl[i].split()[0])
    y[i] = float(fl[i].split()[1])
  return x,y

def write_three_arrays(fname, xa, ya, za):
  fp = open(fname,"w")
  for i in range(len(xa)):
    s = '% 10e\t% 10e\t% 10e\n' % (xa[i],ya[i],za[i])
    fp.write(s)
  fp.close()

def read_three_arrays(fname):
  fp = open(fname,"r")
  fl = fp.readlines()
  n = len(fl)
  x = np.zeros(n)
  y = np.zeros(n)
  z = np.zeros(n)
  for i in range(n):
    x[i] = float(fl[i].split()[0])
    y[i] = float(fl[i].split()[1])
    z[i] = float(fl[i].split()[2])
  return x,y,z

def write_four_arrays(fname, a, b, c, d):
  fp = open(fname,"w")
  for i in range(len(a)):
    s = '% 10e\t% 10e\t% 10e\t% 10e\n' % (a[i],b[i],c[i],d[i])
    fp.write(s)
  fp.close()

def read_four_arrays(fname):
  fp = open(fname,"r")
  fl = fp.readlines()
  n = len(fl)
  a = np.zeros(n)
  b = np.zeros(n)
  c = np.zeros(n)
  d = np.zeros(n)
  for i in range(n):
    a[i] = float(fl[i].split()[0])
    b[i] = float(fl[i].split()[1])
    c[i] = float(fl[i].split()[2])
    d[i] = float(fl[i].split()[3])
  return a,b,c,d

def write_five_arrays(fname, a, b, c, d, f):
  fp = open(fname,"w")
  for i in range(len(a)):
    s = '% 10e\t% 10e\t% 10e\t% 10e\t% 10e\n' % (a[i],b[i],c[i],d[i],f[i])
    fp.write(s)
  fp.close()

def write_five_arrays_nlines(fname, a, b, c, d, f):
  fp = open(fname,"w")
  s = "%s\n" % len(a)
  fp.write(s)
  for i in range(len(a)):
    s = '% 10e\t% 10e\t% 10e\t% 10e\t% 10e\n' % (a[i],b[i],c[i],d[i],f[i])
    fp.write(s)
  fp.close()

def read_five_arrays(fname):
  fp = open(fname,"r")
  fl = fp.readlines()
  n = len(fl)
  a = np.zeros(n)
  b = np.zeros(n)
  c = np.zeros(n)
  d = np.zeros(n)
  f = np.zeros(n)
  for i in range(n):
    a[i] = float(fl[i].split()[0])
    b[i] = float(fl[i].split()[1])
    c[i] = float(fl[i].split()[2])
    d[i] = float(fl[i].split()[3])
    f[i] = float(fl[i].split()[4])
  return a,b,c,d,f

def read_eight_arrays(fname):
  fp = open(fname,"r")
  fl = fp.readlines()
  n = len(fl)
  a = np.zeros(n)
  b = np.zeros(n)
  c = np.zeros(n)
  d = np.zeros(n)
  f = np.zeros(n)
  g = np.zeros(n)
  h = np.zeros(n)
  j = np.zeros(n)
  for i in range(n):
    a[i] = float(fl[i].split()[0])
    b[i] = float(fl[i].split()[1])
    c[i] = float(fl[i].split()[2])
    d[i] = float(fl[i].split()[3])
    f[i] = float(fl[i].split()[4])
    g[i] = float(fl[i].split()[5])
    h[i] = float(fl[i].split()[6])
    j[i] = float(fl[i].split()[7])
  return a,b,c,d,f,g,h,j
