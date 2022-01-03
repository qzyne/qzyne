
from sympy import * ## KHONG XOA
import numpy as np ## KHONG XOA 
global x, y, z, t ## KHONG XOA
x, y, z, t = symbols("x, y, z, t") ## KHONG XOA

def req1(f, g, a): ## KHONG XOA 
    def check_prime_number(f1):
      check = 1
      if isinstance(f1, int) or isinstance(f1, float) != False : 
        check = 0 
      return check   
    try:
      f1a = f + g
      df1a = diff(f1a, x, 1)
      m1a = float(df1a.subs(x, a))
      gt1 = round(m1a,2)
    except:
      gt1 = None
    try: 
      f1b = f * g
      df1b = diff(f1b, x, 1)
      m1b = float(df1b.subs(x, a))
      gt2 = round(m1b,2)    
    except:
      gt2 = None 
    try:
      if check_prime_number(f) == 0:
        fcx = lambdify(x,f)
        f1c_a = fcx(g)
      else:
        f1c_a = f.subs(x,g)    
      df1c = diff(f1c_a, x, 1)
      m1c = float(df1c.subs(x, a))
      gt3 = round(m1c,2) 
    except:
      gt3 = None
    try:
      f1d = (f / g)
      df1d = diff(f1d, x, 1)
      m1d = float(df1d.subs(x, a))
      gt4 = round(m1d,2) 
    except:
      gt4 = None

    def resultf1(ham):
      co = 1
      if ham == nan or ham == zoo:
        co = 0    
      return co
    if resultf1(gt1) == 1:
      m1ax = gt1
    else:
      m1ax = None
    if resultf1(gt2) == 1:
      m1bx = gt2
    else:
      m1bx = None
    if resultf1(gt3) == 1:
      m1cx = gt3 
    else:
      m1cx = None
    if resultf1(gt4) == 1:
      m1dx = gt4
    else:
      m1dx= None
    return m1ax, m1bx, m1cx, m1dx    
    

def req2(f, a, b, c):  ## KHONG XOA
   try:
    f2 = lambdify((x,y,z), f)

    k_x = float(((diff(f, x, 1).subs(x, a)).subs(y, b)).subs(z, c))

    k_y = float(((diff(f, y, 1).subs(x, a)).subs(y, b)).subs(z, c))

    k_z = float(((diff(f, z, 1).subs(x, a)).subs(y, b)).subs(z, c))
    if f2(a, b, c) + k_x*(x - a) + k_y*(y - b) + k_z*(z - c) == nan:
      return None
    else:
      return f2(a, b, c) + k_x*(x - a) + k_y*(y - b) + k_z*(z - c)
   except:
     return None


def req3(w, f1, f2, f3, a):  ## KHONG XOA
  try:
    fw3 = ((w.subs(x, f1)).subs(y, f2)).subs(z, f3)
    df3 = diff(fw3, t, 1)
    result3 = df3.subs(t, a)
    if result3 == nan:
      return None
    else:
      return float(result3) 
  except:
    return None


def req4(a, b, n):  ## KHONG XOA
    result4 = 0.0
 
    for i in range(0, n + 1):
      f4 = (factorial(float(n))/(factorial(i)*factorial(n-i)))*(a**(n-i))*b**i
      result4 = result4 + f4
    return result4


def req5(f):  ## KHONG XOA
    f5_x = diff(f, x, 1)
    f5_y = diff(f, y, 1)
    M = solve((f5_x, f5_y), (x,y), dict = True)
    f5_xx = diff(f5_x, x, 1)
    f5_yy = diff(f5_y, y, 1)
    f5_xy = diff(f5_x, y, 1)
    D_a = f5_xx * f5_yy - (f5_xy*f5_xy)
    D = lambdify((x, y), D_a)
    GTNN = []
    GTLN = []
    yenngua = []
    rong = []

    for i in M:
      if(i[x].is_real and i[y].is_real == True):
        i[x] = float(i[x])
        i[y] = float(i[y])
        A = D(i[x], i[y])
        if A < 0:
          yenngua = i[x],i[y]
        elif A > 0 and f5_xx.subs(x, i[x]).subs(y, i[y])> 0:
          GTNN = i[x],i[y]
        elif A > 0 and f5_xx.subs(x, i[x]).subs(y, i[y]) < 0:
          GTLN = i[x],i[y]
        else:
          rong = [], [], []
          return rong
    return GTNN, GTLN, yenngua


def req6(message, x, y, z):  ## KHONG XOA
    f6 = abs(x**2 - y**2 - z)
    result6 = ''
    for i in range(0, len(message)):
      cirpher = ord(message[i])
      plain = chr(f6^cirpher) 
      result6 += plain
    return result6 


def req7(xp, yp, xc):  ## KHONG XOA
    n_xp = len(xp)
    n_yp = len(yp)    
    sum_x = 0
    sum_y = 0
    sum_xy = 0
    sum_xx = 0
    for i in range(n_xp):
      sum_x = sum_x + xp[i]
      sum_y = sum_y + yp[i]
      sum_xy = sum_xy + yp[i]*xp[i] 
      sum_xx = sum_xx + xp[i]*xp[i]
    sumxx = sum_x*sum_x
    m = (sum_x*sum_y - (n_xp*sum_xy)) / (sumxx - n_xp*sum_xx)
    b = (1/n_xp)*(sum_y - m*sum_x)
    result7 = float(m*xc + b)
    return round(result7, 2) 


def req8(f, eta, xi, tol): ## KHONG XOA
    try:
      def f8(f, value):
        df8 = diff(f, x, 1)
        df8_value = lambdify(x, df8)(value)
        return df8_value
      
      for i in range(1000):
        xi += - eta * f8(f, xi)
        if abs(f8(f, xi)) < tol:
          break
      result8 = float(xi)
      return round(result8, 2)
    except:
      return None



