#!/usr/bin/env python

"""Calculates CLs for simple one bin counting experiments.

CLs is calculated using the standard LHC techniques using a profile likelihood 
for the test statistic.
"""

__author__ = "Daniel Schmeier, Jamie Tattersall, Sebastian Belkner"
__copyright__ = "Copyright 2014, CheckMATE"
__credits__ = ["Daniel Schmeier"]
__license__ = "GPL"
__version__ = "1.0.2"
__maintainer__ = "Daniel Schmeier"
__email__ = "daschm@uni-bonn.de"
__status__ = "Prototype"


from math import *
import os, sys, scipy
import numpy as np
from math import log, exp, sqrt, factorial
import os, sys
import scipy.optimize as opt
#import cls_nocov as CN
from scipy.stats import poisson, norm
from iminuit import Minuit, describe, Struct
# Note: Scipy packages are imported within the code, using a try-except setup
  


def nexp(mu, nu_b, nu_s, b, ddb, s, dds):
  """
  Returns expected number of events given mu and nu
  expected = signal + background.
  Signal is multiplied with modifier mu.
  Signal is 'smeared' by corresponding nu_s.
  Background is 'smeared' by corresponding nu_b.
  The smearing is lognormal, which is performed by taking
    the exponential of a normal distributed random number (nu).
  Note that b0, s0, db and ds are global variables defined at input
  """
  return b*exp(ddb/b*nu_b) + mu*s*exp(dds/s*nu_s)


def reduced_logp(nuTil, nu):
  """ reduced_logp = log(p) without normalisation """
  return -(nuTil-nu)**2./2.

def func_to_max(nu_b_til, nu_s_til,nu_b_til2, nu_s_til2, nu_b_til3, nu_s_til3, nu_b_til4, nu_s_til4, nu_b_til5, nu_s_til5, nu_b_til6, nu_s_til6, mu, nu_s, nu_b,nu_s2, nu_b2, nu_s3, nu_b3, nu_s4, nu_b4, nu_s5, nu_b5, nu_s6, nu_b6):
  xx=n
  result =0

  print "n=" + str(n) 
  print "b=" + str(b0) 
  print "s=" + str(s0) 
  print "db=" + str(ds) 
  print "db=" + str(db) 
  print "n=" + str(n) 
  zwischen = [0 for i in range(len(xx))]
  lam = [0 for i in range(len(xx))]
  nu_b = [nu_b, nu_b2, nu_b3, nu_b4, nu_b5, nu_b6]
  nu_s = [nu_s, nu_s2, nu_s3, nu_s4, nu_s5, nu_s6]
  nu_b_til = [nu_b_til, nu_b_til2, nu_b_til3, nu_b_til4, nu_b_til5, nu_b_til6]
  nu_s_til = [nu_s_til, nu_s_til2, nu_s_til3, nu_s_til4, nu_s_til5, nu_s_til6]
  for i in range(len(xx)):
    lam[i]= b0[i]*exp(db[i]/b0[i]*nu_b[i]) + mu*s0[i]*exp(ds[i]/s0[i]*nu_s[i])

    try:
      zwischen[i] = -log(poisson.pmf(xx[i], lam[i] )*exp(-(nu_b_til[i]-nu_b[i])**2/2.)*exp(-(nu_s_til[i]-nu_s[i])**2/2.))
    except:
      rama = xx[i]* log(xx[i]) - xx[i] + ( log(xx[i]*(1+4*xx[i]*(1+2*xx[i])))/6.  )+ (log(pi)/2. )
      zwischen[i] = - xx[i] * log(lam[i]) + rama + lam[i] + log( exp(-(nu_b_til[i]-nu_b[i])**2/2.)*exp(-(nu_s_til[i]-nu_s[i])**2/2.) )

    
    result += zwischen[i]
  
  return result


def func_to_max_pseu_sig(nu_b_til, nu_s_til, nu_b_til2, nu_s_til2, nu_b_til3, nu_s_til3, nu_b_til4, nu_s_til4, nu_b_til5, nu_s_til5, nu_b_til6, nu_s_til6, mu, nu_s, nu_b, nu_s2, nu_b2, nu_s3, nu_b3, nu_s4, nu_b4, nu_s5, nu_b5, nu_s6, nu_b6):
  xx=n_p_sig
  result =0 
  zwischen = [0 for i in range(len(xx))]
  lam = [0 for i in range(len(xx))]
  nu_b = [nu_b, nu_b2, nu_b3, nu_b4, nu_b3, nu_b4]
  nu_s = [nu_s, nu_s2, nu_s3, nu_s4, nu_s3, nu_s4]
  nu_b_til = [nu_b_til, nu_b_til2, nu_b_til3, nu_b_til4, nu_b_til5, nu_b_til6]
  nu_s_til = [nu_s_til, nu_s_til2, nu_s_til3, nu_s_til4, nu_s_til5, nu_s_til6]
  for i in range(len(xx)):
    print "xx =" + str(xx[i])
    print "Lange xx" + str(len(xx))
    lam[i]= b0[i]*exp(db[i]/b0[i]*nu_b[i]) + mu*s0[i]*exp(ds[i]/s0[i]*nu_s[i])
    print "lam= " + str(lam[i])    
    print "poisson= " + str(poisson.pmf(xx[i], lam[i]))    
    
    try:
      zwischen[i] = -log(poisson.pmf(xx[i], lam[i] )*exp(-(nu_b_til[i]-nu_b[i])**2/2.)*exp(-(nu_s_til[i]-nu_s[i])**2/2.))
    except:
      rama = xx[i]* log(xx[i]) - xx[i] + ( log(xx[i]*(1+4*xx[i]*(1+2*xx[i])))/6.  )+ (log(pi)/2. )
      zwischen[i] = - xx[i] * log(lam[i]) + rama + lam[i] + log( exp(-(nu_b_til[i]-nu_b[i])**2/2.)*exp(-(nu_s_til[i]-nu_s[i])**2/2.) )
    result += zwischen[i]
  
  return result



def func_to_max_pseu_bkg(nu_b_til, nu_s_til, nu_b_til2, nu_s_til2, nu_b_til3, nu_s_til3, nu_b_til4, nu_s_til4, nu_b_til5, nu_s_til5, nu_b_til6, nu_s_til6, mu, nu_s, nu_b, nu_s2, nu_b2, nu_s3, nu_b3, nu_s4, nu_b4, nu_s5, nu_b5, nu_s6, nu_b6):
  xx=n_p_bkg
  result =0 
  zwischen = [0 for i in range(len(xx))]
  lam = [0 for i in range(len(xx))]
  nu_b = [nu_b, nu_b2, nu_b3, nu_b4, nu_b5, nu_b6]
  nu_s = [nu_s, nu_s2, nu_s3, nu_s4, nu_s5, nu_s6]
  nu_b_til = [nu_b_til, nu_b_til2, nu_b_til3, nu_b_til4, nu_b_til5, nu_b_til6]
  nu_s_til = [nu_s_til, nu_s_til2, nu_s_til3, nu_s_til4, nu_s_til5, nu_s_til6]
  for i in range(len(xx)):
    print "xx =" + str(xx[i])
    print "Lange xx" + str(len(xx))
    lam[i]= b0[i]*exp(db[i]/b0[i]*nu_b[i]) + mu*s0[i]*exp(ds[i]/s0[i]*nu_s[i])
    print "lam= " + str(lam[i])    
    print "poisson= " + str(poisson.pmf(xx[i], lam[i]))    
    try:
      zwischen[i] = -log(poisson.pmf(xx[i], lam[i] )*exp(-(nu_b_til[i]-nu_b[i])**2/2.)*exp(-(nu_s_til[i]-nu_s[i])**2/2.))
    except:
      rama = xx[i]* log(xx[i]) - xx[i] + ( log(xx[i]*(1+4*xx[i]*(1+2*xx[i])))/6.  )+ (log(pi)/2. )
      zwischen[i] = - xx[i] * log(lam[i]) + rama + lam[i] + log( exp(-(nu_b_til[i]-nu_b[i])**2/2.)*exp(-(nu_s_til[i]-nu_s[i])**2/2.) )
    result += zwischen[i]
  
  return result


def max_L_iminuit(nu_b_til, nu_s_til, nu_b_til2, nu_s_til2, nu_b_til3, nu_s_til3, nu_b_til4, nu_s_til4, nu_b_til5, nu_s_til5, nu_b_til6, nu_s_til6, nu_s, nu_b, nu_s2, nu_b2, nu_s3, nu_b3, nu_s4, nu_b4, nu_s5, nu_b5, nu_s6, nu_b6,  mu):
  m=Minuit(func_to_max, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True,nu_s_til2= nu_s_til2, fix_nu_s_til2=True,  nu_b_til2= nu_b_til2, fix_nu_b_til2=True,nu_s_til3= nu_s_til3, fix_nu_s_til3=True,  nu_b_til3= nu_b_til3, fix_nu_b_til3=True,nu_s_til4= nu_s_til4, fix_nu_s_til4=True,  nu_b_til4= nu_b_til4, fix_nu_b_til4=True,nu_s_til5= nu_s_til5, fix_nu_s_til5=True,  nu_b_til5= nu_b_til5, fix_nu_b_til5=True,nu_s_til6= nu_s_til6, fix_nu_s_til6=True,  nu_b_til6= nu_b_til6, fix_nu_b_til6=True, mu=0.0, fix_mu = False, limit_mu=(0,4), error_mu = 1.0,  nu_s= nu_s, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-3,3), nu_b= nu_b, fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-3,3),nu_s2= nu_s2, fix_nu_s2=False, error_nu_s2=1.0, limit_nu_s2=(-3,3), nu_b2= nu_b2, fix_nu_b2=False, error_nu_b2= 1.0, limit_nu_b2=(-3,3),nu_s3= nu_s3, fix_nu_s3=False, error_nu_s3=1.0, limit_nu_s3=(-3,3), nu_b3= nu_b3, fix_nu_b3=False, error_nu_b3= 1.0, limit_nu_b3=(-3,3),nu_s4= nu_s4, fix_nu_s4=False, error_nu_s4=1.0, limit_nu_s4=(-3,3), nu_b4= nu_b4, fix_nu_b4=False, error_nu_b4= 1.0, limit_nu_b4=(-3,3),nu_s5= nu_s5, fix_nu_s5=False, error_nu_s5=1.0, limit_nu_s5=(-3,3), nu_b5= nu_b5, fix_nu_b5=False, error_nu_b5= 1.0, limit_nu_b5=(-3,3),nu_s6= nu_s6, fix_nu_s6=False, error_nu_s6=1.0, limit_nu_s6=(-3,3), nu_b6= nu_b6, fix_nu_b6=False, error_nu_b6= 1.0, limit_nu_b6=(-3,3),print_level=0)
  m.migrad()
  #m.hesse()
  #m.minos()
  #m.migrad()
  global f_value
  f_value = m.fval
 # print 'Errors=' + str(m.errors)
#  print m.get_merrors()
 # print m.get_fmin()
 # if m.get_fmin().is_valid == False:
 #   print 'Fehler_valid'
 # if m.matrix_accurate() == False:
 #   print 'Fehler_accurate'
 # Errors= (m.get_merrors())
 # if Errors["nu_b"]["is_valid"] == False:
 #   print "Fehler_nu_b_valid"
 # if Errors["nu_s"]["is_valid"] == False:
 #   print "Fehler_nu_s_valid"
 # if Errors["mu"]["is_valid"] == False:
 #   print "Fehler_mu_valid"
#  f_value = m.args
#  print('args',m.args)
#  print('value', m.values)
#  print "Funktionswert" + str(func(*f_value))
  return m.values

def max_L_iminuit_mu1(nu_b_til, nu_s_til, nu_b_til2, nu_s_til2,nu_b_til3, nu_s_til3,nu_b_til4, nu_s_til4,nu_b_til5, nu_s_til5,nu_b_til6, nu_s_til6,nu_s, nu_b, nu_s2, nu_b2, nu_s3, nu_b3,nu_s4, nu_b4, nu_s5, nu_b5,nu_s6, nu_b6, mu):
  m=Minuit(func_to_max, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, nu_s_til2= nu_s_til2, fix_nu_s_til2=True, nu_b_til2= nu_b_til2, fix_nu_b_til2=True,nu_s_til3= nu_s_til3, fix_nu_s_til3=True, nu_b_til3= nu_b_til3, fix_nu_b_til3=True,nu_s_til4= nu_s_til4, fix_nu_s_til4=True, nu_b_til4= nu_b_til4, fix_nu_b_til4=True,nu_s_til5= nu_s_til5, fix_nu_s_til5=True, nu_b_til5= nu_b_til5, fix_nu_b_til5=True,nu_s_til6= nu_s_til6, fix_nu_s_til6=True, nu_b_til6= nu_b_til6, fix_nu_b_til6=True, mu=1.0, fix_mu = True, nu_s = nu_s_til, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-3,3), nu_b= nu_b_til, fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-3,3), nu_s2 = nu_s_til2, fix_nu_s2=False, error_nu_s2=1.0, limit_nu_s2=(-3,3), nu_b2= nu_b_til2, fix_nu_b2=False, error_nu_b2= 1.0, limit_nu_b2=(-3,3),nu_s3 = nu_s_til3, fix_nu_s3=False, error_nu_s3=1.0, limit_nu_s3=(-3,3), nu_b3= nu_b_til3, fix_nu_b3=False, error_nu_b3= 1.0, limit_nu_b3=(-3,3),nu_s4 = nu_s_til4, fix_nu_s4=False, error_nu_s4=1.0, limit_nu_s4=(-3,3), nu_b4= nu_b_til4, fix_nu_b4=False, error_nu_b4= 1.0, limit_nu_b4=(-3,3),nu_s5 = nu_s_til5, fix_nu_s5=False, error_nu_s5=1.0, limit_nu_s5=(-3,3), nu_b5= nu_b_til5, fix_nu_b5=False, error_nu_b5= 1.0, limit_nu_b5=(-3,3),nu_s6 = nu_s_til6, fix_nu_s6=False, error_nu_s6=1.0, limit_nu_s6=(-3,3), nu_b6= nu_b_til6, fix_nu_b6=False, error_nu_b6= 1.0, limit_nu_b6=(-3,3),print_level=0)
  m.migrad()
  #m.minos()
  #m.hesse()
  #m.migrad()
  #dic =m.get_fmin()

  #if dic["has_accurate_covar"]== False:
  #  print "Fehler"

  #print 'Errors= ' + str(m.errors)
#  print m.get_merrors()
  #if m.get_fmin().is_valid == False:
  #  print 'Fehler_valid'
  #if m.matrix_accurate() == False:
  #  print 'Fehler_accurate'
  global f_value_mu1
  f_value_mu1 = m.fval
 # Errors= (m.get_merrors())
 # if Errors["nu_b"]["is_valid"] == False:
 #   print "Fehler_nu_b"
 # if Errors["nu_s"]["is_valid"] == False:
 #   print "Fehler_nu_s"
#  global f_value_mu1
#  f_value_mu1 = m.fval
#  print('parameters', m.parameters)
#  print('args',m.args)
#  print('value', m.values)

  return m.values



def max_L_iminuit_mu0(nu_b_til, nu_s_til, nu_b_til2, nu_s_til2,nu_b_til3, nu_s_til3,nu_b_til4, nu_s_til4,nu_b_til5, nu_s_til5,nu_b_til6, nu_s_til6,nu_s, nu_b, nu_s2, nu_b2, nu_s3, nu_b3,nu_s4, nu_b4, nu_s5, nu_b5,nu_s6, nu_b6, mu):
  m=Minuit(func_to_max, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, nu_s_til2= nu_s_til2, fix_nu_s_til2=True, nu_b_til2= nu_b_til2, fix_nu_b_til2=True,nu_s_til3= nu_s_til3, fix_nu_s_til3=True, nu_b_til3= nu_b_til3, fix_nu_b_til3=True,nu_s_til4= nu_s_til4, fix_nu_s_til4=True, nu_b_til4= nu_b_til4, fix_nu_b_til4=True,nu_s_til5= nu_s_til5, fix_nu_s_til5=True, nu_b_til5= nu_b_til5, fix_nu_b_til5=True,nu_s_til6= nu_s_til6, fix_nu_s_til6=True, nu_b_til6= nu_b_til6, fix_nu_b_til6=True, mu=0.0, fix_mu = True, nu_s = nu_s_til, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-3,3), nu_b= nu_b_til, fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-3,3), nu_s2 = nu_s_til2, fix_nu_s2=False, error_nu_s2=1.0, limit_nu_s2=(-3,3), nu_b2= nu_b_til2, fix_nu_b2=False, error_nu_b2= 1.0, limit_nu_b2=(-3,3),nu_s3 = nu_s_til3, fix_nu_s3=False, error_nu_s3=1.0, limit_nu_s3=(-3,3), nu_b3= nu_b_til3, fix_nu_b3=False, error_nu_b3= 1.0, limit_nu_b3=(-3,3),nu_s4 = nu_s_til4, fix_nu_s4=False, error_nu_s4=1.0, limit_nu_s4=(-3,3), nu_b4= nu_b_til4, fix_nu_b4=False, error_nu_b4= 1.0, limit_nu_b4=(-3,3),nu_s5 = nu_s_til5, fix_nu_s5=False, error_nu_s5=1.0, limit_nu_s5=(-3,3), nu_b5= nu_b_til5, fix_nu_b5=False, error_nu_b5= 1.0, limit_nu_b5=(-3,3),nu_s6 = nu_s_til6, fix_nu_s6=False, error_nu_s6=1.0, limit_nu_s6=(-3,3), nu_b6= nu_b_til6, fix_nu_b6=False, error_nu_b6= 1.0, limit_nu_b6=(-3,3),print_level=0)
  m.migrad()
 # m.minos()
 # m.hesse()
 # m.migrad()
 # dic =m.get_fmin()

 # if dic["has_accurate_covar"]== False:
 #   print "Fehler"

 # print 'Errors= ' + str(m.errors)
#  print m.get_merrors()
 # if m.get_fmin().is_valid == False:
 #   print 'Fehler_valid'
 # if m.matrix_accurate() == False:
 #   print 'Fehler_accurate'
  global f_value_mu0
  f_value_mu0 = m.fval
 # Errors= (m.get_merrors())
 # if Errors["nu_b"]["is_valid"] == False:
 #   print "Fehler_nu_b"
 # if Errors["nu_s"]["is_valid"] == False:
 #   print "Fehler_nu_s"
#  global f_value_mu1
#  f_value_mu1 = m.fval
#  print('parameters', m.parameters)
#  print('args',m.args)
#  print('value', m.values)

  return m.values


def max_L_iminuit_pseu_sig(nu_b_til, nu_s_til, nu_b_til2, nu_s_til2,nu_b_til3, nu_s_til3,nu_b_til4, nu_s_til4,nu_b_til5, nu_s_til5,nu_b_til6, nu_s_til6,nu_s, nu_b, nu_s2, nu_b2,nu_s3, nu_b3,nu_s4, nu_b4,nu_s5, nu_b5,nu_s6, nu_b6,mu):
  m=Minuit(func_to_max_pseu_sig, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, nu_s_til2= nu_s_til2, fix_nu_s_til2=True,  nu_b_til2= nu_b_til2, fix_nu_b_til2=True, nu_s_til3= nu_s_til3, fix_nu_s_til3=True,  nu_b_til3= nu_b_til3, fix_nu_b_til3=True,nu_s_til4= nu_s_til4, fix_nu_s_til4=True,  nu_b_til4= nu_b_til4, fix_nu_b_til4=True, nu_s_til5= nu_s_til5, fix_nu_s_til5=True,  nu_b_til5= nu_b_til5, fix_nu_b_til5=True,nu_s_til6= nu_s_til6, fix_nu_s_til6=True,  nu_b_til6= nu_b_til6, fix_nu_b_til6=True,mu=1.0, fix_mu = False, limit_mu=(0.0,4), error_mu= 1.0,nu_s= nu_s, fix_nu_s=False, error_nu_s=0.1, limit_nu_s=(-3,3), nu_b= nu_b, fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-3,3), nu_s2= nu_s2, fix_nu_s2=False, error_nu_s2=0.1, limit_nu_s2=(-3,3), nu_b2= nu_b2, fix_nu_b2=False, error_nu_b2= 1.0, limit_nu_b2=(-3,3),nu_s3= nu_s3, fix_nu_s3=False, error_nu_s3=0.1, limit_nu_s3=(-3,3), nu_b3= nu_b3, fix_nu_b3=False, error_nu_b3= 1.0, limit_nu_b3=(-3,3),nu_s4= nu_s4, fix_nu_s4=False, error_nu_s4=0.1, limit_nu_s4=(-3,3), nu_b4= nu_b4, fix_nu_b4=False, error_nu_b4= 1.0, limit_nu_b4=(-3,3),nu_s5= nu_s5, fix_nu_s5=False, error_nu_s5=0.1, limit_nu_s5=(-3,3), nu_b5= nu_b5, fix_nu_b5=False, error_nu_b5= 1.0, limit_nu_b5=(-3,3),nu_s6= nu_s6, fix_nu_s6=False, error_nu_s6=0.1, limit_nu_s6=(-3,3), nu_b6= nu_b6, fix_nu_b6=False, error_nu_b6= 1.0, limit_nu_b6=(-3,3),print_level=1)
  m.migrad(ncall = 50000)
  #m.minos()
  #m.hesse()
#  m.migrad()
#  dic =m.get_fmin()

#  if dic["has_accurate_covar"]== False:
#    print "Fehler_accurate"

#  print 'Errors= ' + str(m.errors)
#  print m.get_merrors()
#  if m.get_fmin().is_valid == False:
#    print 'Fehler_valid'
#  if m.matrix_accurate() == False:
#    print 'Fehler'
 # global f_value_mu1
 # f_value_mu1 = m.fval
 # Errors= (m.get_merrors())
 # if Errors["nu_b"]["is_valid"] == False:
  #  print "Fehler_nu_b"
#  if Errors["nu_s"]["is_valid"] == False:
#    print "Fehler_nu_s"
 # if Errors["mu"]["is_valid"] == False:
 #   print "Fehler_mu"
  global f_value_pseu_sig
  f_value_pseu_sig = m.fval
#  print('args',m.args)
#  print('value', m.values)
#  print "Funktionswert" + str(func(*f_value))
  return m.values

def max_L_iminuit_mu1_pseu_sig(nu_b_til, nu_s_til, nu_b_til2, nu_s_til2,nu_b_til3, nu_s_til3,nu_b_til4, nu_s_til4,nu_b_til5, nu_s_til5,nu_b_til6, nu_s_til6,nu_s, nu_b, nu_s2, nu_b2,nu_s3, nu_b3,nu_s4, nu_b4,nu_s5, nu_b5,nu_s6, nu_b6,mu):
  m=Minuit(func_to_max_pseu_sig, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, nu_s_til2= nu_s_til2, fix_nu_s_til2=True,  nu_b_til2= nu_b_til2, fix_nu_b_til2=True,nu_s_til3= nu_s_til3, fix_nu_s_til3=True,  nu_b_til3= nu_b_til3, fix_nu_b_til3=True,nu_s_til4= nu_s_til4, fix_nu_s_til4=True,  nu_b_til4= nu_b_til4, fix_nu_b_til4=True,nu_s_til5= nu_s_til5, fix_nu_s_til5=True,  nu_b_til5= nu_b_til5, fix_nu_b_til5=True,nu_s_til6= nu_s_til6, fix_nu_s_til6=True,  nu_b_til6= nu_b_til6, fix_nu_b_til6=True,mu=1.0, fix_mu = True, nu_s= nu_s_til, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-4,4), nu_b= nu_b_til , fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-4,4), nu_s2= nu_s_til2, fix_nu_s2=False, error_nu_s2=1.0, limit_nu_s2=(-4,4), nu_b2= nu_b_til2 , fix_nu_b2=False, error_nu_b2= 1.0, limit_nu_b2=(-4,4),nu_s3= nu_s_til3, fix_nu_s3=False, error_nu_s3=1.0, limit_nu_s3=(-4,4), nu_b3= nu_b_til3 , fix_nu_b3=False, error_nu_b3= 1.0, limit_nu_b3=(-4,4),nu_s4= nu_s_til4, fix_nu_s4=False, error_nu_s4=1.0, limit_nu_s4=(-4,4), nu_b4= nu_b_til4 , fix_nu_b4=False, error_nu_b4= 1.0, limit_nu_b4=(-4,4),nu_s5= nu_s_til5, fix_nu_s5=False, error_nu_s5=1.0, limit_nu_s5=(-4,4), nu_b5= nu_b_til5 , fix_nu_b5=False, error_nu_b5= 1.0, limit_nu_b5=(-4,4),nu_s6= nu_s_til6, fix_nu_s6=False, error_nu_s6=1.0, limit_nu_s6=(-4,4), nu_b6= nu_b_til6 , fix_nu_b6=False, error_nu_b6= 1.0, limit_nu_b6=(-4,4),print_level=1)
  m.migrad(ncall = 50000)
  #m.minos()
#  m.hesse()
#  m.migrad()
#  print 'Errors=' + str(m.errors)
#  print m.get_merrors()
  print m.get_fmin()
#  if m.get_fmin().is_valid == False:
 #   print 'Fehler_valid'
#  if m.matrix_accurate() == False:
#    print 'Fehler_accurate'
#  Errors= (m.get_merrors())
#  if Errors["nu_b"]["is_valid"] == False:
#    print "Fehler_nu_b"
#  if Errors["nu_s"]["is_valid"] == False:
#    print "Fehler_nu_s"
  global f_value_mu1_pseu_sig
  f_value_mu1_pseu_sig = m.fval

  return m.values



def max_L_iminuit_pseu_bkg(nu_b_til, nu_s_til, nu_b_til2, nu_s_til2,nu_b_til3, nu_s_til3,nu_b_til4, nu_s_til4,nu_b_til5, nu_s_til5,nu_b_til6, nu_s_til6,nu_s, nu_b, nu_s2, nu_b2,nu_s3, nu_b3,nu_s4, nu_b4,nu_s5, nu_b5,nu_s6, nu_b6,mu):
  m=Minuit(func_to_max_pseu_bkg, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, nu_s_til2= nu_s_til2, fix_nu_s_til2=True,  nu_b_til2= nu_b_til2, fix_nu_b_til2=True, nu_s_til3= nu_s_til3, fix_nu_s_til3=True,  nu_b_til3= nu_b_til3, fix_nu_b_til3=True,nu_s_til4= nu_s_til4, fix_nu_s_til4=True,  nu_b_til4= nu_b_til4, fix_nu_b_til4=True, nu_s_til5= nu_s_til5, fix_nu_s_til5=True,  nu_b_til5= nu_b_til5, fix_nu_b_til5=True,nu_s_til6= nu_s_til6, fix_nu_s_til6=True,  nu_b_til6= nu_b_til6, fix_nu_b_til6=True,mu=1.0, fix_mu = False, limit_mu=(0.0,4), error_mu= 1.0, nu_s= nu_s, fix_nu_s=False, error_nu_s=0.1, limit_nu_s=(-3,3), nu_b= nu_b, fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-3,3), nu_s2= nu_s2, fix_nu_s2=False, error_nu_s2=0.1, limit_nu_s2=(-3,3), nu_b2= nu_b2, fix_nu_b2=False, error_nu_b2= 1.0, limit_nu_b2=(-3,3),nu_s3= nu_s3, fix_nu_s3=False, error_nu_s3=0.1, limit_nu_s3=(-3,3), nu_b3= nu_b3, fix_nu_b3=False, error_nu_b3= 1.0, limit_nu_b3=(-3,3),nu_s4= nu_s4, fix_nu_s4=False, error_nu_s4=0.1, limit_nu_s4=(-3,3), nu_b4= nu_b4, fix_nu_b4=False, error_nu_b4= 1.0, limit_nu_b4=(-3,3),nu_s5= nu_s5, fix_nu_s5=False, error_nu_s5=0.1, limit_nu_s5=(-3,3), nu_b5= nu_b5, fix_nu_b5=False, error_nu_b5= 1.0, limit_nu_b5=(-3,3),nu_s6= nu_s6, fix_nu_s6=False, error_nu_s6=0.1, limit_nu_s6=(-3,3), nu_b6= nu_b6, fix_nu_b6=False, error_nu_b6= 1.0, limit_nu_b6=(-3,3),print_level=1)
  m.migrad()
  #m.minos()
  #m.hesse()
  #m.migrad()
  #dic =m.get_fmin()

  #if dic["has_accurate_covar"]== False:
  #  print "Fehler_accurate"

  #print 'Errors= ' + str(m.errors)
#  print m.get_merrors()
  #if m.get_fmin().is_valid == False:
  #  print 'Fehler_valid'
  #if m.matrix_accurate() == False:
  #  print 'Fehler'
  #Errors= (m.get_merrors())
  #if Errors["nu_b"]["is_valid"] == False:
  #  print "Fehler_nu_b"
  #if Errors["nu_s"]["is_valid"] == False:
  #  print "Fehler_nu_s"
  #if Errors["mu"]["is_valid"] == False:
  #  print "Fehler_mu"
  global f_value_pseu_bkg
  f_value_pseu_bkg = m.fval
#  print('args',m.args)
#  print('value', m.values)
#  print "Funktionswert" + str(func(*f_value))
  return m.values



def max_L_iminuit_mu1_pseu_bkg(nu_b_til, nu_s_til, nu_b_til2, nu_s_til2,nu_b_til3, nu_s_til3,nu_b_til4, nu_s_til4,nu_b_til5, nu_s_til5,nu_b_til6, nu_s_til6,nu_s, nu_b, nu_s2, nu_b2,nu_s3, nu_b3,nu_s4, nu_b4,nu_s5, nu_b5,nu_s6, nu_b6,mu):
  m=Minuit(func_to_max_pseu_bkg, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, nu_s_til2= nu_s_til2, fix_nu_s_til2=True,  nu_b_til2= nu_b_til2, fix_nu_b_til2=True,nu_s_til3= nu_s_til3, fix_nu_s_til3=True,  nu_b_til3= nu_b_til3, fix_nu_b_til3=True,nu_s_til4= nu_s_til4, fix_nu_s_til4=True,  nu_b_til4= nu_b_til4, fix_nu_b_til4=True,nu_s_til5= nu_s_til5, fix_nu_s_til5=True,  nu_b_til5= nu_b_til5, fix_nu_b_til5=True,nu_s_til6= nu_s_til6, fix_nu_s_til6=True,  nu_b_til6= nu_b_til6, fix_nu_b_til6=True,mu=1.0, fix_mu = True, nu_s= nu_s_til, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-3,3), nu_b= nu_b_til , fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-3,3), nu_s2= nu_s_til2, fix_nu_s2=False, error_nu_s2=1.0, limit_nu_s2=(-3,3), nu_b2= nu_b_til2 , fix_nu_b2=False, error_nu_b2= 1.0, limit_nu_b2=(-3,3),nu_s3= nu_s_til3, fix_nu_s3=False, error_nu_s3=1.0, limit_nu_s3=(-3,3), nu_b3= nu_b_til3 , fix_nu_b3=False, error_nu_b3= 1.0, limit_nu_b3=(-3,3),nu_s4= nu_s_til4, fix_nu_s4=False, error_nu_s4=1.0, limit_nu_s4=(-3,3), nu_b4= nu_b_til4 , fix_nu_b4=False, error_nu_b4= 1.0, limit_nu_b4=(-3,3),nu_s5= nu_s_til5, fix_nu_s5=False, error_nu_s5=1.0, limit_nu_s5=(-3,3), nu_b5= nu_b_til5 , fix_nu_b5=False, error_nu_b5= 1.0, limit_nu_b5=(-3,3),nu_s6= nu_s_til6, fix_nu_s6=False, error_nu_s6=1.0, limit_nu_s6=(-3,3), nu_b6= nu_b_til6 , fix_nu_b6=False, error_nu_b6= 1.0, limit_nu_b6=(-3,3),print_level=1)
  m.migrad()
  #m.hesse()
  #m.minos()
  #m.migrad()
  #print 'Errors=' + str(m.errors)
#  print m.get_merrors()
 # print m.get_fmin()
 # if m.get_fmin().is_valid == False:
 #   print 'Fehler_valid'
 # if m.matrix_accurate() == False:
 #   print 'Fehler_accurate'
 # Errors= (m.get_merrors())
 # if Errors["nu_b"]["is_valid"] == False:
 #   print "Fehler_nu_b"
 # if Errors["nu_s"]["is_valid"] == False:
  #  print "Fehler_nu_s"
  global f_value_mu1_pseu_bkg
  f_value_mu1_pseu_bkg = m.fval

  return m.values


def qMu(n, nuTil_b, nuTil_b2, nuTil_b3, nuTil_b4, nuTil_s, nuTil_s2, nuTil_s3,nuTil_s4, nuTil_s5,nuTil_s6, b, db , s, ds, xx,nr):  
  """ Test statisctic for signal+background hypothesis """
  if xx =='obs':
   
    L_max_mu1 = max_L_iminuit_mu1(nuTil_b,nuTil_s,nuTil_b2, nuTil_s2, nuTil_b3, nuTil_s3, nuTil_b4, nuTil_s4, nuTil_b5, nuTil_s5, nuTil_b6, nuTil_s6,0,0,0,0,0,0,0,0,0,0,0,0,1)
    nu_b_mu_hat, nu_s_mu_hat = L_max_mu1['nu_b'], L_max_mu1['nu_s']
    
    nu_b_mu_hat2, nu_s_mu_hat2 = L_max_mu1['nu_b2'], L_max_mu1['nu_s2']
    nu_b_mu_hat3, nu_s_mu_hat3 = L_max_mu1['nu_b3'], L_max_mu1['nu_s3']
    nu_b_mu_hat4, nu_s_mu_hat4 = L_max_mu1['nu_b4'], L_max_mu1['nu_s4']
    nu_b_mu_hat5, nu_s_mu_hat5 = L_max_mu1['nu_b5'], L_max_mu1['nu_s5']
    nu_b_mu_hat6, nu_s_mu_hat6 = L_max_mu1['nu_b6'], L_max_mu1['nu_s6']
#    print "nu_b_mu_hat, " + str(nu_b_mu_hat) + "  "+ str( nu_s_mu_hat)
    L_max = max_L_iminuit(nuTil_b,nuTil_s, nuTil_b, nuTil_s, nuTil_b, nuTil_s, nuTil_b, nuTil_s, nuTil_b, nuTil_s, nuTil_b, nuTil_s, nu_s_mu_hat, nu_b_mu_hat, nu_s_mu_hat2, nu_b_mu_hat2,nu_s_mu_hat3, nu_b_mu_hat3,nu_s_mu_hat4, nu_b_mu_hat4,nu_s_mu_hat5, nu_b_mu_hat5,nu_s_mu_hat6, nu_b_mu_hat6, 1.0)
    mu_hat = L_max['mu']
#    print "mu_hat_obs=" + str(mu_hat)
    nu_b_hat = L_max['nu_b'] 
    nu_s_hat = L_max['nu_s']
  #print nu_b_mu_hat
    nu_b_hat2 = L_max['nu_b2'] 
    nu_s_hat2 = L_max['nu_s2']
  #print nu_s_mu_hat  
  # Having this, formulate the profile likelihood  
    nu_b_hat3 = L_max['nu_b3'] 
    nu_s_hat3 = L_max['nu_s3']
    nu_b_hat4 = L_max['nu_b4'] 
    nu_s_hat4 = L_max['nu_s4']
    nu_b_hat5 = L_max['nu_b5'] 
    nu_s_hat5 = L_max['nu_s5']
    nu_b_hat6 = L_max['nu_b6'] 
    nu_s_hat6 = L_max['nu_s6']
  # First, evaluate lambda for numerator and denominator
    if nr ==0:
      lam_hat_hat = nexp(mu_hat, nu_b_hat, nu_s_hat, b, db, s, ds)
      lam_hat = nexp(1., nu_b_mu_hat, nu_s_mu_hat, b, db, s,ds)
    elif nr ==1:
  
      lam_hat_hat = nexp(mu_hat, nu_b_hat2, nu_s_hat2, b, db, s, ds)
      lam_hat = nexp(1., nu_b_mu_hat2, nu_s_mu_hat2, b, db, s,ds)
    elif nr ==2:
  
      lam_hat_hat = nexp(mu_hat, nu_b_hat3, nu_s_hat3, b, db, s, ds)
      lam_hat = nexp(1., nu_b_mu_hat3, nu_s_mu_hat3, b, db, s,ds)
    elif nr ==3:
  
      lam_hat_hat = nexp(mu_hat, nu_b_hat4, nu_s_hat4, b, db, s, ds)
      lam_hat = nexp(1., nu_b_mu_hat4, nu_s_mu_hat4, b, db, s,ds)
    elif nr ==4:
  
      lam_hat_hat = nexp(mu_hat, nu_b_hat5, nu_s_hat5, b, db, s, ds)
      lam_hat = nexp(1., nu_b_mu_hat5, nu_s_mu_hat5, b, db, s,ds)
    elif nr ==5:
  
      lam_hat_hat = nexp(mu_hat, nu_b_hat6, nu_s_hat6, b, db, s, ds)
      lam_hat = nexp(1., nu_b_mu_hat6, nu_s_mu_hat6, b, db, s,ds)
    else:
      print "Problem, keine SR"
    #lam_hat_hat = nexp(mu_hat, nu_b_hat, nu_s_hat, b, db, s, ds)
    #lam_hat = nexp(1., nu_b_mu_hat, nu_s_mu_hat, b, db, s,ds)

  if xx =='pseu_sig':
    
  # Now, Maximise with mu(here: 1) fixed
    L_max_mu1 = max_L_iminuit_mu1_pseu_sig(nuTil_b,nuTil_s,nuTil_b2, nuTil_s2, nuTil_b3, nuTil_s3, nuTil_b4, nuTil_b5, nuTil_s5, nuTil_b6, nuTil_s6,0,0,0,0,0,0,0,0,0,0,0,1)
    nu_b_mu_hat, nu_s_mu_hat = L_max_mu1['nu_b'], L_max_mu1['nu_s']
#    print "nu_b_mu_hat" + str(nu_b_mu_hat) + str( nu_s_mu_hat)
    nu_b_mu_hat2, nu_s_mu_hat2 = L_max_mu1['nu_b2'], L_max_mu1['nu_s2']
    nu_b_mu_hat3, nu_s_mu_hat3 = L_max_mu1['nu_b3'], L_max_mu1['nu_s3']
    nu_b_mu_hat4, nu_s_mu_hat4 = L_max_mu1['nu_b4'], L_max_mu1['nu_s4']
    nu_b_mu_hat5, nu_s_mu_hat5 = L_max_mu1['nu_b5'], L_max_mu1['nu_s5']
    nu_b_mu_hat6, nu_s_mu_hat6 = L_max_mu1['nu_b6'], L_max_mu1['nu_s6']
    L_max = max_L_iminuit_pseu_sig(nuTil_b,nuTil_s,nuTil_b2,nuTil_s2,nuTil_b3,nuTil_s3,nuTil_b4,nuTil_s4,nuTil_b5,nuTil_s5,nuTil_b6,nuTil_s6,nu_s_mu_hat,nu_b_mu_hat,nu_s_mu_hat2,nu_b_mu_hat2,nu_s_mu_hat3,nu_b_mu_hat3,nu_s_mu_hat4,nu_s_mu_hat5,nu_b_mu_hat5,nu_s_mu_hat6,nu_b_mu_hat6,1)
    mu_hat = L_max['mu']
#    print "mu_hat="+ str(mu_hat)
    nu_b_hat = L_max['nu_b'] 
    nu_s_hat = L_max['nu_s']
  #print nu_b_mu_hat, nu_s_mu_hat
    nu_b_hat2 = L_max['nu_b2'] 
    nu_s_hat2 = L_max['nu_s2']
  #print nu_b_mu_hat
    nu_b_hat3 = L_max['nu_b3'] 
    nu_s_hat3 = L_max['nu_s3']
    nu_b_hat4 = L_max['nu_b4'] 
    nu_s_hat4 = L_max['nu_s4']
    nu_b_hat5 = L_max['nu_b5'] 
    nu_s_hat5 = L_max['nu_s5']
    nu_b_hat6 = L_max['nu_b6'] 
    nu_s_hat6 = L_max['nu_s6']
  #print nu_s_mu_hat  
  # Having this, formulate the profile likelihood  
  # First, evaluate lambda for numerator and denominator
    lam_hat_hat = 0
    lam_hat = 0
    if nr ==0:
      lam_hat_hat = nexp(mu_hat, nu_b_hat, nu_s_hat, b, db, s, ds)
      lam_hat = nexp(1., nu_b_mu_hat, nu_s_mu_hat, b, db, s,ds)
    elif nr ==1:
  
      lam_hat_hat = nexp(mu_hat, nu_b_hat2, nu_s_hat2, b, db, s, ds)
      lam_hat = nexp(1., nu_b_mu_hat2, nu_s_mu_hat2, b, db, s,ds)
    elif nr ==2:
  
      lam_hat_hat = nexp(mu_hat, nu_b_hat3, nu_s_hat3, b, db, s, ds)
      lam_hat = nexp(1., nu_b_mu_hat3, nu_s_mu_hat3, b, db, s,ds)
    elif nr ==3:
  
      lam_hat_hat = nexp(mu_hat, nu_b_hat4, nu_s_hat4, b, db, s, ds)
      lam_hat = nexp(1., nu_b_mu_hat4, nu_s_mu_hat4, b, db, s,ds)
    elif nr ==4:
  
      lam_hat_hat = nexp(mu_hat, nu_b_hat5, nu_s_hat5, b, db, s, ds)
      lam_hat = nexp(1., nu_b_mu_hat5, nu_s_mu_hat5, b, db, s,ds)
    elif nr ==5:
  
      lam_hat_hat = nexp(mu_hat, nu_b_hat6, nu_s_hat6, b, db, s, ds)
      lam_hat = nexp(1., nu_b_mu_hat6, nu_s_mu_hat6, b, db, s,ds)
    else:
      print "Problem, keine SR"
#    print "lam_hat_hat=" + str(lam_hat_hat)
#    print "lam_hat=" + str(lam_hat)
  
  # If any lambda is almost zero, it should be set to zero
  
#  if (lam_hat < 0) and (fabs(lam_hat) < 1E-5):
#    lam_hat = 0
#  if (lam_hat_hat < 0) and (fabs(lam_hat_hat) < 1E-5):
#    lam_hat_hat = 0
    
  # If any lambda is really negative, there has been a problem somewhere
#  if (lam_hat < 0) or (lam_hat_hat < 0):
#    exit("ERROR: In S95 calculation, profile likelihood asked for nexp < 0!")
    
  result1 = 0
  # If any lambda is zero,...
  if lam_hat == 0 or lam_hat_hat == 0:
    # ... both lambdas must be zero (or so small that their product is 0)
    if lam_hat != 0 or lam_hat_hat != 0:
      exit("ERROR: In S95 calculation, profile likelihood asked for impossible maximisation parameters! ("+str((b0, s0, db, ds))+")")
    # In that case, the result is only the likelihood w.r.t the nuisance parameters
    result = -2.*(reduced_logp(nuTil_b, nu_b_mu_hat) + reduced_logp(nuTil_s, nu_s_mu_hat) - reduced_logp(nuTil_b, nu_b_hat) - reduced_logp(nuTil_s, nu_s_hat))
  else:
    if nr ==0:
    # Otherwise, return the normal likelihood    
      result1 = -2.*(n*log(lam_hat/lam_hat_hat) - (lam_hat-lam_hat_hat) + reduced_logp(nuTil_b, nu_b_mu_hat) + reduced_logp(nuTil_s, nu_s_mu_hat) - reduced_logp(nuTil_b, nu_b_hat) - reduced_logp(nuTil_s, nu_s_hat))
    elif nr ==1:
    # Otherwise, return the normal likelihood    
      result1 = -2.*(n*log(lam_hat/lam_hat_hat) - (lam_hat-lam_hat_hat) + reduced_logp(nuTil_b2, nu_b_mu_hat2) + reduced_logp(nuTil_s2, nu_s_mu_hat2) - reduced_logp(nuTil_b2, nu_b_hat2) - reduced_logp(nuTil_s2, nu_s_hat2))
    elif nr ==2:
    # Otherwise, return the normal likelihood    
      result1 = -2.*(n*log(lam_hat/lam_hat_hat) - (lam_hat-lam_hat_hat) + reduced_logp(nuTil_b3, nu_b_mu_hat3) + reduced_logp(nuTil_s3, nu_s_mu_hat3) - reduced_logp(nuTil_b3, nu_b_hat3) - reduced_logp(nuTil_s3, nu_s_hat3))
    elif nr ==3:
    # Otherwise, return the normal likelihood    
      result1 = -2.*(n*log(lam_hat/lam_hat_hat) - (lam_hat-lam_hat_hat) + reduced_logp(nuTil_b4, nu_b_mu_hat4) + reduced_logp(nuTil_s4, nu_s_mu_hat4) - reduced_logp(nuTil_b4, nu_b_hat4) - reduced_logp(nuTil_s4, nu_s_hat4))
    elif nr ==4:
    # Otherwise, return the normal likelihood    
      result1 = -2.*(n*log(lam_hat/lam_hat_hat) - (lam_hat-lam_hat_hat) + reduced_logp(nuTil_b5, nu_b_mu_hat5) + reduced_logp(nuTil_s5, nu_s_mu_hat5) - reduced_logp(nuTil_b5, nu_b_hat5) - reduced_logp(nuTil_s5, nu_s_hat5))
    elif nr ==5:
    # Otherwise, return the normal likelihood    
      result1 = -2.*(n*log(lam_hat/lam_hat_hat) - (lam_hat-lam_hat_hat) + reduced_logp(nuTil_b6, nu_b_mu_hat6) + reduced_logp(nuTil_s6, nu_s_mu_hat6) - reduced_logp(nuTil_b6, nu_b_hat6) - reduced_logp(nuTil_s6, nu_s_hat6))
    else:
      print "Problem, keine SR"
  if result1 < 0 and result1 > -0.002:
    result1 = 0
  # We have to return a rounded result as E-15 differences in the floats can
  # have weird effects on the CLs+b and CLb results for special cases
  return round(result1, 10)

def calc_CLs(n_obs, nuTil_b_obs, b0_in, db_in, nuTil_s_obs, s0_in, ds_in, nPseudo, nr_SR, bla1, bla2, bla3):
#def calc_CLs( b0_in, db_in, s0_in, ds_in):
  """ Calcuates the CLs value given the measurement and signal/background estimates by running nPseudo experiments """

  nr = nr_SR

  NSIG = nPseudo*5  # You need more signal MC as only few will pass
  global b0, s0, db, ds, n # Passes this to nexp function
  
  b0 = b0_in
  s0 = s0_in
  db = db_in
  ds = ds_in
  n = n_obs
  print "Test"
  print "n=" + str(n)
  NBKG = nPseudo
  qMu_obs=[0 for a in range(nr)]
#  nu_b_hat_obs_mu = [0 for a in range(n)]
#  nu_s_hat_obs_mu = [0 for a in range(n)]
# nu_b_hat_obs_0 = [0 for a in range(n)]
# nu_s_hat_obs_0 = [0 for a in range(n)]

  nexp_sig = [0 for a in range (nr)]
  n_pseu_sig = [0 for a in range(nr)]
  nuTil_b_pseu_sig = [0 for a in range(6)]
  nuTil_s_pseu_sig = [0 for a in range(6)]

  nexp_bkg = [0 for a in range(nr)]
  n_pseu_bkg = [0 for a in range(nr)]
  nuTil_b_pseu_bkg = [0 for a in range(6)]
  nuTil_s_pseu_bkg = [0 for a in range(6)]


  print "n=" + str(n) 
  L_max_mu1 = max_L_iminuit_mu1(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0,0,0,0,0,0,0,0,0,0,0,0,1.)
  L1_obs = f_value_mu1
  print "L1_obs= " + str(L1_obs)
  L_max = max_L_iminuit(0,0,0,0,0,0,0,0,0,0,0,0,L_max_mu1["nu_s"],L_max_mu1["nu_b"],L_max_mu1["nu_s2"],L_max_mu1["nu_b2"],L_max_mu1["nu_s3"],L_max_mu1["nu_b3"],L_max_mu1["nu_s4"],L_max_mu1["nu_b4"],L_max_mu1["nu_s5"],L_max_mu1["nu_b5"],L_max_mu1["nu_s6"],L_max_mu1["nu_b6"],1.)
  L2_obs = f_value
  print "L2_obs= " + str(L2_obs)

  if L_max["mu"] > 1:
    qMu_tot_obs = 0
  else:
    qMu_tot_obs = (-2*(L2_obs-L1_obs))


  #qMu_tot_obs = (-2*(L2_obs-L1_obs))


  nu_b_hat_obs_mu, nu_s_hat_obs_mu = L_max_mu1['nu_b'], L_max_mu1['nu_s']
  nu_b_hat_obs_mu2, nu_s_hat_obs_mu2 = L_max_mu1['nu_b2'], L_max_mu1['nu_s2']
  nu_b_hat_obs_mu3, nu_s_hat_obs_mu3 = L_max_mu1['nu_b3'], L_max_mu1['nu_s3']
  nu_b_hat_obs_mu4, nu_s_hat_obs_mu4 = L_max_mu1['nu_b4'], L_max_mu1['nu_s4']
  nu_b_hat_obs_mu5, nu_s_hat_obs_mu5 = L_max_mu1['nu_b5'], L_max_mu1['nu_s5']
  nu_b_hat_obs_mu6, nu_s_hat_obs_mu6 = L_max_mu1['nu_b6'], L_max_mu1['nu_s6']


  nu_b_hat_obs_mu_v = [nu_b_hat_obs_mu, nu_b_hat_obs_mu2, nu_b_hat_obs_mu3, nu_b_hat_obs_mu4, nu_b_hat_obs_mu5, nu_b_hat_obs_mu6]
  nu_s_hat_obs_mu_v = [nu_s_hat_obs_mu, nu_s_hat_obs_mu2, nu_s_hat_obs_mu3, nu_s_hat_obs_mu4, nu_s_hat_obs_mu5, nu_s_hat_obs_mu6]


  print "nu_b_hat_obs" + str(nu_b_hat_obs_mu)
  nuTil_b_pseu_sig[0] = np.random.normal(nu_b_hat_obs_mu, 1., NSIG)  
  nuTil_s_pseu_sig[0] = np.random.normal(nu_s_hat_obs_mu, 1., NSIG)  
  nuTil_b_pseu_sig[1] = np.random.normal(nu_b_hat_obs_mu2, 1., NSIG)  
  nuTil_s_pseu_sig[1] = np.random.normal(nu_s_hat_obs_mu2, 1., NSIG)  
  nuTil_b_pseu_sig[2] = np.random.normal(nu_b_hat_obs_mu3, 1., NSIG)  
  nuTil_s_pseu_sig[2] = np.random.normal(nu_s_hat_obs_mu3, 1., NSIG)  
  nuTil_b_pseu_sig[3] = np.random.normal(nu_b_hat_obs_mu4, 1., NSIG)  
  nuTil_s_pseu_sig[3] = np.random.normal(nu_s_hat_obs_mu4, 1., NSIG)  

  nuTil_b_pseu_sig[4] = np.random.normal(nu_b_hat_obs_mu5, 1., NSIG)  
  nuTil_s_pseu_sig[4] = np.random.normal(nu_s_hat_obs_mu5, 1., NSIG)  

  nuTil_b_pseu_sig[5] = np.random.normal(nu_b_hat_obs_mu6, 1., NSIG)  
  nuTil_s_pseu_sig[5] = np.random.normal(nu_s_hat_obs_mu6, 1., NSIG)  
  L_max_mu0 = max_L_iminuit_mu0(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0,0,0,0,0,0,0,0,0,0,0,0,1.)

  nu_b_hat_obs_0, nu_s_hat_obs_0 = L_max_mu0['nu_b'], L_max_mu0['nu_s']
  nu_b_hat_obs_02, nu_s_hat_obs_02 = L_max_mu0['nu_b2'], L_max_mu0['nu_s2']
  nu_b_hat_obs_03, nu_s_hat_obs_03 = L_max_mu0['nu_b3'], L_max_mu0['nu_s3']
  nu_b_hat_obs_04, nu_s_hat_obs_04 = L_max_mu0['nu_b4'], L_max_mu0['nu_s4']
  nu_b_hat_obs_05, nu_s_hat_obs_05 = L_max_mu0['nu_b5'], L_max_mu0['nu_s5']
  nu_b_hat_obs_06, nu_s_hat_obs_06 = L_max_mu0['nu_b6'], L_max_mu0['nu_s6']


  nu_b_hat_obs_0_v = [nu_b_hat_obs_0, nu_b_hat_obs_02, nu_b_hat_obs_03, nu_b_hat_obs_04, nu_b_hat_obs_05, nu_b_hat_obs_06]
  nu_s_hat_obs_0_v = [nu_s_hat_obs_0, nu_s_hat_obs_02, nu_s_hat_obs_03, nu_s_hat_obs_04, nu_s_hat_obs_05, nu_s_hat_obs_06]


  print "nu_b_hat_obs_0 "+ str(nu_b_hat_obs_0)
  nuTil_b_pseu_bkg[0] = np.random.normal(nu_b_hat_obs_0, 1., NBKG)  
  nuTil_s_pseu_bkg[0] = np.random.normal(nu_s_hat_obs_0, 1., NBKG)  
  nuTil_b_pseu_bkg[1] = np.random.normal(nu_b_hat_obs_02, 1., NBKG)  
  nuTil_s_pseu_bkg[1] = np.random.normal(nu_s_hat_obs_02, 1., NBKG)  
  nuTil_b_pseu_bkg[2] = np.random.normal(nu_b_hat_obs_03, 1., NBKG)  
  nuTil_s_pseu_bkg[2] = np.random.normal(nu_s_hat_obs_03, 1., NBKG)  
  nuTil_b_pseu_bkg[3] = np.random.normal(nu_b_hat_obs_04, 1., NBKG)  
  nuTil_s_pseu_bkg[3] = np.random.normal(nu_s_hat_obs_04, 1., NBKG)  
  nuTil_b_pseu_bkg[4] = np.random.normal(nu_b_hat_obs_05, 1., NBKG)  
  nuTil_s_pseu_bkg[4] = np.random.normal(nu_s_hat_obs_05, 1., NBKG)  
  nuTil_b_pseu_bkg[5] = np.random.normal(nu_b_hat_obs_06, 1., NBKG)  
  nuTil_s_pseu_bkg[5] = np.random.normal(nu_s_hat_obs_06, 1., NBKG)  



  for i in range(nr): 
   
   global minimize, brentq, fsolve
   from scipy.optimize import minimize, brentq, fsolve

  # Find observed test statistic
#   qMu_obs[i] = qMu(n[i], nuTil_b_obs, nuTil_s_obs, nuTil_b_obs, nuTil_s_obs, nuTil_s_obs, nuTil_b_obs, b0_in[i], db_in[i], s0_in[i], ds_in[i], 'obs',i)
  # print b0

  
 #  print nu_b_hat_obs_mu
  
 #  print nu_s_hat_obs_mu
  # Run pseudo experiments for the nuisance parameters. 
   qMu_pseu_sig_bla = []
   qMu_pseu_bkg_bla = []

  ###print "<<<<NSIG: "+str(NSIG)
  # - Retoss n_obs, assuming the most likely values for the nu  
  # - Retoss nu_tilde_observed assuming the most likely values for the nu
   NBKG = nPseudo
   nexp_sig[i] = nexp(1., nu_b_hat_obs_mu_v[i], nu_s_hat_obs_mu_v[i], b0_in[i], db_in[i], s0_in[i], ds_in[i])
   n_pseu_sig[i] = np.random.poisson(nexp_sig[i],NSIG)


   nexp_bkg[i] = nexp(0., nu_b_hat_obs_0_v[i], nu_s_hat_obs_0_v[i], b0_in[i], db_in[i], s0_in[i], ds_in[i])
   n_pseu_bkg[i] = np.random.poisson(nexp_bkg[i],NBKG)
      
   positive_sig = 0
   positive_bkg = 0
 
   
  n_pseu_sig = zip(*n_pseu_sig)
  print n_pseu_sig
  n_pseu_sig = np.asarray(n_pseu_sig)
  n_pseu_sig = n_pseu_sig.astype(float)
  qMu_tot_sig = [0 for a in range(NSIG)]
  qMu_tot_bkg = [0 for a in range(NBKG)]

  for i in range(NSIG):
    global n_p_sig
    n_p_sig = n_pseu_sig[i]
    print "n_p_sig= " + str(n_p_sig)
    L_max_mu1 = max_L_iminuit_mu1_pseu_sig(nuTil_b_pseu_sig[0][i], nuTil_s_pseu_sig[0][i], nuTil_b_pseu_sig[1][i], nuTil_s_pseu_sig[1][i], nuTil_b_pseu_sig[2][i], nuTil_s_pseu_sig[2][i], nuTil_b_pseu_sig[3][i], nuTil_s_pseu_sig[3][i], nuTil_b_pseu_sig[4][i], nuTil_s_pseu_sig[4][i], nuTil_b_pseu_sig[5][i], nuTil_s_pseu_sig[5][i],nuTil_s_pseu_sig[0][i], nuTil_b_pseu_sig[0][i], nuTil_s_pseu_sig[1][i], nuTil_b_pseu_sig[1][i], nuTil_s_pseu_sig[2][i], nuTil_b_pseu_sig[2][i], nuTil_s_pseu_sig[3][i], nuTil_b_pseu_sig[3][i], nuTil_s_pseu_sig[4][i], nuTil_b_pseu_sig[4][i], nuTil_s_pseu_sig[5][i], nuTil_b_pseu_sig[5][i],1. )
    L1 = f_value_mu1_pseu_sig
    L_max = max_L_iminuit_pseu_sig(nuTil_b_pseu_sig[0][i], nuTil_s_pseu_sig[0][i], nuTil_b_pseu_sig[1][i], nuTil_s_pseu_sig[1][i], nuTil_b_pseu_sig[2][i], nuTil_s_pseu_sig[2][i], nuTil_b_pseu_sig[3][i], nuTil_s_pseu_sig[3][i], nuTil_b_pseu_sig[4][i], nuTil_s_pseu_sig[4][i], nuTil_b_pseu_sig[5][i], nuTil_s_pseu_sig[5][i], L_max_mu1["nu_s"], L_max_mu1["nu_b"], L_max_mu1["nu_s2"], L_max_mu1["nu_b2"], L_max_mu1["nu_s3"], L_max_mu1["nu_b3"], L_max_mu1["nu_s4"], L_max_mu1["nu_b4"], L_max_mu1["nu_s5"], L_max_mu1["nu_b5"], L_max_mu1["nu_s6"], L_max_mu1["nu_b6"],1. )
    L2 = f_value_pseu_sig

    if L_max["mu"]>1:
       qMu_tot_sig[i] = 0
    else:
       qMu_tot_sig[i] = -2*(L2-L1)

     

    #qMu_tot_sig[i] = -2*(L2-L1)

    
    if qMu_tot_sig[i] < 0:
       print "Problem"
    if qMu_tot_sig[i]  >=  qMu_tot_obs:
      positive_sig += 1


  n_pseu_bkg = zip(*n_pseu_bkg)
  print n_pseu_bkg
  n_pseu_bkg = np.asarray(n_pseu_bkg)
  n_pseu_bkg = n_pseu_bkg.astype(float)

  for i in range(NBKG):
    global n_p_bkg
    print n_pseu_bkg[i]
    n_p_bkg = n_pseu_bkg[i]
    print "n_p_bkg= " + str(n_p_bkg)
    L_max_mu1 = max_L_iminuit_mu1_pseu_bkg(nuTil_b_pseu_bkg[0][i], nuTil_s_pseu_bkg[0][i], nuTil_b_pseu_bkg[1][i], nuTil_s_pseu_bkg[1][i], nuTil_b_pseu_bkg[2][i], nuTil_s_pseu_bkg[2][i], nuTil_b_pseu_bkg[3][i], nuTil_s_pseu_bkg[3][i], nuTil_b_pseu_bkg[4][i], nuTil_s_pseu_bkg[4][i], nuTil_b_pseu_bkg[5][i], nuTil_s_pseu_bkg[5][i],nuTil_s_pseu_bkg[0][i], nuTil_b_pseu_bkg[0][i], nuTil_s_pseu_bkg[1][i], nuTil_b_pseu_bkg[1][i], nuTil_s_pseu_bkg[2][i], nuTil_b_pseu_bkg[2][i], nuTil_s_pseu_bkg[3][i], nuTil_b_pseu_bkg[3][i], nuTil_s_pseu_bkg[4][i], nuTil_b_pseu_bkg[4][i], nuTil_s_pseu_bkg[5][i], nuTil_b_pseu_bkg[5][i],1. )
    L1 = f_value_mu1_pseu_bkg
    L_max = max_L_iminuit_pseu_bkg(nuTil_b_pseu_bkg[0][i], nuTil_s_pseu_bkg[0][i], nuTil_b_pseu_bkg[1][i], nuTil_s_pseu_bkg[1][i], nuTil_b_pseu_bkg[2][i], nuTil_s_pseu_bkg[2][i], nuTil_b_pseu_bkg[3][i], nuTil_s_pseu_bkg[3][i], nuTil_s_pseu_bkg[4][i], nuTil_b_pseu_bkg[4][i], nuTil_s_pseu_bkg[5][i], nuTil_b_pseu_bkg[5][i], L_max_mu1["nu_s"], L_max_mu1["nu_b"], L_max_mu1["nu_s2"], L_max_mu1["nu_b2"], L_max_mu1["nu_s3"], L_max_mu1["nu_b3"], L_max_mu1["nu_s4"], L_max_mu1["nu_b4"], L_max_mu1["nu_s5"], L_max_mu1["nu_b5"], L_max_mu1["nu_s6"], L_max_mu1["nu_b6"],1. )
    L2 = f_value_pseu_bkg

    qMu_tot_bkg[i] = -2*(L2-L1)

    
    if qMu_tot_bkg[i] < 0:
       print "Problem"
    if qMu_tot_bkg[i]  >=  qMu_tot_obs:
      positive_bkg += 1
   
  ps = positive_sig/float(NSIG)
   
  oneMinusPb = positive_bkg/float(NBKG)

  if oneMinusPb==0:
    CLs = 0

  else:
    CLs = ps/oneMinusPb
 
  del n_pseu_sig
  del nuTil_b_pseu_sig
  del nuTil_s_pseu_sig
  result = [qMu_tot_obs, ps, qMu_tot_sig, oneMinusPb, qMu_tot_bkg, CLs]
  result_save_2 = result[5]
  result_save_1 = result[0]
  print result_save_2

  with open ("cls_value.txt", "w") as f:
    f.write("%f %f" % (result_save_1, result_save_2))

  
  return (qMu_tot_obs, ps, qMu_tot_sig, oneMinusPb, qMu_tot_bkg, CLs) 
"""
if __name__ == '__main__':
  if len(sys.argv) == 6:
    result = cls_obs(float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), 25000)    
    print "CLs: "+str(result[0])+" +- "+str(result[1])    
  else :
    print "Syntax for CLs:   ./cls.py   nobs   bkg   bkgerr   sig   sigerr" 
"""

if __name__ == '__main__':
  '''
  obsa= 5
  obsa= 10
  ba= 3
  ba= 4
  berra= 0
  berra= 0
  sa= 1
  sa= 2
  serra= 0
  serra= 0



  '''
#  n = [10]
#  nuTil = [0,0,0]
#  nPseudo = 10
#  obs = []
  #b0 = [9]
  #db = [1]
  #s0 = [7]
  #ds = [1]
  
  global minimize, brentq, fsolve
  from scipy.optimize import minimize, brentq, fsolve
  global b0, db, s0,ds, n
  b0 =[90.0,90.0]
  db=[10.0, 10.0]
  s0=[30.0, 30.0]
  ds =[10., 10.0]
  n = [100.0, 100.0]
#def calc_CLs(n_obs, nuTil_b_obs, b0_in, db_in, nuTil_s_obs, s0_in, ds_in, nPseudo,nr_Sr,i):
#  result = qMu(14,0.5,0.5)

#def qMu(n, nuTil_b_mu, nuTil_s_mu, nuTil_b, nuT_s, mu_hat):  
#def max_L_iminuit(n,b,da,s,ds,nu_b_til, nu_s_til, nu_s, nu_b, mu):
# result = qMu(150.0, -0.4, -0.4, )
  result = calc_CLs(n, 0.0, b0, db, 0.0, s0, ds, 2500,2)
##  result= max_L_iminuit(0 ,0,0,0,1)
#  result=max_L_iminuit()
  #result=qMu(10, 0,0 )
  print result
  print result[5]
  #result = calc_CLs(n, nuTil,b, berr, nuTil,s, serr, nPseudo, 3)
  #print result 
"""
  if len(sys.argv) == 6:
    result = cls_obs(float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), 25000)    
    print "CLs: "+str(result[0])+" +- "+str(result[1])    
  else :
    print "Syntax for CLs:   ./cls.py   nobs   bkg   bkgerr   sig   sigerr" 
"""

