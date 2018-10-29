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

def func_to_max(nu_b_til, nu_s_til,nu_b_til2, nu_s_til2, nu_b_til3, nu_s_til3, nu_b_til4, nu_s_til4, nu_b_til5, nu_s_til5, nu_s_tilsys, mu, nu_s, nu_b,nu_s2, nu_b2, nu_s3, nu_b3, nu_s4, nu_b4, nu_s5, nu_b5, nu_sys):
  xx=n
  result =0

#  print "n=" + str(n) 
#  print "b=" + str(b0) 
#  print "s=" + str(s0) 
#  print "db=" + str(ds) 
#  print "db=" + str(db) 
#  print "n=" + str(n) 
  zwischen = [0 for i in range(len(xx))]
  lam = [0 for i in range(len(xx))]
  nu_b = [nu_b, nu_b2, nu_b3, nu_b4, nu_b5]
  nu_s = [nu_s, nu_s2, nu_s3, nu_s4, nu_s5]
  nu_b_til = [nu_b_til, nu_b_til2, nu_b_til3, nu_b_til4, nu_b_til5]
  nu_s_til = [nu_s_til, nu_s_til2, nu_s_til3, nu_s_til4, nu_s_til5]
  for i in range(len(xx)):
    try:
      if s0[i] == 0.0:
        lam[i]= b0[i]*exp(db[i]/b0[i]*nu_b[i])
      else: 
        lam[i]= b0[i]*exp(db[i]/b0[i]*nu_b[i]) + mu*s0[i]*exp(ds[i]/s0[i]*nu_s[i])*exp(ds_sys[i]/s0[i]*nu_sys)
    except OverflowError as e:
      lam[i]=100000000
  
  for i in range(len(xx)):
#    try:
#    print lam[i]
#    print xx[i]
#    print poisson.pmf(xx[i], lam[i])
    try:
      zwischen[i] = -log(poisson.pmf(xx[i], lam[i] )*exp(-(nu_b_til[i]-nu_b[i])**2/2.)*exp(-(nu_s_til[i]-nu_s[i])**2/2.))
    except ValueError as e:
      if xx[i] == 0.001:
        rama = xx[i]* log(xx[i]) - xx[i] + ( log(xx[i]*(1+4*xx[i]*(1+2*xx[i])))/6.  )+ (log(pi)/2. )
        zwischen[i] = - xx[i] * log(lam[i]) + rama + lam[i] + log( exp(-(nu_b_til[i]-nu_b[i])**2/2.)*exp(-(nu_s_til[i]-nu_s[i])**2/2.) )
      else:
        zwischen[i] = log(factorial(xx[i])) - xx[i]*log(lam[i]) + lam[i]+log( exp(-(nu_b_til[i]-nu_b[i])**2/2.)*exp(-(nu_s_til[i]-nu_s[i])**2/2.))

#      rama = xx[i]* log(xx[i]) - xx[i] + ( log(xx[i]*(1+4*xx[i]*(1+2*xx[i])))/6.  )+ (log(pi)/2. )
#      zwischen[i] = - xx[i] * log(lam[i]) + rama + lam[i] + log( exp(-(nu_b_til[i]-nu_b[i])**2/2.)*exp(-(nu_s_til[i]-nu_s[i])**2/2.) )

    
    result += zwischen[i]
  result += (nu_s_tilsys - nu_sys)**2/2.
  return result




def max_L_iminuit(nu_b_til, nu_s_til, nu_b_til2, nu_s_til2, nu_b_til3, nu_s_til3, nu_b_til4, nu_s_til4, nu_b_til5, nu_s_til5, nu_s_tilsys, nu_s, nu_b, nu_s2, nu_b2, nu_s3, nu_b3, nu_s4, nu_b4, nu_s5, nu_b5, nu_sys, mu):
  m=Minuit(func_to_max, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True,nu_s_til2= nu_s_til2, fix_nu_s_til2=True,  nu_b_til2= nu_b_til2, fix_nu_b_til2=True,nu_s_til3= nu_s_til3, fix_nu_s_til3=True,  nu_b_til3= nu_b_til3, fix_nu_b_til3=True,nu_s_til4= nu_s_til4, fix_nu_s_til4=True,  nu_b_til4= nu_b_til4, fix_nu_b_til4=True,nu_s_til5= nu_s_til5, fix_nu_s_til5=True,  nu_b_til5= nu_b_til5, fix_nu_b_til5=True,nu_s_tilsys= nu_s_tilsys, fix_nu_s_tilsys=True, mu=0.0, fix_mu = False, limit_mu=(0,6), error_mu = 1.0,  nu_s= nu_s, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-7,7), nu_b= nu_b, fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-7,7),nu_s2= nu_s2, fix_nu_s2=False, error_nu_s2=1.0, limit_nu_s2=(-7,7), nu_b2= nu_b2, fix_nu_b2=False, error_nu_b2= 1.0, limit_nu_b2=(-7,7),nu_s3= nu_s3, fix_nu_s3=False, error_nu_s3=1.0, limit_nu_s3=(-7,7), nu_b3= nu_b3, fix_nu_b3=False, error_nu_b3= 1.0, limit_nu_b3=(-7,7),nu_s4= nu_s4, fix_nu_s4=False, error_nu_s4=1.0, limit_nu_s4=(-7,7), nu_b4= nu_b4, fix_nu_b4=False, error_nu_b4= 1.0, limit_nu_b4=(-7,7),nu_s5= nu_s5, fix_nu_s5=False, error_nu_s5=1.0, limit_nu_s5=(-7,7), nu_b5= nu_b5, fix_nu_b5=False, error_nu_b5= 1.0, limit_nu_b5=(-7,7), nu_sys= 0.0, fix_nu_sys=False, error_nu_sys= 1.0, limit_nu_sys=(-7,7),print_level=0)
  m.migrad()
  m.minos()
  m.hesse()
  m.migrad()
  global f_value
  f_value = m.fval
#  print 'Errors=' + str(m.errors)
#  print m.get_merrors()
  print m.get_fmin()
  if m.get_fmin().is_valid == False:
    print 'Fehler_valid'
  if m.matrix_accurate() == False:
    print 'Fehler_accurate'
  Errors= (m.get_merrors())
#  if Errors["nu_b"]["is_valid"] == False:
#    print "Fehler_nu_b_valid"
#  if Errors["nu_s"]["is_valid"] == False:
#    print "Fehler_nu_s_valid"
#  if Errors["mu"]["is_valid"] == False:
#    print "Fehler_mu_valid"
#  f_value = m.args
#  print('args',m.args)
#  print('value', m.values)
#  print "Funktionswert" + str(func(*f_value))
  return m.values

def max_L_iminuit_mu1(nu_b_til, nu_s_til, nu_b_til2, nu_s_til2,nu_b_til3, nu_s_til3,nu_b_til4, nu_s_til4, nu_s_til5, nu_b_til5, nu_s_tilsys,nu_s, nu_b, nu_s2, nu_b2, nu_s3, nu_b3,nu_s4, nu_b4, nu_s5, nu_b5, nu_sys,mu):
  m=Minuit(func_to_max, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, nu_s_til2= nu_s_til2, fix_nu_s_til2=True, nu_b_til2= nu_b_til2, fix_nu_b_til2=True,nu_s_til3= nu_s_til3, fix_nu_s_til3=True, nu_b_til3= nu_b_til3, fix_nu_b_til3=True,nu_s_til4= nu_s_til4, fix_nu_s_til4=True, nu_b_til4= nu_b_til4, fix_nu_b_til4=True,nu_s_til5= nu_s_til5, fix_nu_s_til5=True, nu_b_til5= nu_b_til5, fix_nu_b_til5=True, nu_s_tilsys= nu_s_tilsys, fix_nu_s_tilsys=True, mu=1.0, fix_mu = True, nu_s = nu_s_til, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-7,7), nu_b= nu_b_til, fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-7,7), nu_s2 = nu_s_til2, fix_nu_s2=False, error_nu_s2=1.0, limit_nu_s2=(-7,7), nu_b2= nu_b_til2, fix_nu_b2=False, error_nu_b2= 1.0, limit_nu_b2=(-7,7),nu_s3 = nu_s_til3, fix_nu_s3=False, error_nu_s3=1.0, limit_nu_s3=(-7,7), nu_b3= nu_b_til3, fix_nu_b3=False, error_nu_b3= 1.0, limit_nu_b3=(-7,7),nu_s4 = nu_s_til4, fix_nu_s4=False, error_nu_s4=1.0, limit_nu_s4=(-7,7), nu_b4= nu_b_til4, fix_nu_b4=False, error_nu_b4= 1.0, limit_nu_b4=(-7,7),nu_s5 = nu_s_til5, fix_nu_s5=False, error_nu_s5=1.0, limit_nu_s5=(-7,7), nu_b5= nu_b_til5, fix_nu_b5=False, error_nu_b5= 1.0, limit_nu_b5=(-7,7), nu_sys= nu_s_tilsys, fix_nu_sys=False, error_nu_sys= 1.0, limit_nu_sys=(-7,7),print_level=0)
  m.migrad()
  m.minos()
  m.hesse()
  m.migrad()
  dic =m.get_fmin()

  if dic["has_accurate_covar"]== False:
    print "Fehler"

  print 'Errors= ' + str(m.errors)
#  print m.get_merrors()
  if m.get_fmin().is_valid == False:
    print 'Fehler_valid'
  if m.matrix_accurate() == False:
    print 'Fehler_accurate'
  global f_value_mu1
  f_value_mu1 = m.fval
  Errors= (m.get_merrors())
#  if Errors["nu_b"]["is_valid"] == False:
#    print "Fehler_nu_b"
#  if Errors["nu_s"]["is_valid"] == False:
#    print "Fehler_nu_s"
#  global f_value_mu1
#  f_value_mu1 = m.fval
#  print('parameters', m.parameters)
#  print('args',m.args)
#  print('value', m.values)

  return m.values





def calc_CLs(n_obs, nuTil_b_obs, b0_in, db_in, nuTil_s_obs, s0_in, ds_in, ds_insys, nPseudo, nr_SR):
#def calc_CLs( b0_in, db_in, s0_in, ds_in):
  """ Calcuates the CLs value given the measurement and signal/background estimates by running nPseudo experiments """

  nr = nr_SR

  global b0, s0, db, ds, ds_sys, n # Passes this to nexp function
  
  b0 = b0_in
  s0 = s0_in
  db = db_in
  ds = ds_in
  ds_sys = ds_insys
  n = n_obs

  high_L1 = False
  print "n=" + str(n)
  try:
    L_max_mu1 = max_L_iminuit_mu1(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0,0,0,0,0,0,0.,0.,0.,0.,1.)
    L1_obs = f_value_mu1
  except:
    L1_obs = 500
  print "L1_obs= " + str(L1_obs)
  if L1_obs > 1000000000: 
    high_L1 = True 
  #  print "HILFE"
  #  print n
  #  print b0
  #  print db
  #  print s0
  #  print ds
  #  print ds_sys

  try:
    L_max = max_L_iminuit(0,0,0,0,0,0,0,0,0,0,0,L_max_mu1["nu_s"],L_max_mu1["nu_b"],L_max_mu1["nu_s2"],L_max_mu1["nu_b2"],L_max_mu1["nu_s3"],L_max_mu1["nu_b3"],L_max_mu1["nu_s4"],L_max_mu1["nu_b4"],L_max_mu1["nu_s5"],L_max_mu1["nu_b5"],L_max_mu1["nu_sys"],1.)
    L2_obs = f_value
  except:
    L2_obs= 500
  print "L2_obs= " + str(L2_obs)

  if L2_obs == 500 or L1_obs == 500: 
    qMu_tot_obs = 500
  else:
    qMu_tot_obs = (-2.*(L2_obs-L1_obs))
  print "qMu_tot_obs= " + str(qMu_tot_obs)    

  return (qMu_tot_obs, -L1_obs, L_max["mu"], high_L1, 0, 0) 
