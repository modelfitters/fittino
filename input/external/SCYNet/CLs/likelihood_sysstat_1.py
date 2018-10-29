#!/usr/bin/env python



from math import *
import os, sys, scipy
import numpy as np
from math import log, exp, sqrt, factorial, pi
import os, sys
import scipy.optimize as opt
#import cls_nocov as CN 
from scipy.stats import poisson
from iminuit import Minuit, describe, Struct
# Note: Scipy packages are imported within the code, using a try-except setup

def nuWithMaxL(n, nuTil_b, nuTil_s, mu, b0, db, s0, ds):
  """ Finds root of eqForMaxL to maximise Likelihood   """  
  global fsolve
  return fsolve(eqForMaxL, (nuTil_b, nuTil_s), (n, nuTil_b, nuTil_s, mu, b0, db, s0, ds))
 

def eqForMaxL((nu_b, nu_s), n, nuTil_b, nuTil_s, mu, b0, db, s0, ds):
  """
  Helper function to maximise likelihood 
  The nu_b, nu_s for which these functions are zero maximise the likelihood 
  """
  print b0
  print s0
  lam = nexp(mu, nu_b, nu_s, b0, db, s0, ds)
  return ((n-lam)*db*exp(db/b0*nu_b)-lam*(nu_b-nuTil_b), (n-lam)*ds*mu*exp(ds/s0*nu_s)-lam*(nu_s-nuTil_s))
  #return ((n-lam)*db*exp(db/b0*nu_b)-lam*(nu_b-nuTil_b), db*exp(db/b0*nu_b-ds/s0*nu_s)*(nu_s-nuTil_s)-ds*mu*(nu_b-nuTil_b))
 

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
  print b
  print ddb
  print s
  print dds
  
  return b*exp(ddb/b*nu_b) + mu*s*exp(dds/s*nu_s)


def reduced_logp(nuTil, nu):
  """ reduced_logp = log(p) without normalisation """
  return -(nuTil-nu)**2./2.

def func_to_max(nu_b_til, nu_s_til,mu, nu_s, nu_b):
  """ describes the function which is minimized for the calculation of q_Mu_obs"""

  xx=n
  result =0
  zwischen = [0 for i in range(len(xx))]
  lam = [0 for i in range(len(xx))]
  for i in range(len(xx)):
    print b0[i]
    print s0[i]
    print db[i]
    print ds[i]
    print xx[i]
    try:
      if s0[i] == 0.0:
        lam[i]= b0[i]*exp(db[i]/b0[i]*nu_b) 
      else:
        lam[i]= b0[i]*exp(db[i]/b0[i]*nu_b) + mu*s0[i]*exp(ds[i]/s0[i]*nu_s)
        
    except OverflowError as e:
      lam[i]=100000000

  for i in range(len(xx)):
#    try:
    print lam[i]
    print xx[i]
    print poisson.pmf(xx[i], lam[i])
    try:
      zwischen[i] = -log(poisson.pmf(xx[i], lam[i] ))
    except ValueError as e:
      try:
        zwischen[i] = log(factorial(xx[i])) - xx[i]*log(lam[i]) + lam[i]
      except:
        rama = xx[i]* log(xx[i]) - xx[i] + ( log(xx[i]*(1+4*xx[i]*(1+2*xx[i])))/6.  )+ (log(pi)/2. )
        zwischen[i] = - xx[i] * log(lam[i]) + rama + lam[i] 





    
    result += zwischen[i]
  if s0[i] == 0.0:
    result += (nu_b_til-nu_b)**2/2. 
  else:
    result += (nu_b_til-nu_b)**2/2. + (nu_s_til- nu_s)**2/2. 
  return result


def func_to_max_pseu_sig(nu_b_til, nu_s_til,mu, nu_s, nu_b):
  """ describes the function which is minimized for the calculation in the toy experiments""" 
  xx=n_p_sig
  result =0 
  zwischen = [0 for i in range(len(xx))]
  lam = [0 for i in range(len(xx))]
  for i in range(len(xx)):
    
    lam[i]= b0[i]*(exp( db[i]/b0[i]*nu_b)) + mu*s0[i]*(exp( ds[i]/s0[i]*nu_s ))
    print "lam= " + str(lam[i])
    print "xx= " + str(xx[i])
    print "poisson= " + str(poisson.pmf(xx[i], lam[i]))
  #  zwischen[i] = -log(poisson.pmf(xx[i], lam[i]) * (exp(-nu_b_til**2/2.) + nu_b*nu_b_til*exp(-nu_b_til**2/2.)) * (exp(-nu_s_til**2/2.) + nu_s*nu_s_til*exp(-nu_s_til**2/2.)))
    try:
      zwischen[i] = -log(poisson.pmf(xx[i], lam[i]) )
    except:
      rama = xx[i]* log(xx[i]) - xx[i] + ( log(xx[i]*(1+4*xx[i]*(1+2*xx[i])))/6.  )+ (log(pi)/2. )
      zwischen[i] = - xx[i] * log(lam[i]) + rama + lam[i]
    result += zwischen[i]
  result +=  (nu_b_til-nu_b)**2/2. + (nu_s_til- nu_s)**2/2.
  return result

def func_to_max_pseu_bkg(nu_b_til, nu_s_til,mu, nu_s, nu_b):
  """ describes the function which is minimized for the calculation in the toy experiments""" 
  xx=n_p_bkg
  result =0 
  zwischen = [0 for i in range(len(xx))]
  lam = [0 for i in range(len(xx))]
  for i in range(len(xx)):
    print "xx= " + str(xx[i])
    print "lam= " + str(lam[i])
    print "poisson= " + str(poisson.pmf(xx[i], lam[i]))
    lam[i]= b0[i]*(exp( db[i]/b0[i]*nu_b)) + mu*s0[i]*(exp( ds[i]/s0[i]*nu_s ))
  #  zwischen[i] = -log(poisson.pmf(xx[i], lam[i]))
  #  zwischen[i] = -log(poisson.pmf(xx[i], lam[i]) * (exp(-nu_b_til**2/2.) + nu_b*nu_b_til*exp(-nu_b_til**2/2.)) * (exp(-nu_s_til**2/2.) + nu_s*nu_s_til*exp(-nu_s_til**2/2.)))
    try:
      zwischen[i] = -log(poisson.pmf(xx[i], lam[i]) )
    except:
      rama = xx[i]* log(xx[i]) - xx[i] + ( log(xx[i]*(1+4*xx[i]*(1+2*xx[i])))/6.  )+ (log(pi)/2. )
      zwischen[i] = - xx[i] * log(lam[i]) + rama + lam[i]
    result += zwischen[i]
  result +=  (nu_b_til-nu_b)**2/2. + (nu_s_til- nu_s)**2/2.
  return result



def max_L_iminuit(nu_b_til, nu_s_til, nu_s, nu_b, mu):
  m=Minuit(func_to_max, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, mu= 0 , fix_mu = False, limit_mu=(0,6), error_mu = 1.0,  nu_s= 0, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-10,10), nu_b= 0, fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-10,10), print_level=2)
  m.migrad()
  m.hesse()
  m.minos()
  m.migrad()
  global f_value
  f_value = m.fval
  checks = True
  print 'Errors=' + str(m.errors)
  print m.get_merrors()
  print m.get_fmin()
  if m.get_fmin().is_valid == False:
    print 'Fehler'
    checks = False
  if m.matrix_accurate() == False:
    print 'Fehler'
    checks = False
  Errors= (m.get_merrors())
#  if Errors["nu_b"]["is_valid"] == False:
#    print "Fehler_nub"
#  if Errors["nu_s"]["is_valid"] == False:
#    print "Fehler_nus"
  if Errors["mu"]["is_valid"] == False:
    print "Fehler"

  return (m.values, checks)


def max_L_iminuit_mu1(nu_b_til, nu_s_til, nu_s, nu_b, mu):
  m=Minuit(func_to_max, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, mu=1.0, fix_mu = True, nu_s = nu_s, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-10,10), nu_b= nu_b, fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-10,10), print_level=2)
  m.migrad()
  m.minos()
  m.hesse()
  m.migrad()
  dic =m.get_fmin()
  checks = True
  if dic["has_accurate_covar"]== False:
    print "Fehler"
    
#  print 'Errors= ' + str(m.errors)
#  print m.get_merrors()
  if m.get_fmin().is_valid == False:
    print 'Fehler'
    checks = False
  if m.matrix_accurate() == False:
    print 'Fehler'
    checks = False
  global f_value_mu1
  f_value_mu1 = m.fval
  Errors= (m.get_merrors())
  if Errors["nu_b"]["is_valid"] == False:
    print "Fehler_nub"
    checks=False
  if Errors["nu_s"]["is_valid"] == False:
    print "Fehler_nus"
    checks = False
  return (m.values, checks)


def max_L_iminuit_mu0(nu_b_til, nu_s_til, nu_s, nu_b, mu):
  m=Minuit(func_to_max, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, mu=0.0, fix_mu = True, nu_s = nu_s, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-3,3), nu_b= nu_b, fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-3,3), print_level=2)
  m.migrad()
  m.minos()
  m.hesse()
  m.migrad()
  dic =m.get_fmin()
  
  if dic["has_accurate_covar"]== False:
    print "Fehler"
    
  print 'Errors= ' + str(m.errors)
  print m.get_merrors()
  if m.get_fmin().is_valid == False:
    print 'Fehler'
  if m.matrix_accurate() == False:
    print 'Fehler'
  global f_value_mu0
  f_value_mu0 = m.fval
  Errors= (m.get_merrors())
  if Errors["nu_b"]["is_valid"] == False:
    print "Fehler_nub"
  if Errors["nu_s"]["is_valid"] == False:
    print "Fehler_nus"
  return m.values



def max_L_iminuit_pseu_sig(nu_b_til, nu_s_til, nu_s, nu_b, mu):
  m=Minuit(func_to_max_pseu_sig, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, mu=mu, fix_mu = False, limit_mu=(0,4), error_mu= 1.0,nu_s= nu_s, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-3,3), nu_b= nu_b, fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-3,3), print_level=3)
  m.migrad()
  m.minos()
  m.hesse()
  m.migrad()
#  print m.get_fmin()
  dic =m.get_fmin()
  
  if dic["has_accurate_covar"]== False:
    print "Fehler"
  print 'Errors= ' + str(m.errors)
  print m.get_merrors()
  if m.get_fmin().is_valid == False:
    print 'Fehler'
  if m.matrix_accurate() == False:
    print 'Fehler'
  global f_value_pseu_sig
  f_value_pseu_sig = m.fval
  Errors= (m.get_merrors())
  if Errors["nu_b"]["is_valid"] == False:
    print "Fehler"
  if Errors["nu_s"]["is_valid"] == False:
    print "Fehler"
  if Errors["mu"]["is_valid"] == False:
    print "Fehler"
  return m.values

def max_L_iminuit_mu1_pseu_sig(nu_b_til, nu_s_til, nu_s, nu_b, mu):
  m=Minuit(func_to_max_pseu_sig, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, mu=1.0, fix_mu = True, nu_s= nu_s, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-3,3), nu_b= nu_b , fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-3,3), print_level=3)
  m.migrad()
  m.minos()
  m.hesse()
  m.migrad() 
#  print m.get_fmin()
  print 'Errors ' + str(m.errors)
  print m.get_merrors()
  dic =m.get_fmin()
  
  if dic["has_accurate_covar"]== False:
    print "Fehler"
  if m.get_fmin().is_valid == False:
    print 'Fehler'
  if m.matrix_accurate() == False:
    print 'Fehler'
  global f_value_mu1_pseu_sig
  f_value_mu1_pseu_sig = m.fval
  Errors= (m.get_merrors())
  if Errors["nu_b"]["is_valid"] == False:
    print "Fehler"
  if Errors["nu_s"]["is_valid"] == False:
    print "Fehler"
  return m.values


def max_L_iminuit_pseu_bkg(nu_b_til, nu_s_til, nu_s, nu_b, mu):
  m=Minuit(func_to_max_pseu_bkg, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, mu=mu, fix_mu = False, limit_mu=(0,4), error_mu= 1.0,nu_s= nu_s, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-3,3), nu_b= nu_b, fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-3,3), print_level=3)
  m.migrad()
  m.minos()
  m.hesse()
  m.migrad()
  print m.get_fmin()
  dic =m.get_fmin()
  
  if dic["has_accurate_covar"]== False:
    print "Fehler"
  print 'Errors= ' + str(m.errors)
  print m.get_merrors()
  if m.get_fmin().is_valid == False:
    print 'Fehler'
  if m.matrix_accurate() == False:
    print 'Fehler'
  global f_value_pseu_bkg
  f_value_pseu_bkg = m.fval
  Errors= (m.get_merrors())
  if Errors["nu_b"]["is_valid"] == False:
    print "Fehler"
  if Errors["nu_s"]["is_valid"] == False:
    print "Fehler"
  if Errors["mu"]["is_valid"] == False:
    print "Fehler"
  return m.values




def max_L_iminuit_mu1_pseu_bkg(nu_b_til, nu_s_til, nu_s, nu_b, mu):
  m=Minuit(func_to_max_pseu_bkg, nu_s_til= nu_s_til, fix_nu_s_til=True,  nu_b_til= nu_b_til, fix_nu_b_til=True, mu=1.0, fix_mu = True, nu_s= nu_s, fix_nu_s=False, error_nu_s=1.0, limit_nu_s=(-3,3), nu_b= nu_b , fix_nu_b=False, error_nu_b= 1.0, limit_nu_b=(-3,3), print_level=3)
  m.migrad()
  m.minos()
  m.hesse()
  m.migrad() 
  print m.get_fmin()
  print 'Errors ' + str(m.errors)
  print m.get_merrors()
  dic =m.get_fmin()
  
  if dic["has_accurate_covar"]== False:
    print "Fehler"
  if m.get_fmin().is_valid == False:
    print 'Fehler'
  if m.matrix_accurate() == False:
    print 'Fehler'
  global f_value_mu1_pseu_bkg
  f_value_mu1_pseu_bkg = m.fval
  Errors= (m.get_merrors())
  if Errors["nu_b"]["is_valid"] == False:
    print "Fehler"
  if Errors["nu_s"]["is_valid"] == False:
    print "Fehler"
  return m.values



def calc_CLs(n_obs, nuTil_b_obs, b0_in, db_in, nuTil_s_obs, s0_in, ds_in, nPseudo, nr_SR, nu_b, nu_s, mu):

  nr = nr_SR
  NBKG = nPseudo
  NSIG = nPseudo * 5 # You need more signal MC as only few will pass
  global b0, s0, db, ds, n # Passes this to nexp function

  b0 = b0_in
  s0 = s0_in
  db = db_in
  ds = ds_in
  n = n_obs
  for i in range(len(n)): 
    if n[i] ==0: n[i] = 0.0001
  """ Calculation of q_Mu_obs and throwing dices for nuicances parameters"""

  qMu_obs=[0 for a in range(nr)]

  nexp_sig = [0 for a in range (nr)]
  n_pseu_sig = [0 for a in range(nr)]

  nexp_bkg = [0 for a in range (nr)]
  n_pseu_bkg = [0 for a in range(nr)]

  #L_max_mu1 = max_L_iminuit_mu1( nuTil_b_obs, nuTil_s_obs, nu_s, nu_b,1.)
  #nu_b_hat_obs_mu, nu_s_hat_obs_mu = L_max_mu1['nu_b'], L_max_mu1['nu_s']


  #L_max_mu0 = max_L_iminuit_mu0( nuTil_b_obs, nuTil_s_obs, nu_s, nu_b,0.)
  #nu_b_hat_obs_0, nu_s_hat_obs_0 = L_max_mu0['nu_b'], L_max_mu0['nu_s']


  nuTil_b_pseu_sig = [0 for i in range(NSIG)]
  nuTil_s_pseu_sig = [0 for i in range(NSIG)]

  #nuTil_b_pseu_sig = np.random.normal(nu_b_hat_obs_mu, 1., NSIG)  
  #nuTil_s_pseu_sig = np.random.normal(nu_s_hat_obs_mu, 1., NSIG)  


  nuTil_b_pseu_bkg = [0 for i in range(NBKG)]
  nuTil_s_pseu_bkg = [0 for i in range(NBKG)]

  #nuTil_b_pseu_bkg = np.random.normal(nu_b_hat_obs_0, 1., NBKG)  
  #nuTil_s_pseu_bkg = np.random.normal(nu_s_hat_obs_0, 1., NBKG)  

  """nominator of qMu_obs """
  try:
    L_max_mu1 = max_L_iminuit_mu1(nuTil_b_obs, nuTil_s_obs, nu_b, nu_s, 1.0)
    L1_obs = f_value_mu1
    print "L1=" + str(L1_obs)
    print L_max_mu1
    if L_max_mu1[1]==False:
      L1_obs = 500.
  
      """denominator of qMu_obs """
  except:
    L1_obs = 500.
  try:
    L_max = max_L_iminuit(nuTil_b_obs,nuTil_s_obs, L_max_mu1[0]["nu_b"], L_max_mu1[0]["nu_s"],0.0)
    L2_obs = f_value
    print "L2=" + str(L2_obs)
    if L_max[1] ==False:
      L2_obs = 500.
  except:
    L2_obs = 500.0
  if L2_obs == 500.0 or L1_obs==500.0:
    qMu_tot_obs = 500.
#    qMu_tot_obs = (-2.*(L2_obs-L1_obs))
  else:
    qMu_tot_obs = (-2.*(L2_obs-L1_obs))
  
  print "qMu_obs=" + str(qMu_tot_obs)
  if qMu_tot_obs < 0.0:
     print "Achtung"
     print (b0)
     print (db)
     print (s0)
     print (ds)
     print (n)

  """thorwing dices for n_exp, individual per SR """

#  for i in range(nr):
#    nexp_sig[i] = nexp(1., nu_b_hat_obs_mu, nu_s_hat_obs_mu, b0_in[i], db_in[i], s0_in[i], ds_in[i])
#    n_pseu_sig[i] = np.random.poisson(nexp_sig[i], NSIG)


#  for i in range(nr):
#    nexp_bkg[i] = nexp(0., nu_b_hat_obs_0, nu_s_hat_obs_0, b0_in[i], db_in[i], s0_in[i], ds_in[i])
#    n_pseu_bkg[i] = np.random.poisson(nexp_bkg[i], NBKG)


  """Calculation of the Test Statistics through Toy Experiments"""
  '''
  positive_sig = 0
  positive_bkg = 0
  n_pseu_sig = zip(*n_pseu_sig)
  print n_pseu_sig
  n_pseu_sig = np.asarray(n_pseu_sig)
  qMu_tot_sig = [0 for a in range(NSIG)]


  n_pseu_bkg = zip(*n_pseu_bkg)
  print n_pseu_bkg
  n_pseu_bkg = np.asarray(n_pseu_bkg)
  qMu_tot_bkg = [0 for a in range(NBKG)]


  for i in range(NSIG):
      """lists of N_exp for all SR per toy experiment i """

      global n_p_sig
      n_p_sig = n_pseu_sig[i]
      """ Calulation of toy experiments"""
 #     if n_p_sig[0] ==0:
 #       L1 = 0
 #       L2 = 0
 #       L_max["mu"] = 4
  #    else:

      L_max_mu1 = max_L_iminuit_mu1_pseu_sig(nuTil_b_pseu_sig[i],nuTil_s_pseu_sig[i], nu_b, nu_s,1.0)
      L1 = f_value_mu1_pseu_sig
      try:
        L_max = max_L_iminuit_pseu_sig(nuTil_b_pseu_sig[i],nuTil_s_pseu_sig[i], nu_b, nu_s,mu)
        L2 = f_value_pseu_sig
      except:
        L2 = 500.

      if L2 == 500:
        qMu_tot_sig[i] = 500.

      else:
        qMu_tot_sig[i]  = -2.*(L2-L1)

      print "qMu_tot_sig =" + str(qMu_tot_sig[i])


      if qMu_tot_sig[i] < 0:
         print "Problem"
   #    print "n_pseu_sig =" + str(n_pseu_sig[i])
   #    print "nuTil_b_pseu_sig[i]" + str(nuTil_b_pseu_sig[i])
   #    print "nuTil_s_pseu_sig[i]" + str(nuTil_s_pseu_sig[i])
   #    print "n_p_sig" + str(n_p_sig)
  #     test = CN.qMu(n_pseu_sig[i][0], nuTil_b_pseu_sig[i], nuTil_s_pseu_sig[i], b0_in[0], db_in[0], s0_in[0], ds_in[0])
  #     print "qMu1_nocov=" + str(test)

   #    print "qMu2=" + str(qMu(n_pseu_sig[i][1], nuTil_b_pseu_sig[i], nuTil_s_pseu_sig[i], b0_in[1], db_in[1], s0_in[1], ds_in[1], 'pseu_sig'))
   #    print "n_pseu_sig =" + str(n_pseu_sig[i][1])
   #    print "nuTil_b_pseu_sig[i]" + str(nuTil_b_pseu_sig[i])
   #    print "nuTil_s_pseu_sig[i]" + str(nuTil_s_pseu_sig[i])
   #    print "n_p_sig" + str(n_p_sig)

   #    test1 = CN.qMu(n_pseu_sig[i][1], nuTil_b_pseu_sig[i], nuTil_s_pseu_sig[i], b0_in[1], db_in[1], s0_in[1], ds_in[1])
   #    print "qMu2_nocov=" + str(test1)
      if qMu_tot_sig[i]  >=  qMu_tot_obs:
          positive_sig += 1


  for i in range(NBKG):
      """lists of N_exp for all SR per toy experiment i """

      global n_p_bkg
      n_p_bkg = n_pseu_bkg[i]
      
      print "n_p_bkg= " + str(n_p_bkg)
      """ Calulation of toy experiments"""
     # if n_p_bkg[0]==0:
     #   L1 = 0
     #   L2 = 0
     #   L_max["mu"] = 4
     # else:

      L_max_mu1 = max_L_iminuit_mu1_pseu_bkg(nuTil_b_pseu_bkg[i],nuTil_s_pseu_bkg[i], nu_b, nu_s,1.0)
      L1 = f_value_mu1_pseu_bkg
      L_max = max_L_iminuit_pseu_bkg(nuTil_b_pseu_bkg[i],nuTil_s_pseu_bkg[i], nu_b, nu_s,mu)
      L2 = f_value_pseu_bkg

     # if L_max["mu"] > 1:
     #   qMu_tot_bkg[i] = 0

     # else:
      qMu_tot_bkg[i]  = -2*(L2-L1)

      print "qMu_tot_bkg =" + str(qMu_tot_bkg[i])


      if qMu_tot_bkg[i]  >=  qMu_tot_obs:
          positive_bkg += 1
  '''
#  oneMinusPb = positive_bkg/float(NBKG)


#  ps = positive_sig/float(NSIG)

#  if oneMinusPb == 0:
#    CLs = 1.0
#    dCLs = 1.0
  
#  else: 
#    CLs = ps/oneMinusPb


# del n_pseu_sig
# del nuTil_b_pseu_sig
# del nuTil_s_pseu_sig

  #result = [qMu_tot_obs, ps, qMu_tot_sig, oneMinusPb, qMu_tot_bkg, CLs]
#  result = [qMu_tot_obs, 0, 0, 0,0,0]
#  result_save_2 = result[5]
#  result_save_1 = result[0]
#  print result_save_2

#  with open ("cls_value.txt", "w") as f:
#    f.write("%f %f" % (result_save_1, result_save_2))



#  return (qMu_tot_obs, ps, qMu_tot_sig, oneMinusPb, qMu_tot_bkg, CLs)
  return (qMu_tot_obs, -L1_obs,L_max[0]["mu"],0,0,mu)


def cls_exp(b_in, db_in, o_in, s_in, ds_in, randomseed):
    # Attempt to import scipy with minimize. If this does not exist, CLs cannot be evaluated
    try:
        global minimize, brentq, fsolve
        from scipy.optimize import minimize, brentq, fsolve
    except ImportError:
        print "!!! IMPORTANT WARNING: CLS CANNOT BE EVALUATED SINCE SCIPY.MINIMIZE CANNOT BE LOADED! RETURNING -1"
        return (-1, -1)

    if randomseed != 0:
        np.random.seed = randomseed
        
#    global b0, s0, db, ds # Passes this to nexp function
#    b0 = b_in
#    s0 = s_in
#    db = db_in
#    ds = ds_in
    # find best fit nuisance parameters and corresponding nexp
    
    nu_b_hat_obs_0, nu_s_hat_obs_0 = nuWithMaxL(o_in, 0., 0., 0., b_in, db_in, s_in, ds_in)
    nexp_bkg = nexp(0., nu_b_hat_obs_0, nu_s_hat_obs_0, b_in, db_in, s_in, ds_in)
    # calculate CLs with these
    return calc_CLs([nexp_bkg], nu_b_hat_obs_0, [b_in], [db_in], nu_s_hat_obs_0, [s_in], [ds_in], 10000, 1,  1., 1., 1.)
#b=[109.9000015258789, 1.5, 4.409999847412109, 0.9800000190734863, 8.5]
#db=[11.399999618530273, 0.4000000059604645, 1.1399999856948853, 0.5, 3.5999999046325684]
#s=[7.478017807006836, 1.4449514150619507, 235.7992706298828, 185.80776977539062, 71.80697631835938]
#ds=[0.49771007895469666, 0.26363080739974976, 6.895144462585449, 6.115604877471924, 4.633857727050781]
#n=[141.0, 3.0, 2.0, 0.0001, 8.0]


b=[ 0.9800000190734863]
db=[0.5]
s=[235.7992706298828]
ds=[6.895144462585449]
n=[0.0001]



test =  calc_CLs(n, 0.0,b , db, 0.0, s, ds, 0, 1, 1, 1, 1)
#print test[0]
