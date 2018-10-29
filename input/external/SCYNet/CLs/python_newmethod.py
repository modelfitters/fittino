import sys
import math
import likelihood_sysstat as CM
import likelihood_sysstat_2 as CM2
import likelihood_sysstat_3 as CM3
import likelihood_sysstat_1 as CN
import numpy as np

def main(obs,bkg, bkg_err,si, si_err_stat, si_err_sys ):
    
    nr_SR=4
#    SR1=[0]
#    SR2=[0 for a in range(nr_poi)]
#    SR3=[0 for a in range(nr_poi)]
#    SR4=[0 for a in range(nr_poi)]

    result_vector = [[]]
    result_vector_2 = [[]]
    count = 0
    count1 = 0
    if si[0] <=0.00001 : 
     print "kein Signal"
     print si[i]
     SR1=[0,100,100,100,100, 0]  
    elif si[1] <=0.00001:
     if si_err[0]>= 3*si[0]: si_err[0] = 3*si[0]
     nr_SR =1
     SR1= CN.calc_CLs([obs[0]], 0.0, [bkg[0]], [bkg_err[0]] , 0.0, [si[0]], [si_err[0]], count, nr_SR, 1., 1., 1.)

    elif si[2] <=0.0001:
     if si_err[0]>= 3*si[0]: si_err[0] = 3*si[0]
     if si_err[1]>= 3*si[1]: si_err[1] = 3*si[1]
     if si_err_stat[0]>= 3*si[0]: si_err_stat[0] = 3*si[0]
     if si_err_stat[1]>= 3*si[1]: si_err_stat[1] = 3*si[1]
     if si_err_sys[0]>= 3*si[0]: si_err_sys[0] = 3*si[0]
     if si_err_sys[1]>= 3*si[1]: si_err_sys[1] = 3*si[1]
     nr_SR =2
     SR1= CM2.calc_CLs([obs[0],obs[1]], 0.0, [bkg[0], bkg[1]], [bkg_err[0], bkg_err[1]] , 0.0, [si[0], si[1]], [si_err_stat[0], si_err_stat[1]], [si_err_sys[0], si_err_sys[1]], count, nr_SR)
    elif si[3] <=0.0001:
     if si_err[0]>= 3*si[0]: si_err[0] = 3*si[0]
     if si_err[1]>= 3*si[1]: si_err[1] = 3*si[1]
     if si_err[2]>= 3*si[2]: si_err[2] = 3*si[2]
     if si_err_stat[0]>= 3*si[0]: si_err_stat[0] = si[0]
     if si_err_stat[1]>= 3*si[1]: si_err_stat[1] = si[1]
     if si_err_stat[2]>= 3*si[2]: si_err_stat[2] = si[2]
     if si_err_sys[0]>= 3*si[0]: si_err_sys[0] = si[0]
     if si_err_sys[1]>= 3*si[1]: si_err_sys[1] = si[1]
     if si_err_sys[2]>= 3*si[2]: si_err_sys[2] = si[2]
     nr_SR =3
     SR1= CM3.calc_CLs([obs[0],obs[1],obs[2]], 0.0, [bkg[0], bkg[1], bkg[2]], [bkg_err[0], bkg_err[1], bkg_err[2]] , 0.0, [si[0], si[1], si[2]], [si_err_stat[0], si_err_stat[1], si_err_stat[2]], [si_err_sys[0], si_err_sys[1], si_err_sys[2]], count, nr_SR)
    else:
     if si_err_stat[0]>= 3*si[0]: si_err_stat[0] = si[0]
     if si_err_stat[1]>= 3*si[1]: si_err_stat[1] = si[1]
     if si_err_stat[2]>= 3*si[2]: si_err_stat[2] = si[2]
     if si_err_stat[3]>= 3*si[3]: si_err_stat[3] = si[3]
     if si_err_sys[0]>= 3*si[0]: si_err_sys[0] = si[0]
     if si_err_sys[1]>= 3*si[1]: si_err_sys[1] = si[1]
     if si_err_sys[2]>= 3*si[2]: si_err_sys[2] = si[2]
     if si_err_sys[3]>= 3*si[3]: si_err_sys[3] = si[3]
     nr_SR =4
     SR1= CM.calc_CLs([obs[0],obs[1],obs[2],obs[3]], 0.0, [bkg[0], bkg[1], bkg[2], bkg[3]], [bkg_err[0], bkg_err[1], bkg_err[2], bkg_err[3]] , 0.0, [si[0], si[1], si[2], si[3]], [si_err_stat[0], si_err_stat[1], si_err_stat[2], si_err_stat[3]], [si_err_sys[0], si_err_sys[1], si_err_sys[2], si_err_sys[3]], count, nr_SR) 


#    for i in range(start, nr_poi):
#      print "-----------------------"
#      if preselection_fail[i] == False:
    SR2=SR1[0]
    SR3=SR1[1]
    SR4=SR1[2]
#   result = [qMu_tot_obs, ps, qMu_tot_sig, oneMinusPb, qMu_tot_bkg, CLs]
#   result_save_2 = result[5]
#   result_save_1 = result[0]
#   print result_save_2
    print SR2
    print SR3
    print SR4
    with open ("cls_value.txt", "w") as f:
     f.write("%f %f %f" % (SR2, SR3, SR4))

      
#main(0,10,0)
