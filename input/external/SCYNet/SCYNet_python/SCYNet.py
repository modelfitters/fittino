from copy import deepcopy
from Preprocessor_lena import pmssm, chi2, log_norm, div_max, min_max, sub_mean_div_std
from keras.models import model_from_json, load_model
#import misc
import numpy as np
import sys
import time
import datetime
from keras.utils.generic_utils import get_custom_objects
import keras.backend as K
import os


def load_result_file(path):
    data = []
    with open(path,'r') as file:
        lines = file.readlines()
    for l in lines:
        data.append(result_string_to_dict(l,verbose=True))
    return data




def result_string_to_dict(line, verbose=True):
    hp = {}
    L = line.split(';')
    for l in L:
        key, value = l.split('?')
        if verbose: print (key, '\t', value)
        hp.update({str(key):str(value)})
        try:
            print hp[key]
            hp[key] = eval(hp[key])
        except NameError:
            pass
    return hp

class SCYNet:
    '''generate a fast chi2 prediction for pmssm points'''
    def __init__(self, args, inp ,model='SCYNET.h5', hp='hyperpoint.txt'):
              
        self.hp = load_result_file(hp)[0] #hyperparameter+information on the preprocessing. we need this to apply to the input pmssm the user gives.

        for k,v in zip(self.hp.keys(), self.hp.values()):
          print (k,v)
        for attribute in self.hp['chi2_back_info']:
            value = self.hp['chi2_back_info'][attribute]
            if attribute == 'std': y_dev = value
            if attribute == 'mean': y_mean = value
            print ('setting chi2.%s to %s'%(attribute,value))
        def abspl(y_true, y_pred):
          diff=(K.square(K.abs((y_true-y_pred)*y_dev/K.clip(K.abs((y_true-4.0)*y_dev+y_mean),K.epsilon(), None))))
#  diff=(K.abs(y_true-y_pred)*y_dev)
          return 100. * K.mean(diff,axis=-1) 
        get_custom_objects().update({"abspl": abspl})
        self.inp = inp 
        self.model = load_model(model) #generate keras model
        self.hp = load_result_file(hp)[0] #hyperparameter+information on the preprocessing. we need this to apply to the input pmssm the user gives.
        for i in range(len(args)):
            if args[i] =='-data' or args[i] == '-d':
                self.datapath = args[i+1]
            if args[i] == '-m' or args[i] == '-mask':
                #for example, masks = [[60,90],[91,92]] will output only chi2 in this ranges
                self.masks = [map(float, arg.split('-')) for arg in args[i+1].split(',')]
            if args[i] == '-o' or args[i] == '-output':
                self.outputfile = str(args[i+1])
            if args[i] == '-v' or args[i] == '-verbose':
                self.verbose=True
            if args[i] == '-i' or args[i] == '-include':
                try:
                    self.include_additionally = np.load(args[i+1])
                except:
                    print ('need valid numpy file, could not load your desired data')
            if args[i] == '-w' or args[i] == '-write':
                self.mode = args[i+1]
        if not hasattr(self, 'outputfile'):
            self.outputfile = 'SCYNet.out'
        if not hasattr(self, 'masks'):
            self.masks = [[0.,20.]]
        if not hasattr(self, 'verbose'):
            self.verbose = False
        if not hasattr(self, 'mode'):
            self.mode = 'all'
        if hasattr(self, 'datapath'):
            self.read_data(self.datapath)
        if not hasattr(self, 'datapath'):
            self.read_data(self.inp)
        print ('initialized SCYNet, your options:')
        print ('\tmasks:', self.masks)
        print ('\toutputfile:', self.outputfile)
        if hasattr(self, 'datapath'):
            print ('\tinput data:', self.datapath)
        print ('\tverbose:', self.verbose)

    def predict(self, data=None):
        if data == None:
            if hasattr(self, 'pp_data'):
                self.pred = self.model.predict(self.pp_data)
            else:
                raise Exception('SCYNet has no attribute "pp_data" expected at this point')
                data = self.inp
                print data
                data = np.array(data)
                print data[4]
                print data.shape[0]
                self.pred = self.model.predict(data)
                
        elif isinstance(data, np.ndarray) and data.shape[1] == 11:
            print "Hallo2"
            self.data = data
            self.preprocess() #preprocesses the just arrived self.data to self.pp_data
            self.pred = self.model.predict(self.pp_data)
        else:
            raise Exception('SCYNet has no pmssm data with shape(N,11) to predict chi2')
        self.pred = self.pred.flatten()
        return self.backtransform_chi2() #need to backtransform chi2 into 0..100 range

    def preprocess(self):
        '''preprocesses the p11mssm points in exactly the same way done for the net as it was trained'''
        if self.verbose:
            for k,v in zip(self.hp.keys(), self.hp.values()):
                print (k,v)
        self.pp_data = deepcopy(self.data)
        for func in self.hp['pp_pmssm']:
            print func
            if self.verbose: print ('\n*****func:',func)
            preproc_info = self.hp['pmssm_back_info'][func]
            if self.verbose: print ('\n*****preproc_info', preproc_info)
            self.pp_data = eval(func)(self.pp_data, preproc_info) #these functions got imported from Preprocessor above
       
    def backtransform_chi2(self):
        '''do the same transformations that were applied to the chi2 during the training backwards'''
        params=[self.hp['cut'], self.hp['delta']]

        #chi2 class has inbuilt functionality for the backtrafo.
        #It matters only that we init chi2 class with params,rest irrelevant
        y = chi2(data=np.array([]), preproc=[],params=params,split=0,verbose=True)
        y.back = self.hp['chi2_back_trafos']
        if self.verbose: print ('y.back:',y.back)

        #the chi2 class needs accesses the attributes in this way.
        for attribute in self.hp['chi2_back_info']:
            value = self.hp['chi2_back_info'][attribute]
            if self.verbose: print ('setting chi2.%s to %s'%(attribute,value))
            setattr(y, attribute, value)

        #actual backtransformation
        self.pred = y.backtrafo(self.pred) #applys functions in backward order
        return self.pred

   
    def read_data(self, inputfile):
        '''reading in the data. test if numpy file ('something.npy') or txtfile'''
        print "Hier" 
        self.data = np.array(self.inp)
        self.preprocess()


    def write_output(self,mode='all'):
        '''writes output to file
        '''
        self.mode = mode
        if self.verbose: print ('writing output ... mode: %s' % self.mode)
        output = ''
        written = 0
        self.masks_ = []
        self.total_mask = np.zeros(len(self.pred))
        for m in self.masks:
            chi_squared = deepcopy(self.pred)
            mask = (chi_squared >= m[0]) == (chi_squared < m[1])
            chi_squared = chi_squared[mask]
            self.total_mask += mask.astype(int)
            self.masks_.append([mask,chi_squared])
            print ('points in this chi2 range:', sum(mask.astype(int)))
            if self.mode == 'chi2_only':
                for i in range(len(chi_squared)):
                    output += str(chi_squared[i]) + '\n'
            elif self.mode == 'all':
                for i in range(len(chi_squared)):
                    line=' '.join(map(str, self.data[mask][i]))+' '+str(chi_squared[i])+'\n'
                    output += line
            elif self.mode == 'pmssm_only':
                for i in range(len(chi_squared)):
                    line=' '.join(map(str, self.data[mask][i]))+'\n'
                    output += line
            else:
                raise IOError('wrong mode, can only have "chi2_only","pmssm_only", or "all"')
            written += len(chi_squared)
            if self.verbose: print ('will write %s lines in mask %s' % (len(chi_squared), m))
        
        with open(self.outputfile, 'w') as file:
            file.write(output)
            if self.verbose: print ('wrote %s lines' % written)
def main(a):

    net = '1.1565_03Aug'
    print (net+'.txt')
    SN = SCYNet(args= sys.argv, inp = a, model=net+'.h5', hp=net+'.txt')
    #SN = SCYNet(args='/home/lb948295/SCY_Net/data/points_test.npy', model=net+'.h5', hp=net+'.txt')
    pred = SN.predict()
    print pred
    SN.write_output(mode = 'all')

if __name__ == '__main__':
    #chosing SCYNet!


#    net = '1.1110_21Jun'
#    print (net+'.txt')
#    SN = SCYNet(args= sys.argv, inp = [[-2568.98183931,   703.93268716,   -1166.66417905,   723.10329132,  253.2964175,  1553.93306017, 2575.6808874, 3404.42270274, -1787.54625896, 3262.60239357, 19.64002462]], model=net+'.h5', hp=net+'.txt')
    #SN = SCYNet(args='/home/lb948295/SCY_Net/data/points_test.npy', model=net+'.h5', hp=net+'.txt')
#    pred = SN.predict()
#    print pred
#    SN.write_output(mode = 'all')
    main([[-2568.98183931,   703.93268716,   -1166.66417905,   723.10329132,  253.2964175,  1553.93306017, 2575.6808874, 3404.42270274, -1787.54625896, 3262.60239357, 19.64002462]])
     
