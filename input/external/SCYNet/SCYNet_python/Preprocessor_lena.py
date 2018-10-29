import os
import numpy as np
from copy import deepcopy

from sklearn.utils import shuffle
def shuffle_data(x, y):
    x.x, y.chi2 = shuffle(x.x, y.chi2) #option: random_state=0 controlling seed
    y.train_test_init_sq() #so we compare prediction with correct ytrain,test
    x.split_train_data() #new train split
    y.split_train_data() #new test split

class fulldata:
    check = False
    def __init__(self, path = os.environ['HOME'] + '/13TeV_chi2_disjoint_2.npy',use=range(2,14)):
        self.path = path
        #self.data = np.genfromtxt(self.path[:-4])[:,1:]
        self.full = np.load(self.path)
        self.data = self.full[:,use]
        #print 'loaded pmssm11 data with shape', self.data.shape
        if self.check==True:
            self.quality_ensurance()

    def shuffle(self,seed=None):
        if seed == None:
           self.randomseed = np.random.randint(0,65536**2-1)
        if seed != None:
           self.randomseed = seed

        np.random.seed(self.randomseed)
        np.random.shuffle(self.data)

        np.random.seed(self.randomseed)
        np.random.shuffle(self.full)
        
    def quality_ensurance(self):
        '''does not work properly the way i intend it...
        checks for nans'''
        print ('checking if a "nan" is in the data....')
        if np.isnan(sum(sum(self.data))):
            raise ValueError('found a "nan" in the input data')
        else:
            self.check = False

class pmssm:
    def __init__(self, data, preproc, split, use_only=range(11)):
        self.x = data
        self.use_only = range(self.x.shape[1])
        self.split = int(split*len(self.x))
        self.preproc = preproc
        self.back_info = {}
        self.preprocess() 

    def sub_mean_div_std(self):
        self.mod = np.zeros((self.x.shape[1],2))
        backinfo = []
        for i in range(self.x.shape[1]):
            mean = self.x[:,i].mean()
            std = self.x[:,i].std()
            self.x[:,i] = (self.x[:,i] - mean) / std
            self.mod[i] = [mean, std]
            backinfo.append([mean, std])
        self.back_info.update({'sub_mean_div_std':backinfo})

    def div_max(self):
        '''divide by maximum'''
        self.maxis = np.zeros(self.x.shape[1])
        for i in range(self.x.shape[1]):
            maxi = max(abs(self.x[:,i]))
            self.x[:,i] /= maxi 
            self.maxis[i] = maxi
        self.back_info.update({'div_max':list(self.maxis)})

    def min_max(self):
        '''normalize between zero and one'''
        self.minmaxis = np.zeros((self.x.shape[1],2))
        mimaxis = []
        for i in range(self.x.shape[1]):
            mini = min(self.x[:,i])
            maxi = max(self.x[:,i])
            self.minmaxis[i] = [mini, maxi]
            mimaxis.append([mini,maxi])
            self.x[:,i] = (self.x[:,i] - mini) / (maxi - mini)
        self.back_info.update({'min_max':mimaxis})        

    def log_norm(self):
        '''scale by log'''
        self.logis = np.zeros(self.x.shape[1])
        for i in range(self.x.shape[1]):
            p = self.x[:,i]
            mini = np.min(np.abs(p))
            mask = p < 0
            p[mask] *= -1
            p = np.log(p - mini + 2) #+2 regularizes, we want no log(0)
            self.logis[i] = mini
            #print 'logis', i, self.logis[i], mini, p[:10]
            p[mask] *= -1
            self.x[:,i] = p
        self.back_info.update({'log_norm':list(self.logis)})

    def abs_norm(self):
        '''take the absolute value of everything'''
        for i in range(self.x.shape[1]):
            self.x[:,i] = np.abs(self.x[:,i])
        self.back_info.update({'abs_norm':1})

    def split_train_data(self):
        '''simple split function'''
        use_only = self.use_only #never does anything. If you want to use only certain pmssm parameters, you may take advantage of this.
        self.train = self.x[:,use_only][:self.split]
        self.test = self.x[:,use_only][self.split:]

    def preprocess(self):
        self.test_init = self.x[:,self.use_only][self.split:]
        self.train_init = self.x[:,self.use_only][:self.split]
        for p in self.preproc:
            getattr(self, p)() #retrieves function name, calls it with ()
        self.split_train_data()

class chi2:
    '''reads in the chi2, has functionality for evaluation of model'''
    def __init__(self, data, preproc, params, split=6.0/8, verbose=False):
        self.verbose = verbose #just for some printing
        self.preproc = preproc
        self.params = params
        self.chi2 = data
        self.split = int(split*len(self.chi2))
        self.start = deepcopy(self.chi2[:50])
        self.back = [] #training+SCYNet
        self.back_info = {}#not needed during training. but for SCYNet.
        self.preprocess() 

    def split_train_data(self):
        self.train = self.chi2[:self.split]
        self.test = self.chi2[self.split:]
 
    '''divide chi2 by maximum'''
    def div_max(self):
        self.maximum = 1.0*max(self.chi2)
        self.chi2 /= self.maximum
        self.back.append('mult_max')#training + SCYNet
        self.back_info.update({'maximum':self.maximum})#used only for SCYNet
    def mult_max(self, chi2):
        return self.maximum * chi2

    '''substract mean and divide by std'''
    def sub_mean_div_std(self):
        self.mean = self.chi2.mean()
        self.std = self.chi2.std()
        self.chi2 = (self.chi2 - self.mean)/self.std
        self.back.append('mult_std_add_mean')
        self.back_info.update({'mean':self.mean, 'std':self.std})
    def mult_std_add_mean(self, chi2):
        print self.std
        print self.mean
        return self.std*chi2 + self.mean

    '''transform arbitrary thing to 0..1'''
    def min_max(self):
        self.mini = min(self.chi2)
        self.maxi = max(self.chi2)
        self.chi2 = (self.chi2 - self.mini)/(self.maxi - self.mini)
        self.back.append('max_min')
        self.back_info.update({'mini':self.mini, 'maxi':self.maxi})
    def max_min(self, chi2):
        return (self.maxi - self.mini)*chi2 + self.mini

    def hard_cut(self):
        cut = self.params[0]
        self.chi2[self.chi2 > cut] = cut

    '''smoothly cut chi2 data at self.cut'''
    def square_cut(self):
        cut, delta = self.params
        y = self.chi2
        self.mask = (y > (cut - delta)) == (y < (cut + delta)) #true for max-cut <= chi2 <= max+cut, false else
        y[self.mask] = cut - (y[self.mask] - cut - delta)**2/(4*delta)
        y[y > (cut + delta)] = cut
        self.chi2 = y   
        self.back.append('back_square')
        self.back_info.update({'cut':cut, 'delta':delta})
    def back_square(self, chi2):
        cut, delta = self.params
        self.cut = cut
        self.delta = delta
        y = chi2
        mask = (y > cut-delta) == (y < cut) #backtrafo needs a different mask than trafo
        y[mask]= (cut+delta)-2*np.sqrt(delta*(cut-y[mask]))
        y[y>cut] = cut
        return y

    
    def train_test_init_sq(self):
        #only works for square cut! #if tanh, add another!#TODO
        self.cut = self.params[0]
        self.train_init = deepcopy(self.chi2)[:self.split].flatten()
        self.train_init[self.train_init > self.cut] = self.cut

        self.test_init= deepcopy(self.chi2)[self.split:].flatten()
        self.test_init[self.test_init > self.cut] = self.cut

    def preprocess(self):
        if True: #if I implement tanh cut, change this#TODO
            self.train_test_init_sq()
        for p in self.preproc:
            prefunc = getattr(self, p)
            prefunc()
        self.split_train_data()   

    def backtrafo(self, chi2):
        for func in reversed(self.back): #apply in reverse order
            backfunc = getattr(self, func)
            chi2 = backfunc(chi2)
        return chi2

 
###############################
# preprocessing for SCYNet.py #
###############################
def sub_mean_div_std(data, smds):
    for i in range(data.shape[1]):
        mean, std = smds[i]
        data[:,i] = (data[:,i] - mean)/std
    return data
def log_norm(data, logis):
    for i in range(data.shape[1]):
        p = data[:,i]
        mask = p < 0
        p[mask] *= -1
        p= np.log(p - logis[i] + 2)
        p[mask] *= -1
        data[:,i] = p
    return data
def div_max(data, maxis):
    for i in range(data.shape[1]):
        data[:,i] /= maxis[i]
    return data
def min_max(data, minmaxis):
    for i in range(data.shape[1]):
        mini, maxi = minmaxis[i]
        data[:,i] = (data[:,i] - mini)/(maxi - mini)
    return data


if __name__ == '__main__': 
    y = chi2(['square_cut','div_max'], [100,25], split=7.0/8)

    back = y.backtrafo(y.train)
    print (y.start)

    for a,b,c in zip(y.train[:50], y.train_init[:50], back[:50]):
        print (b,c,c==b)#==c
    #x=pmssm(preproc = ['log_norm','min_max'], split = 7.0/8)
    #y=pmssm(preproc = ['sub_mean_div_std','min_max'], split = 7.0/8)

    #import matplotlib
    #matplotlib.use("GTKAgg")
    #import matplotlib.pyplot as plt

    #print x.x[:,0]
    #print y.x[:,0]

    #plt.hist(y.x[:,0], bins=100)
    #plt.show()
    
