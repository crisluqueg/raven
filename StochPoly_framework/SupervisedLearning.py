'''
Created on Mar 16, 2013

@author: crisr
'''
from __future__ import division, print_function, unicode_literals, absolute_import
import warnings
warnings.simplefilter('default',DeprecationWarning)

from sklearn import svm
import numpy as np
import Datas
import numpy
import h5py
from itertools import product as itprod
try:
  import cPickle as pk
except:
  import pickle as pk

#import DataBases #TODO shouldn't need this, see StochPoly.train() for instance check
'''here we intend ROM as super-visioned learning, 
   where we try to understand the underlying model by a set of labeled sample
   a sample is composed by (feature,label) that is easy translated in (input,output)
   '''
'''module to be loaded from scikitlearn
 Generalized Linear Models
 Support Vector Machines
 Stochastic Gradient Descent
 Nearest Neighbors
 Gaussian Processes
 Partial Least Squares
 Naive Bayes
 Decision Trees
 Ensemble methods
 Multiclass and multilabel algorithms
 Feature selection
 Linear and Quadratic Discriminant Analysis
 Isotonic regression
 '''

class superVisioned():
  def __init__(self,**kwargs):
    self.initializzationOptionDict = kwargs
  
  def fillDist(self,distributions):
    self.distDict = distributions
  
  def train(self,obj):
    '''override this method to train the ROM'''
    return

  def reset(self):
    '''override this method to re-instance the ROM'''
    return

  def evaluate(self):
    '''override this method to get the prediction from the ROM'''
    return

  def returnInitialParamters(self):
    '''override this method to pass the fix set of parameters of the ROM'''
    InitialParamtersDict={}
    return InitialParamtersDict

  def returnCurrentSetting(self):
    '''override this method to pass the set of parameters of the ROM that can change during simulation'''
    CurrentSettingDict={}
    return CurrentSettingDict

class StochasticPolynomials(superVisioned):
  def __init__(self,**kwargs):
    superVisioned.__init__(self,**kwargs)
    self.poly_coeffs={}

  def train(self,inDictionary):
    data=inDictionary['Input'][0]
    print('\n\n\ndata:',data,'\n\n\n')
    self.solns={}
    
    if data.type=='HDF5':
      attr={}
      hists=data.getEndingGroupNames()
      M=[]
      for i,h in enumerate(hists):
        if h=='':continue
        attr['history']=h
        M.append(data.returnHistory(attr))
      
      if data.targetParam:
        self.targetParam = data.targetParam
      else:
        raise IOError('No target Parameter for ROM Stochastic Polynomials')
      # How to get specific values from solution?
#TODO this should be an input on the front end; the user should choose the index
      self.solnIndex=numpy.where(M[0][1]['headers']==self.targetParam)


      # for each run, sampler passes the values (quad pt) to eval at, as well as
      # the partial coefficients for that _quad point_.  Here, for each of those,
      # we simply need to sum over each partCoeff[quad_pt][ord]*soln[quad_pt]
      # to construct poly_coeff[ord]

      
      for history in M:
        self.poly_coeffs[tuple(history[1]['exp_order'])]=0
        for partCoeff in history[1]['partial_coeffs']:
          self.poly_coeffs[tuple(history[1]['exp_order'])]+=\
                    history[0][0][self.solnIndex]*partCoeff
      
      print('StochasticPolynomials ROM successfully trained.')
    else:
      print('Reading from non-HDF5 for StochPolys not supported yet...')
    return

  def evaluate(self,data):
    # valDict is dict of values to evaluate at, keyed on var
    #FIXME these need to be adjusted for changes in train()
    tot=0
    if type(data) == 'dict':
      valDict = data
    else:
      valDict = data.getInpParametersValues()
      # other temporary fix -.-
      for k in valDict.keys():
        for ke in self.distDict.keys():
          if k in ke:
            temp = valDict.pop(k)
            valDict[ke] = temp
      
    
    for ords,coeff in self.poly_coeffs:
#      tot+=coeff*np.prod([self.distDict[var].quad().evNormPoly(\
#              ords[v],self.distDict[var].revertPt(valDict[var])) for v,var in enumerate(valDict)])
      tot+=coeff*np.prod([self.distDict[var].quad().evNormPoly(\
              ords,self.distDict[var].revertPt(valDict[var])) for v,var in enumerate(valDict)])
      #TODO revertPt may not always be straightforward to implement!
    return tot

  def reset(self,*args):
    pass




class SVMsciKitLearn(superVisioned):
  def __init__(self,**kwargs):
    superVisioned.__init__(self,**kwargs)
    #dictionary containing the set of available Support Vector Machine by scikitlearn
    self.availSVM = {}
    self.availSVM['LinearSVC'] = svm.LinearSVC
    self.availSVM['C-SVC'    ] = svm.SVC
    self.availSVM['NuSVC'    ] = svm.NuSVC
    self.availSVM['epsSVR'   ] = svm.SVR
    if not self.initializzationOptionDict['SVMtype'] in self.availSVM.keys():
      raise IOError ('not known support vector machine type ' + self.initializzationOptionDict['SVMtype'])
    self.SVM = self.availSVM[self.initializzationOptionDict['SVMtype']]()
    kwargs.pop('SVMtype')
    self.SVM.set_params(**kwargs)
    return

  def train(self,data):
    ''' The data is always a dictionary'''
    """Fit the model according to the given training data.
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vector, where n_samples in the number of samples and
            n_features is the number of features.
        y : array-like, shape = [n_samples]
            Target vector relative to X
        class_weight : {dict, 'auto'}, optional
            Weights associated with classes. If not given, all classes
            are supposed to have weight one.
        Returns
        -------
        self : object
            Returns self.
        fit( X, y, sample_weight=None):"""
    if(data.type != 'TimePointSet'):
      raise IOError('The SVM type ' + self.initializzationOptionDict['SVMtype'] + 'requires a TimePointSet to be trained')
    self.trainInputs = data.getInpParametersValues().items()
    self.trainTarget = data.getOutParametersValues().items()
    X = np.zeros(shape=(self.trainInputs[0][1].size,len(self.trainInputs)))
    y = np.zeros(shape=(self.trainTarget[0][1].size))
    for i in range(len(self.trainInputs)):
      X[:,i] = self.trainInputs[i][1]
    y = self.trainTarget[0][1]

    print('SVM           : Training ' + self.initializzationOptionDict['SVMtype'])
    self.SVM.fit(X,y)
    print('SVM           : '+ self.initializzationOptionDict['SVMtype'] + ' trained!')
    
    return

  def returnInitialParamters(self):
    return self.SVM.get_params()

  def evaluate(self,data):
    """Perform regression on samples in X.
        For an one-class model, +1 or -1 is returned.
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
        Returns
        -------
        y_pred : array, shape = [n_samples]
        predict(self, X)"""
    if(data.type != 'TimePointSet'):
      raise IOError('The SVM type ' + self.initializzationOptionDict['SVMtype'] + 'requires a TimePointSet to be trained')
    trainInputs = data.getInpParametersValues().items()
    X = np.zeros(shape=(trainInputs[0][1].size,len(trainInputs)))
    for i in range(len(trainInputs)):
      X[:,i] = self.trainInputs[i][1]
    print('SVM           : Predicting by ' + self.initializzationOptionDict['SVMtype'])
    return self.SVM.predict(X)

  def reset(self):
    self.SVM = self.availSVM[self.initializzationOptionDict['SVMtype']](self.initializzationOptionDict)

def returnInstance(Type):
  '''This function return an instance of the request model type'''
  base = 'superVisioned'
  InterfaceDict = {}
  InterfaceDict['SVMscikitLearn'       ] = SVMsciKitLearn
  InterfaceDict['StochasticPolynomials'] = StochasticPolynomials
  try: return InterfaceDict[Type]
  except: raise NameError('not known '+base+' type '+Type)
  
