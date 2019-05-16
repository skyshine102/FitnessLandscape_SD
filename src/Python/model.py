# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 19:52:43 2018

@author: JeremyJ
"""
import numpy as np
import util

import matplotlib.pylab as plt 
from sklearn.metrics import mean_squared_error, r2_score

class Model(object):   
    def __init__(self, name,rep, dic, epi_dic):
        self.name = name
        self.replicate = rep
        self.dic = dic
        self.epi_dic = epi_dic
        self.fitted_expression_list = None
        self.explained_variance_ratio_ = None 
        
    def variance_decomposition_fit(self, order_list = list(range(1,10)), save = False):
        '''
        create "a dictionary of list" with fitted values of order i = 1 ~ 9
        '''
        if self.fitted_expression_list == None:
            print('Fitting values... ')
            predictive_values = {}   
            background = np.nanmean(list(self.dic.values()))
            
            
            for index, order in enumerate(order_list):
                array = np.array([ np.nansum([self.epi_dic[order-1][k] for k in util.SubSeq(i,9-order)]  ) for i in self.dic.keys()])
                predictive_values[index] = array
            print(predictive_values)
            
            predictive_values[0] = predictive_values[0] + background
            for i in range(1,len(order_list)):
                predictive_values[i] = predictive_values[i] + predictive_values[i-1]
            
            
            self.fitted_expression_list = predictive_values
        
        # save temporary npy file for fast access and loading predicted data
        np.save(self.name + "_" + self.replicate + "_predicted_epistasis.npy", self.fitted_expression_list)
        
        if self.explained_variance_ratio_ == None:
            self.explained_variance_ratio_ = np.zeros(len(order_list))

        for order in order_list:
            x = np.array(list(self.dic.values()))
            y = self.fitted_expression_list[order-1]
            ind_x = ~np.isnan(x)
            ind_y = ~np.isnan(y)
            ind = np.logical_and(ind_x, ind_y)
            x = x[ind]
            y = y[ind]
            plt.scatter(x,y, s= 0.05)
            if save == True:
                plt.savefig("scatterplot_fitted-exp"+self.name+"order"+str(order)+".svg")
            plt.show()
            r2 = r2_score(x,y)
            print("r-squared:", r2)
            self.explained_variance_ratio_[order-1] = r2
    
    
    def readIn_predicted_npy_files(self):
        self.fitted_expression_list = np.load(self.name + "_" + self.replicate + "_predicted_epistasis.npy").item() 
        self.explained_variance_ratio_ = np.zeros(9)
        for order in list(range(1,10)):            
            x = np.array(list(self.dic.values()))
            y = self.fitted_expression_list[order-1]
            ind_x = ~np.isnan(x)
            ind_y = ~np.isnan(y)
            ind = np.logical_and(ind_x, ind_y)
            x = x[ind]
            y = y[ind]
            r2 = r2_score(x,y)
            print("r-squared:", r2)
            self.explained_variance_ratio_[order-1] = r2
            '''
            plt.scatter(x,y, s= 0.05)
            plt.savefig("scatterplot_fitted-exp"+self.name+"order"+str(order)+".svg")
            plt.show()
            '''
    
    
    def regression_fit(self):
        pass