# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 22:01:28 2018

@author: JeremyJ
"""

import pandas as pd
import matplotlib.pylab as plt 
import numpy as np
import itertools 
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
import util
from model import Model
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from scipy import stats
from sklearn.metrics import r2_score
#import os


# global path parameters 
readcount_threshold = 25
readcount_u_path = "../../data/Read Count Table/UnionTable/"
readcount_path = "../../data/Read Count Table/"
epi_path = "../../data/Epistasis Table/"
vienna_path = "../../data/Vienna Package Prediction/"

class library(object):
    """ Represents a experimental libraries
    Args:
    Attributes:
    """
    def __init__(self,lib_name,rep_num,dataframe):
        self.name = lib_name
        self.replicate = rep_num
        self.dataframe = dataframe
        self.dic = {i:j for i,j in zip(self.dataframe["SeqID"], self.dataframe["LogMean"])}
        self.bucket_dict = None
        self.epistasis = [{} for i in range( int(np.log(len(self.dataframe["SeqID"]))/np.log(4)))]
        self.min = np.nanmin(self.dataframe["LogMean"])
        self.max = np.nanmax(self.dataframe["LogMean"])
        self.mean = np.nanmean(self.dataframe["LogMean"])
        if 'Rep1' in list(self.dataframe):
            self.std = np.mean(self.dataframe[["SeqID","Rep1","Rep2","Rep3"]].std(axis=1))
        
        if self.replicate[:3] == 'rep':
            self.num_detected = len(self.dataframe[self.dataframe['Counts']>0])
            self.bool_readcount_threshold = (self.dataframe.Counts > readcount_threshold).values
        else:
            self.num_detected = len(self.dataframe[self.dataframe['LogMean'] > 0])
            
        if self.name == 'dmsC':
            self.spacing = 4
        elif self.name == 'arti' or self.name == 'yfp':
            self.spacing = 5 
        elif self.name == 'fepB':
            self.spacing = 6
        else:
            self.spacing = 10
    def __str__(self):
        return "library = {},\n replicate = {},\n total_num_variants = {},\n detected # of variants = {}".format(self.name,self.replicate,len(self.dataframe), self.num_detected)
    
    __repr__ = __str__
    
    
    def model(self):
        self.model = Model(self.name, self.replicate, self.dic, self.epistasis)
        
    
    def get_common_loci_dict(self, start, end):
        '''
        "start" and "end" are included.
        
        To Do: use mutational_set to implement (complete)
        '''

        '''
        allseqs = [''.join for i in list(itertools.product('ACGU', repeat = length))]
        
        list(filter(lambda x: x[start:end+1] == i, allseqs ))
        '''
        length = end - start + 1
        allseqs = [''.join(i) for i in list(itertools.product('ACGU', repeat = length))]
        total_length = len(self.dataframe['SeqID'][0])
        template = list('NNNNNNNNN')
        template[0:start] = ['A'] * (start-0)
        template[end+1:] = ['A'] * (total_length-1 - end)
        template_set = list(util.mutational_set(''.join(template)))
        
        dictionary = { i: np.nanmean([ self.dic[j.replace('N'*length,i)]  for j in template_set ]) for i in allseqs}
        self.common_loci_dic = dictionary
        '''
        x = list(self.common_loci_dic.values())
        unique, counts = np.unique(x, return_counts=True)
        print(len(unique))
        print(counts)
        '''
        
    def readIn_epistasis_table(self,index = list(range(1,10))):
        '''
        Need to be modify the input.. currently read rep1 data
        '''
        if self.name == 'vienna':
            path_epistasis = vienna_path + 'Vienna' +'_SDR_order'
        else:
            path_epistasis = epi_path + self.name +'_SDR_'+ self.replicate +'_order'
        for i in index:
            try:
                print('Now Read in epistasis order {}'.format(i), end="", flush=True)
                tmp = pd.read_csv(path_epistasis + str(i) +'.csv')
                self.epistasis[i-1] = {i:j for i, j in zip(tmp["SeqID"],tmp["LogMean"])}
            except FileNotFoundError:
                print("File of order {} NOT FOUND!!".format(i))        
    
    def generate_motif_contribution_outlier(self, orderlist = [1,2,3,4,5] ,topk = 15):
        table = pd.DataFrame()
        for i in orderlist:
            dic = dict(zip(self.epistasis[i-1].values(),self.epistasis[i-1].keys()))
            tmp_v = sorted(list(self.epistasis[i-1].values()),reverse = True)
            tmp_v = np.append(tmp_v[:topk], tmp_v[-topk:])
            tmp_k = [dic[i] for i in tmp_v]
            table['Order {}'.format(i)] = tmp_k
            table['O{} effect'.format(i)] = tmp_v
        table.to_csv("../data/Result/outliers/" + self.name + '_' + self.replicate+'motif_outlier.csv',index = False)
        
    
    def bucket_classifier(self,mini = 0, maxi = 3.2, number_bins = 20):
        def comparator(input_value, bins = np.arange(mini, maxi+maxi/number_bins, (maxi-mini)/number_bins)):
            for i in range(len(bins)-1):
                if input_value > bins[i] and input_value <= bins[i+1]:
                    return i
        
        values = list(map(lambda x: comparator(self.dic[x]), self.dic.keys()))
        self.bucket_dict = dict(zip(self.dic.keys(), values))      
    
    def cross_bucket_similarity(self, mini=0, maxi =3.2, number_bins = 20):
        '''
        This function is moved to Julia language version.
        '''
        # create a matrix to store the values 
        matrix = np.zeros((number_bins,number_bins))
        
        # classify the variants into bins 
        self.bucket_classifier(mini = mini, maxi = maxi, number_bins = number_bins)
        variants = np.array(list(self.bucket_dict.keys()))
        mapping = np.array(list(self.bucket_dict.values()))
        
        # upper triangle
        for i,j in list(itertools.combinations(range(20), r = 2)):
            print('>>>>> Current progress: >>>>>')
            print(i,j)
            
            v = 0
            inputdata = [ variants[mapping ==i], variants[mapping ==j] ]
            '''
            if len(inputdata[0])*len(inputdata[1]) < 10e6:
                print(len(inputdata[0]),len(inputdata[1]))
                result = list(itertools.product(*inputdata))
                print(result[0])
                for k,l in result:
                    v += util.hamming_distance(k,l)
                matrix[i,j] = v/len(result)
            else:
            '''

            for seq1 in inputdata[0]:
                a = [ [seq1], inputdata[1] ]
                for s1,s2 in list(itertools.product(*a)):
                    v += util.hamming_distance(s1,s2)
            matrix[i,j] = v/len(inputdata[0])/len(inputdata[1])
                

        # diagonal
        for i in range(20):
            print('Current progress:')
            print(i)
            v = 0
            inputdata = [ variants[mapping ==i], variants[mapping ==i] ]
            
            '''
            if len(inputdata[0])*len(inputdata[1]) < 10e6:
                result = list(itertools.product(*inputdata))
                for k,l in result:
                    v += util.hamming_distance(k,l)
                matrix[i,i] = v/len(result)
            '''
            for seq1 in inputdata[0]:
                a = [ [seq1], inputdata[1] ]
                for s1,s2 in list(itertools.product(*a)):
                    v += util.hamming_distance(s1,s2)
            matrix[i,j] = v/len(inputdata[0])/len(inputdata[1])
                
        
        plt.imshow(matrix, cmp = plt.get_cmap('jet'))
        plt.show()

        # minimum spanning tree?



    def RBS_calculator(self, path = "../data/ChouLab_Libraries_1_4_RBSCalcv2/", save = False):
        file = pd.read_csv(path + "ChouLab_Library_simplified_" + self.name + "SD.csv" )
        tmp_dict = dict(zip(file['mRNA_Sequence'],file['Translation_Initiation_Rate']))
        emperical = [self.dic[i] for i in self.dic.keys()]
        RBS = [np.log10(tmp_dict[i]) for i in self.dic.keys()]
        print("Max in RBS = ", np.max(RBS))
        print("Mean in RBS = ", np.mean(RBS))
        print("Min in RBS = ", np.min(RBS))
        
        plt.figure(figsize=(6,6))
        plt.scatter(RBS,emperical, color = 'g', s = 1)
        r = util.correlation_coef(RBS,emperical)
        print("Pearson r is: ", r)
        plt.ylim(0,3.2)
        plt.yticks(np.arange(0,3.2+3.2/4, 3.2/4))
        plt.xlim(-1,7)
        plt.xticks(np.arange(-1,7+8/4, 8/4))
        #plt.xticks(np.arange(-x_max,x_max+x_max/2, x_max/2))    
        if save == True:
            plt.savefig("../data/Result/RBS_calculator/" + self.name + "yRBSxEmperical.tif")
        plt.show()
        return RBS
    
        #plt.plot()
        #print(file)
        
        
        
        





    
    
    def distribution(self, title = '', normed = False , c = 'g',log = False, save = False, mini = 0, maxi = 3.2, number_bins = 20, ticks = 5):
        '''
        A function used to visualize the distribution of the library
        '''
        a = self.dataframe["LogMean"]
        a = a[~np.isnan(a)]
        if normed == True:
            weights = np.ones_like(a)/float(len(a))
        else:
            weights = np.ones_like(a)

        bins = np.arange(mini, maxi+maxi/number_bins, (maxi-mini)/number_bins)
        print('bins =', bins)
        n, bins, patches = plt.hist(a, weights= weights, bins = bins, color = c, log= log)
        plt.title(title)
        plt.xticks(np.arange(mini, maxi + (maxi-0)/(ticks-1) ,(maxi-mini)/(ticks-1) ))
        
        if normed == True:
            plt.yticks(np.arange(0, 0.7,0.1))
        elif log == True:
            plt.yscale('log')
            plt.ylim(0,10**6)
            plt.yticks([1,10,100,1000,10000,100000])
        if save == True:
            plt.savefig("../data/Result/distribution/"+self.name + "_"+ self.replicate +"_distribution"+".svg")          
        plt.show()
        
        print("The min value =",np.nanmin(a))
        print("The max value =",np.nanmax(a))
        print("The mean value =",np.nanmean(a))
        
        print("Total variants =", np.count_nonzero(a))
        print("Coverage = ", "{:.2%}".format(np.count_nonzero(a)/4**9))
        return n, bins, patches
    
    
    def times_detect_distribution(self, num_bins = 20):
        
        if 'rep' in self.replicate:
            print("This function does not support replicate table!!")
        else:
            # missing ratio
            na1 = ~np.isnan(self.dataframe['Rep1'])
            print(np.count_nonzero(na1))
            na2 = ~np.isnan(self.dataframe['Rep2'])
            print(np.count_nonzero(na2))
            na3 = ~np.isnan(self.dataframe['Rep3'])
            print(np.count_nonzero(na3))
            na = np.logical_and(np.logical_and(na1,na2),na3)
            print(np.count_nonzero(na3))
            detected_vec = na.astype(int)
            #count_detected = np.ones(262144)*3 - (na1.astype(int) + na2.astype(int) + na3.astype(int))
            self.bucket_classifier()
            average_detected_bins = [[] for i in range(20)]
            for index, i in enumerate(self.dic.keys()):
                if self.bucket_dict[i]: # due to None type
                    average_detected_bins[ self.bucket_dict[i] ].append(detected_vec[index])  #count_detected[index] ) 
            average_detected_array = [ np.nanmean(j) for j in average_detected_bins]
            
            plt.bar(np.arange(1,21)-0.5 ,average_detected_array)
            plt.xlim(0,20)
            plt.xticks(np.arange(0,21))
            plt.ylim(0,1)
            plt.yticks(np.arange(0,1+0.1,0.1))
            plt.show()
        
    
    
    def nucleotide_trend(self,save = False, ymax = 3.2, yticks = 5):
        mean = np.zeros((4,10))
        std = np.zeros((4,10))
        # non_nan = ~np.isnan(vec_expression)
        
        vec_expression = np.array(list(self.dic.values()))
        vec_A = np.array(list(map(lambda x: x.count('A') , self.dic.keys())))
        vec_U = np.array(list(map(lambda x: x.count('U') , self.dic.keys())))
        vec_C = np.array(list(map(lambda x: x.count('C') , self.dic.keys())))
        vec_G = np.array(list(map(lambda x: x.count('G') , self.dic.keys())))
        for i in range(10):
            A = vec_expression[vec_A == i] 
            U = vec_expression[vec_U == i]
            C = vec_expression[vec_C == i]
            G = vec_expression[vec_G == i]
            mean[0][i] = np.nanmean(A)
            std[0][i] = np.nanstd(A)
            mean[1][i] = np.nanmean(U)
            std[1][i] = np.nanstd(U)
            mean[2][i] = np.nanmean(C)
            std[2][i] = np.nanstd(C)
            mean[3][i] = np.nanmean(G)
            std[3][i] = np.nanstd(G)
        ''' 
        plt.plot(mean[0], c = 'g') 
        plt.plot(mean[1], c = 'r')
        plt.plot(mean[2], c = 'b')
        plt.plot(mean[3], c = 'k')
        '''
        idx = np.arange(0,10)

        plt.errorbar(idx-0.2, mean[0], yerr = std[0], fmt='g-', linewidth = 3, ecolor='g',elinewidth = 2, capsize = 3)
        plt.errorbar(idx+0.1, mean[1], yerr = std[1], fmt='r-', linewidth = 3, ecolor='r',elinewidth = 2, capsize = 3)
        plt.errorbar(idx-0.1, mean[2], yerr = std[2], fmt='b-', linewidth = 3, ecolor='b',elinewidth = 2, capsize = 3)
        plt.errorbar(idx, mean[3], yerr = std[3], fmt='k-', linewidth = 3, ecolor='k',elinewidth = 2, capsize = 3)
        plt.ylim(0,ymax)
        plt.yticks(np.arange(0,ymax + (ymax-0)/(yticks-1), (ymax-0)/(yticks-1) ))
        plt.xticks(np.arange(0,10))
        if save == True:
            plt.savefig("../data/Result/nucleotide_trend/" + 'nucleotide_trend_'+ self.name + "_" +self.replicate + '.svg')
        plt.show()
        
        return mean, std
        '''
        result = stats.spearmanr(vec_expression[non_nan],vec_A[non_nan])
        print('A')      
        print("Spearman-r:", result[0])
        print("P- value:", result[1])
        result = stats.spearmanr(vec_expression[non_nan],vec_U[non_nan])
        print('U')
        print("Spearman-r:", result[0])
        print("P- value:", result[1])
        
        result = stats.spearmanr(vec_expression[non_nan],vec_G[non_nan])
        print('G')
        print("Spearman-r:", result[0])
        print("P- value:", result[1])
        result = stats.spearmanr(vec_expression[non_nan],vec_C[non_nan])
        print('C')
        print("Spearman-r:", result[0])
        print("P- value:", result[1])
        '''
        
        
    def nucleotide_composition(self, save = False, mini = 0 , maxi = 3.2 , number_bins = 20, ticks = 5):
        # set an array to store
        table = np.zeros((20,4))  # columns: A, U, C, G 
        b = np.arange(mini, maxi + (maxi-mini)/number_bins, (maxi-mini)/number_bins)
        print(b,len(b))
        for index, v in enumerate(b[:-1]):
            tmp = list(filter(lambda x: (self.dic[x]>=v and self.dic[x] < v+  (maxi-mini)/number_bins ), self.dic.keys()))
            longstring = ''.join(tmp)
            table[index][0] = longstring.count('A')
            table[index][1] = longstring.count('U')
            table[index][2] = longstring.count('C')
            table[index][3] = longstring.count('G')
            # Get the weight of each bin
            #print(table[index,:])
            if np.sum(table[index,:]) != 0:
                table[index,:] = table[index,:]/np.sum(table[index,:])*100
        # plot
        #from matplotlib import rc
        # y-axis in bold
        #rc('font', weight='bold')    
        N = 20
        ind = np.arange(N)    # the x locations for the groups
        width = 1       # the width of the bars: can also be len(x) sequence
        x_interval = (maxi-mini)/(ticks-1)
    
        plt.bar(ind+0.5, table[:,0], width, color = 'g', edgecolor='white')
        plt.bar(ind+0.5, table[:,1], width, bottom=table[:,0], color = 'r',  edgecolor='white')
        plt.bar(ind+0.5, table[:,2], width, bottom=(table[:,1]+table[:,0]), color ='b', edgecolor='white')
        plt.bar(ind+0.5, table[:,3], width, bottom=(table[:,2]+table[:,1]+table[:,0]), color = 'k', edgecolor='white')
        
        #plt.xticks(r, names, fontweight='bold')
        
        # This line is problematic
        tmp = [N/(ticks-1)*i for i in range(ticks)]
        plt.xticks(tmp,[str(round(i,2)) for i in np.arange(mini,maxi+x_interval,x_interval)])
        #plt.xticks(ind, ['R'+str(i) for i in range(1,21)])
        plt.yticks(np.arange(0, 110, 10))
        if save == True:
            plt.savefig( "../data/Result/nucleotide_composition/"+"nt_composition_"+self.name+"_"+self.replicate+".svg" )
        #plt.tight_layout()
        #ax = plt.gca().yaxis.grid(True)
        plt.title(self.name)
        plt.show()
        
        return table
    
    def read_count_distribution(self, save = False):
        if self.replicate[:3] == 'rep':
            bins = np.array([0,1,26,51,76,101,126,151,176,201,226,251,301,351,401,501,601,801,1000,500000])
            labels = ['0','1-25','26-50','51-75','76-100','101-125','126-150','151-175','176-200','201-225','226-250',
                      '251-300','301-350','351-400','401-500','501-600','601-800','801-1000','1000up'
                      ]
            print('bins =', bins)            
            n, bins, patches = plt.hist(self.dataframe['Counts'], bins= bins)
            plt.close()
            plt.bar(range(len(n)),n/np.sum(n)*100,tick_label=labels,  color = 'forestgreen')
            plt.xticks(rotation=90)
            plt.xlabel('NGS Read Counts')
            plt.ylabel('No. of varients (%)')
            plt.ylim(0,32)
            plt.yticks(np.arange(0,32+8,8))
            if save == True:
                plt.savefig( "../data/Result/read_count_distribution/"+"read_count_distribution_"+self.name+"_"+self.replicate+".svg" )
                plt.savefig( "../data/Result/read_count_distribution/"+"read_count_distribution_"+self.name+"_"+self.replicate+".tif" )
            plt.show()
        else:
            print('No Counts Data!')
        return n,bins,patches

    def max_separation_of_nt_distribution(self, order,num_buckets = 10, x_max = 0.34, y_min = 0, y_max = 6, h_flip = False ,save = False):
        '''
        Notice that num_buckets is the number of buckets on one side of distribution
        '''
        '''
        use color to show the trend ?
        '''
        '''
        Currently ylim set 1~7 (max 7.0x)
        '''
        length_list = np.zeros((2,num_buckets))
        M = x_max #max(max(self.epistasis[order].values()),abs(min(self.epistasis[order].values())))
        boundary = np.arange(0, M+M/num_buckets , M/num_buckets)
        for i in range(len(boundary)-1):
            #print('Bucket boundary {},{}'.format(boundary[i],boundary[i+1]))
            bucket = list(filter(lambda x: self.epistasis[order][x] > boundary[i] and self.epistasis[order][x] <= boundary[i+1], self.epistasis[order].keys()))
            length_list[0][i] = np.mean(list(map(lambda x: util.max_separation_of_nt(x),bucket)))
            bucket2 = list(filter(lambda x: self.epistasis[order][x] < -boundary[i] and self.epistasis[order][x] >= -boundary[i+1], self.epistasis[order].keys()))
            length_list[1][i] = np.mean(list(map(lambda x: util.max_separation_of_nt(x),bucket2)))
    
        new_length_list = np.append(length_list[1][::-1],length_list[0])
        
        if h_flip == True:
            new_length_list = new_length_list[::-1]
        #print("The shape of list:", new_length_list.shape)

        
        pos = (boundary[:-1] + boundary[1:])/2
        pos = np.append(-pos[::-1],pos)
        
        
        x = np.array(list(self.epistasis[order].values()))
        if h_flip == True:
            x = x * (-1)
        y = np.random.rand(len(x))
        plt.figure(figsize=(6,6/8))
        # add errorbar
        a = np.quantile(x, 0.25)
        b = np.quantile(x, 0.75)
        print(a,b)
        plt.errorbar([0],[0.5],xerr=[[abs(a)],[b]],fmt='none',ecolor='red',elinewidth = 4, capsize = 12, capthick = 2)
        plt.scatter(x,y, color = 'k', s = 15) # point larger, orignally 5 -->
        plt.xlim(-x_max,x_max)
        plt.xticks(np.arange(-x_max,x_max+x_max/2, x_max/2))
        plt.ylim(-0.15,1.15)
        if save == True:
            plt.savefig( "../data/Result/max_separation_of_nt_distribution/" + "scatter_distance_"+self.name+"_"+self.replicate+"_"+str(order)+"_"+"layer1"+".tif")
        plt.show()

        plt.figure(figsize=(6,6/8))
        mask = np.isfinite(new_length_list)
        
        

        print('height')
        print(new_length_list[mask])
        
        '''
        for plotting purpouse
        '''
        new_length_list2 = []
        for i in new_length_list[mask]:
            if i <7:
                if i > 1:
                    new_length_list2.append(i-1)
                else:
                    new_length_list2.append(0)
            else:
                new_length_list2.append(6)
                
        new_length_list_new = new_length_list2
        
        plt.bar(pos[mask], new_length_list_new, width = x_max/num_buckets)  # bin width is still subject to change
        plt.xlim(-x_max,x_max)
        plt.ylim(y_min,y_max)
        plt.xticks(np.arange(-x_max,x_max+x_max/2, x_max/2))
        plt.yticks([y_min, y_max])
        if save == True:
            plt.savefig( "../data/Result/max_separation_of_nt_distribution/" + "scatter_distance_"+self.name+"_"+self.replicate+"_"+str(order)+"_"+"layer2"+".svg")
        plt.show()

        '''
        gs = gridspec.GridSpec(2, 1)
        ax1 = plt.subplot(gs[0, :])
        ax2 = plt.subplot(gs[1, :])
        
        # usage:xaxis.set_ticks(np.arange(start, end, stepsize))
        # plot layer 1 : data distribution
        x = list(self.epistasis[order].values())
        y = np.random.rand(len(x))
        ax1.scatter(x,y, color = 'k', s = 1)
        ax1.set_xlim([-x_max,x_max])
        
        # plot layer 2 : bar plot
        mask = np.isfinite(new_length_list)
        ax2.bar(pos[mask], new_length_list[mask], width = x_max/num_buckets)  # bin width is still subject to change
        ax2.set_xlim([-x_max,x_max])
        ax2.set_ylim([0,y_max])
        #ax2.set_ylabel('avg. max separation of n.t.')
        
        if save == True:
            plt.savefig( "../data/Result/max_separation_of_nt_distribution/" + "data_scatter_distance_annotation"+self.name+"_"+self.replicate+"_"+str(order)+".svg")
        plt.show()
        '''
        return pos, new_length_list
    
    def violin_epistasis(self, index = [1,2,3,4], save = False):
        '''
        half violin
        half histogram
        '''
        if len(self.epistasis[0])== 0:
            print('Please Read in epistasis table first!')
        else:            
            tmp = [ list(self.epistasis[i-1].values()) for i in index]
            plt.violinplot(tmp, vert=False, widths=0.8,
                      showmeans=True, showextrema=True, showmedians=True)
            plt.title(self.name + ' violin plot')
            plt.yticks(range(1,len(index)+1),['e'+str(i) for i in index])
            #plt.xticks(ind, ['R'+str(i) for i in range(1,21)])
            #plt.yticks(np.arange(-0.35, 0.35, 0.7))
            plt.xlim(-0.36,0.36)
            if save == True:
                plt.savefig("../data/Result/violin_epistasis/"+ "violinplot_"+self.name+"_"+self.replicate+".svg")
            plt.show()        
    
    def boxplot_epistasis(self,index = [1,2,3,4], xmin = -0.36, xmax = 0.36, save = False ):
        if len(self.epistasis[0])== 0:
            print('Please Read in epistasis table first!')
        else:            
            tmp = [ list(self.epistasis[i-1].values()) for i in index]
            plt.boxplot(tmp, vert = False )
            plt.title(self.name + ' boxplot')
            plt.yticks(range(1,len(index)+1),['e'+str(i) for i in index])
            #plt.xticks(ind, ['R'+str(i) for i in range(1,21)])
            #plt.yticks(np.arange(-0.35, 0.35, 0.7))
            plt.xlim(xmin,xmax)
            if save == True:
                plt.savefig("../data/Result/boxplot_epistasis/"+"boxplot_"+self.name+"_"+self.replicate+".svg")
            plt.show()
            
    def order1_epistasis_viz(self, ylim_low = -0.34, ylim_high = 0.34, upSidedown = False ,save = False):
        '''
        fepb: max epistasis 0.308 --> take 0.32 as ylim
        '''
        x = range(-self.spacing-9,-self.spacing)
        y = np.array(list(self.epistasis[0].values())).reshape(9,4)
        if upSidedown == True:
            y *= -1
        plt.plot(x,y[:,0],'g.-',label = 'A')
        plt.plot(x,y[:,1],'b.-',label = 'C')
        plt.plot(x,y[:,2],'k.-',label = 'G')
        plt.plot(x,y[:,3],'r.-',label = 'U')
        plt.ylim(ylim_low, ylim_high)
        plt.yticks(np.arange(ylim_low, ylim_high+(ylim_high-ylim_low)/4, (ylim_high-ylim_low)/4))
        if self.name == 'vienna':
            plt.xticks(list(range(-11,-20,-1)))
        else:
            plt.xticks(list(range(-4,-16,-1)))
        #plt.legend(title='n.t.')
        if save == True:
            plt.savefig("../data/Result/order1_epistasis_viz/"+ "o1_"+ self.name +"_"+self.replicate+".svg")
        plt.show()
        
    def order2_epistasis_viz(self, save = False):
        tmp = list(self.epistasis[1].values())
        matrix = np.zeros((36,36))
        
        index = 0 
        for sitei in range(8):
            for sitej in range(sitei+1,9):
                for ni in range(0,4):
                    for nj in range(0,4):
                        matrix[sitei*4+ni,sitej*4+nj] = tmp[index]
                        index += 1        
        matrix = matrix + matrix.T
        
        plt.imshow(matrix, cmap = 'bwr', vmax = np.max(abs(matrix)), vmin = -np.max(abs(matrix)))
        plt.colorbar()
        plt.gca().invert_yaxis()
        if save == True:
            plt.savefig("../data/Result/order2_epistasis_viz/"+ "o2_"+ self.name + "_" + self.replicate +".svg")
        plt.show()
        # To Do : diagonal --> gray
    def fitness_graph_out_degree_distribution(self):
        pass
    
    def fitness_landscape_2D(self, seq ,cumulative = False, data_type = 'layer',save =False):
        '''
        data_type = 'layer' (epistasis) or 'sum' (all order effect)
        '''
        
        g  = nx.Graph()
        variant_set = list(util.mutational_set(seq))
        locs = [i for i in range(len(seq)) if seq[i] != 'N']     
        
        # for real  fitness landscape
        if seq.count('N') == 0:
            for ele in variant_set:
                g.add_node(ele, fitness= self.dic[ele])
        # for motif fitness landscape
        else:
            if np.any([len(self.epistasis[i]) > 0 for i in range(9)]) == False:
                print('Have Not Load Epistasis information !!')
            else:
                for ele in variant_set:
                    if data_type == 'layer':
                        g.add_node(ele, fitness = self.epistasis[len(locs)-1][ele])
                    elif data_type == 'sum':
                        f = self.mean
                        for i in range(0,len(locs)):
                            f +=  np.sum([self.epistasis[ len(locs)-i -1 ][subseq] for subseq in util.SubSeq(ele, i)])
                        print(f)
                        g.add_node(ele, fitness = f)
                    else:
                        print('Not valid data type')
                        return 0
        
        for seq in variant_set:
            for i in util.neighbors(seq): # problematic ??
                g.add_edge(seq,i)
        pos = nx.spring_layout(g)
        label_dict = {i: ''.join([i[loc] for loc in locs]) for i in g.nodes()}
        nx.draw_networkx_labels(g, pos, label_dict, font_size=10)
        nx.draw(g,pos)
        plt.show()
        return g, pos , label_dict
    
    def fitness_landscape_3D(self, seq, pos = None, data_type = 'layer', angle = 100 , vmax = 0.2, vmin = -0.2, elev = 3, azim=30 ,save=False):
        if pos == None:
            g, pos , label_dict = self.fitness_landscape_2D(seq, data_type = data_type)
        else:
            g, pos_not_used, label_dict = self.fitness_landscape_2D(seq, data_type = data_type)
            
        #att = np.array(list(nx.get_node_attributes(g, 'fitness').values()))
        #with plt.style.context(('ggplot')):        
        fig = plt.figure(figsize=(10,7))
        ax = Axes3D(fig)
        # Loop on the pos dictionary to extract the x,y,z coordinates of each node
        x_vec = []
        y_vec = []
        z_vec = []
        c_vec = []
        for key, value in pos.items():
            xi = value[0]
            yi = value[1]
            zi = g.node[key]['fitness']
            x_vec.append(xi)
            y_vec.append(yi)
            z_vec.append(zi)
            c_vec.append(zi)
            
            if ('GAG' in key) or ('GUG' in key) or ('GCG' in key) or ('GGG' in key):
            #if key.count('G')>=2:
            #if key[start]=='G' and key[start+num_nonN-1]=='G':
                ax.text(xi,yi,zi+0.01, label_dict[key], color ='red', fontsize= 16 )
            # Scatter plot
            #'viridis' 'inferno'  'YlGnBu' 'gist_rainbow'
            #color = plt.cm.bwr( int(round((zi-min(att))*255/max((att-min(att))))))
            
            
            #print(color)
            #PCM=ax.get_children()[2] #get the mappable, the 1st and the 2nd are the x and y axes
            #plt.colorbar(PCM, ax=ax) 
        #color = plt.cm.bwr( int(round((zi-min(att))*255/max((att-min(att))))))
        norm = mpl.colors.Normalize(vmin=vmin,vmax = vmax)
        
        p = ax.scatter(x_vec, y_vec, z_vec,c=z_vec ,cmap ='bwr', norm=norm, vmax = vmax, vmin = vmin , s=90, edgecolors='k', alpha=0.9)  #, c=zi   
        # vmin = -0.2,vmax = 0.2
        
        # Loop on the list of edges to get the x,y,z, coordinates of the connected nodes
        # Those two points are the extrema of the line to be plotted
        #plt.colorbar()
        for i,j in g.edges():
            x = np.array((pos[i][0], pos[j][0]))
            y = np.array((pos[i][1], pos[j][1]))
            z = np.array((g.node[i]['fitness'], g.node[j]['fitness']))        
        # Plot the connecting lines
            ax.plot(x, y, z, c='tab:gray', alpha=0.2)
        
        ax.grid(False)
        ax.set_zlim(vmin,vmax)
        ax.set_xticks([]) 
        ax.set_yticks([]) 
        ax.set_zticks([])
        
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0))
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0))
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0))

        ax.xaxis._axinfo['juggled'] = (0,0,0)
        ax.yaxis._axinfo['juggled'] = (1,1,0)
        ax.zaxis._axinfo['juggled'] = (2,2,0)
        
        # Set the initial view
        ax.view_init(elev=elev, azim=azim)
        print(sorted(z_vec, reverse = True))
        # Hide the axes
        #ax.set_axis_off()
        fig.colorbar(p, shrink=0.5) # , ticks=np.arange(-0.2,0.2+0.05,0.05)
        
        '''
        cmap = mpl.cm.bwr
        norm = mpl.colors.Normalize(vmin=-0.2, vmax=0.2)
        cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                        norm=norm,
                                        orientation='vertical', ax = ax)
        cb1.set_label('Some Units')
        '''
        
        
        
        if save is not False:
            plt.savefig("../data/Result/fitness_landscape_3D/" + "GNGplot_"+self.name+"_"+ data_type +".tif") # str(angle).zfill(3)
            plt.close('all')
        else:
            #plt.title('seq = '+ seq)
            plt.show()


    
    def out_degree_distribution(self,rand = False):
        m_error = self.std
        print("Threshold =", m_error)
        hist_array = np.zeros(28)
        tmp = list()
        val_count = 0
        
        if rand == True:
            tt = np.random.permutation(self.dataframe['LogMean'])
            dicc = {i:j for i, j in zip(self.dataframe["SeqID"],tt)}
        else:
            dicc = self.dic
        
        for i in dicc.keys():
            count = [0,0,0]
            if (True not in [np.isnan(dicc[j]) for j in util.neighbors(i)]) and (not np.isnan(dicc[i])):
                val_count += 1
                for j in util.neighbors(i):
                    dif = dicc.get(j,np.nan) - dicc[i]
                    if dif > m_error:
                        count[0] += 1
                    elif dif < -m_error:
                        count[2] += 1
                    else:
                        count[1] +=1
                hist_array[count[0]] += 1
                tmp.append(tuple(count))                
            #else:
            #    continue
            
        tmp1 = Counter(tmp)
        print(tmp1)


        h = np.zeros((28,28))
        
        for i in tmp1.keys():
            h[i[0]][i[2]] = np.log10(tmp1[i])
        print("# of variants included = ",val_count)
        print("# of variants detected = ",self.num_detected)
        print('ratio =',val_count/self.num_detected)
        print('# of variants affected by a noninclusion event =', (self.num_detected-val_count)/(4**9-self.num_detected))
        print("Number of Peaks= ",10**h[0][27])
        print("Number of Valleys= ",10**h[27][0])
    
        util.hm(h,[str(i) for i in range(0,28)],[str(i) for i in range(0,28)],title='', cbarlabel = 'log(count)')
        #plt.savefig("roughness_plot_error2_"+df.name+".svg")
        #heatmap(h,[str(i) for i in range(1,28)],[str(i) for i in range(0,28)])
        #plt.scatter(X,Y, s=area, c=colors, alpha=0.5)
        plt.gca()      
        #return h,hist_array

    def out_degree_distribution_trend(self,bins=20, mini=0, maxi = 3.2, k=1 ,threshold = None, error_bar = False, Alpha = 1,ticks = 5, legend = False,save = False):
        self.bucket_classifier(mini = mini, maxi = maxi, number_bins = bins)  # make bucket dict
        try:
            m_error = k * self.std
            print("Threshold =", m_error)
            
        except AttributeError:
            if threshold:
                m_error = threshold
            else:
                print("No threshold !!!!")
        '''
        if rand == True:
            tt = np.random.permutation(self.dataframe['LogMean'])
            dicc = {i:j for i, j in zip(self.dataframe["SeqID"],tt)}
        else:
        '''
        dicc = self.dic
        
        #buckets = np.zeros((3,bins)) # First 3 rows --> up/down/equal
        up_bucket = [[] for i in range(bins)]
        down_bucket = [[] for i in range(bins)]
        equal_bucket = [[] for i in range(bins)]
        val_count = np.zeros(bins) # valid count for each bin
        
        for i in dicc.keys():
            if (True not in [np.isnan(dicc[j]) for j in util.neighbors(i)]) and (not np.isnan(dicc[i])):
                ind = self.bucket_dict[i]
                val_count[ind] += 1
                tmp = [0,0,0] #up/down/equal
                for j in util.neighbors(i):
                    dif = dicc.get(j,np.nan) - dicc[i]
                    if dif > m_error:
                        tmp[0] += 1
                    elif dif < -m_error:
                        tmp[1] += 1
                    else:
                        tmp[2] +=1
                
                up_bucket[ind].append(tmp[0])
                down_bucket[ind].append(tmp[1])
                equal_bucket[ind].append(tmp[2])
                
            #else:
            #    continue
        
        mean_bucket = np.zeros((3,bins))
        std_bucket = np.zeros((3,bins))
        
        
        for i in range(bins):
            mean_bucket[0][i] = np.mean(up_bucket[i])
            mean_bucket[1][i] = np.mean(down_bucket[i])
            mean_bucket[2][i] = np.mean(equal_bucket[i])
            std_bucket[0][i] = np.std(up_bucket[i]) 
            std_bucket[1][i] = np.std(down_bucket[i])
            std_bucket[2][i] = np.std(equal_bucket[i])

        # For debugging

        
        #print(mean_bucket[0][17])
        #print(mean_bucket[0][18])
        #print(mean_bucket[0][19])
            
        # For debugging
        #print(mean_bucket)
        #print(std_bucket)
        
        x = np.arange(mini,maxi+maxi/bins,maxi/bins)
        x = (x[:-1] + x[1:])/2
        print('length of x =', len(x))
        print('x', x)
        
        y1 = mean_bucket[0]
        mask1 = np.isfinite(y1)
        plt.plot(x[mask1], y1[mask1], color='red', linestyle='-',linewidth=3, label = 'up') #marker='o'
                
        y2 = mean_bucket[1]
        mask2 = np.isfinite(y2)        
        plt.plot(x[mask2], y2[mask2], color='blue',linestyle='-',linewidth=3, label = 'down') #  marker='o'

        y3 = mean_bucket[2]
        mask3 = np.isfinite(y3)
        plt.plot(x[mask3], y3[mask3], color='black',linestyle='-',linewidth=3, label = 'equal') #  marker='o'  
        
        if error_bar == True:
            # fmt='.k'  #color='black', xecolor='lightgray', elinewidth=3, capsize=0  
            # fmt = '.'
            plt.errorbar(x[mask1]+0.03, y1[mask1], yerr = std_bucket[0][mask1], fmt='none' ,ecolor='red',elinewidth = 2, capsize = 3, alpha = Alpha)
            plt.errorbar(x[mask2]-0.03, y2[mask2], yerr = std_bucket[1][mask2], fmt='none' ,ecolor='blue', elinewidth= 2, capsize = 3,alpha = Alpha)
            plt.errorbar(x[mask3], y3[mask3], yerr = std_bucket[2][mask3], fmt='none' ,ecolor='black',elinewidth = 2, capsize = 3, alpha = Alpha)        

        plt.yticks(np.arange(0,27+1,3))
        
        plt.xticks(np.arange(mini,maxi+(maxi-mini)/(ticks-1),(maxi-mini)/(ticks-1)))
        if legend == True:
            plt.legend()
        if save == True:
            plt.savefig("../data/Result/out_degree_distribution_trend/"+"roughness_plot_trend_"+self.name+ "_" +self.replicate + "_" + str(k) + ".svg")
        plt.show()
        
        print("# of variants included = ", np.sum(val_count))
        print("# of variants detected = ",self.num_detected)
        print('ratio =', np.sum(val_count)/self.num_detected)
        print('# of variants affected by a noninclusion event =', (self.num_detected-np.sum(val_count))/(4**9-self.num_detected))
              
        return mean_bucket, val_count

    def out_neighbor_moment_trend(self,bins=20, mini=0, maxi = 3.2, s=1, moment = 1, Alpha = 1,ticks = 5, rand = False, legend = False, regress_line = False,save = False):        
        if rand == True:
            tt = np.random.permutation(self.dataframe['LogMean'])
            dicc = {i:j for i, j in zip(self.dataframe["SeqID"],tt)}
        else:
            dicc = self.dic
        '''
        elif test == True:
            #from scipy.stats import skewnorm
            #a = 4
            #mean, var, skew, kurt = skewnorm.stats(a, moments='mvsk')
            #xx = skewnorm.rvs(a, size=262144)
            
            #mu, sigma = 0.8, 0.1
            #xx = np.random.normal(mu, sigma, 262144)
            #dicc = {i:j for i, j in zip(self.dataframe["SeqID"],xx)}
            #plt.hist(xx)
            #plt.show()

            #xx = np.random.uniform(0.99, 1.0, 262144)
            #dicc = {i:j for i, j in zip(self.dataframe["SeqID"],xx)}
        '''
        

        
        #buckets = np.zeros((3,bins)) # First 3 rows --> up/down/equal
        #up_bucket = [[] for i in range(bins)]
        #down_bucket = [[] for i in range(bins)]
        #equal_bucket = [[] for i in range(bins)]
        ori = []
        neighbor_mean = []
        neighbor_std = []
        val_count = 0
        
        for i in dicc.keys():
            if (True not in [np.isnan(dicc[j]) for j in util.neighbors(i)]) and (not np.isnan(dicc[i])):
                val_count += 1
                tmp = np.zeros(27)
                for index,j in enumerate(util.neighbors(i)):
                    tmp[index] = dicc[j]
                #print(tmp)
                ori.append(dicc[i])
                neighbor_mean.append(np.nanmean(tmp))
                neighbor_std.append(np.nanstd(tmp))
        
        #print("max in mean =", max(neighbor_mean))
        #print("min in mean =", min(neighbor_mean))
        plt.figure(figsize=(6,6))
        if moment == 1:
            #    plt.scatter(f,a, c = 'g',s = 9, alpha=1)
            #plt.hist(neighbor_mean)
            #plt.show()
            plt.scatter(ori,neighbor_mean, c = 'g', s = s, alpha = Alpha)
            X = np.array(ori).reshape(-1, 1)
            Y = np.array(neighbor_mean).reshape(-1, 1)
        elif moment == 2:
            plt.scatter(ori,neighbor_std, c = 'g', s = s, alpha = Alpha)
            X = np.array(ori)
            Y = np.array(neighbor_std).reshape(-1, 1)
            
        plt.xlim(mini,maxi)
        plt.ylim(mini,maxi)
        plt.xticks(np.arange(mini,maxi+(maxi-mini)/(ticks-1),(maxi-mini)/(ticks-1)))
        plt.yticks(np.arange(mini,maxi+(maxi-mini)/(ticks-1),(maxi-mini)/(ticks-1)))
        
        if regress_line == True:
            from sklearn.linear_model import LinearRegression
            linear_regressor = LinearRegression()  # create object for the class
            linear_regressor.fit(X, Y)  # perform linear regression
            #Y_pred = linear_regressor.predict(X) # make predictions
            print("m = ", linear_regressor.coef_)
            print("b = ", linear_regressor.intercept_)
            print("R^2 score = ", linear_regressor.score(X,Y))

            plt.plot([0,maxi],[linear_regressor.intercept_,linear_regressor.coef_*maxi+linear_regressor.intercept_],'r--')
            
        
        if legend == True:
            plt.legend()
        if save == True:
            if rand == True:
                plt.savefig("../data/Result/out_degree_moment_trend/"+"out_degree_moment_trend_"+self.name+ "_" +self.replicate + "_" +"size" +str(s) + "_moment" + str(moment) + "_rand"+".tif")
            else:
                plt.savefig("../data/Result/out_degree_moment_trend/"+"out_degree_moment_trend_"+self.name+ "_" +self.replicate + "_" +"size" +str(s) + "_moment" + str(moment) +".tif")
                
        plt.show()
        
        print(">>> Correlation coefficient >>>")
        if moment == 1:
            cor = util.correlation_coef(ori,neighbor_mean)
            print(cor)
            cor2, p = stats.pearsonr(ori,neighbor_mean)
            print(cor2, p)
        elif moment ==2:
            cor = util.correlation_coef(ori,neighbor_std)
            print(cor)
            cor2, p = stats.pearsonr(ori,neighbor_std)
            print(cor2, p)            

        
        print("# of variants included = ", val_count)
        print("# of variants detected = ",self.num_detected)
        print('ratio =', val_count/self.num_detected)
        print('# of variants affected by a noninclusion event =', (self.num_detected- val_count)/(4**9-self.num_detected))
        return cor2            
        #return ori, neighbor_mean, neighbor_std




    def plot_all_shortest_paths(self,s1,s2,save = False):
        
        def non_decreasing(L):
            return all(x<=y for x, y in zip(L, L[1:]))
    
        difference_index = [i for i in range(len(s1)) if s1[i] != s2[i]]
        max_steps = len(difference_index)
        counter = 0
        all_shortest_paths = util.all_shortest_paths(s1,s2)
        total_shortest_paths = len(all_shortest_paths)
        print("Total # of shortest paths = ", total_shortest_paths)
        X = range(max_steps+1)
        for path in all_shortest_paths:
            # check whether each path satisfies non-decreasing condition
            tmp = [self.dic[i] for i in path]    
            plt.plot(X,tmp)
            if non_decreasing(tmp):
                counter += 1
        plt.xticks(np.arange(max_steps+1))
        if save == True:
            plt.savefig("../data/Result/shortest_paths/" +self.name+"shortest_paths.tif")
        plt.show()
        
        print("# of viable shortest paths =", counter)
        print("ratio = ", counter/total_shortest_paths)
        
    def out_degree_moment_trend_plot_all(self,bins=20, mini=0, maxi = 3.2, s=1, moment = 1, Alpha = 1,ticks = 5, rand = False, legend = False, regress_line = False,density_estimate = False,save = False):        
        if rand == True:
            tt = np.random.permutation(self.dataframe['LogMean'])
            dicc = {i:j for i, j in zip(self.dataframe["SeqID"],tt)}
        else:
            dicc = self.dic

        ori = []
        neighbor = []
        val_count = 0
        
        for i in dicc.keys():
            if (True not in [np.isnan(dicc[j]) for j in util.neighbors(i)]) and (not np.isnan(dicc[i])):
                val_count += 1
                for index,j in enumerate(util.neighbors(i)):
                    ori.append(dicc[i])
                    neighbor.append(dicc[j])
        plt.figure(figsize=(6,6))
        if density_estimate == True:
            
            # Calculate the point density
            from scipy.stats import gaussian_kde
            # random sample 1/10 of the population to plot
            x = np.array(ori)
            y = np.array(neighbor)
            print(len(x),len(y))
            idxx = x>0.8
            x = x[idxx>0.8]
            y = y[idxx>0.8]
            seeds =  np.random.choice(len(x), len(x)//50)
            x = x[seeds]
            y = y[seeds]
            xy = np.vstack([x,y])
            z = gaussian_kde(xy)(xy)
            idx = z.argsort()
            x, y, z = x[idx], y[idx], z[idx]
            plt.scatter(x, y, c=z, s=0.1, alpha = Alpha, edgecolor='')
        else:
            plt.scatter(ori,neighbor, c = 'g', s = s, alpha = Alpha)
        
        X = np.array(ori).reshape(-1, 1)
        Y = np.array(neighbor).reshape(-1, 1)
            
        plt.xlim(mini,maxi)
        plt.ylim(mini,maxi)
        plt.xticks(np.arange(mini,maxi+(maxi-mini)/(ticks-1),(maxi-mini)/(ticks-1)))
        plt.yticks(np.arange(mini,maxi+(maxi-mini)/(ticks-1),(maxi-mini)/(ticks-1)))
        
        if regress_line == True:
            from sklearn.linear_model import LinearRegression
            linear_regressor = LinearRegression()  # create object for the class
            linear_regressor.fit(X, Y)  # perform linear regression
            #Y_pred = linear_regressor.predict(X) # make predictions
            print("m = ", linear_regressor.coef_)
            print("b = ", linear_regressor.intercept_)
            print("R^2 score = ", linear_regressor.score(X,Y))

            plt.plot([0,maxi],[linear_regressor.intercept_,linear_regressor.coef_*maxi+linear_regressor.intercept_],'r--')
        
        if legend == True:
            plt.legend()
        if save == True:
            if rand == True:
                plt.savefig("../data/Result/out_degree_moment_trend/"+"out_degree_moment_trend_plotall"+self.name+ "_" +self.replicate + "_" +"size" +str(s) + "_moment" + str(moment) + "_rand"+".tif")
            else:
                plt.savefig("../data/Result/out_degree_moment_trend/"+"out_degree_moment_trend_plotall"+self.name+ "_" +self.replicate + "_" +"size" +str(s) + "_moment" + str(moment) +".tif")
                
        plt.show()
        
        print(">>> Correlation coefficient >>>")
        if moment == 1:
            cor = util.correlation_coef(ori,neighbor)
            print(cor)
            cor2, p = stats.pearsonr(ori,neighbor)
            print(cor2, p)

        
        print("# of variants included = ", val_count)
        print("# of variants detected = ",self.num_detected)
        print('ratio =', val_count/self.num_detected)
        print('# of variants affected by a noninclusion event =', (self.num_detected- val_count)/(4**9-self.num_detected))
        return cor2   
    
def main():
    pass
        
