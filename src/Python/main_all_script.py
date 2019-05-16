#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 09:24:30 2019

@author: rljahn
"""

from library import *



arti = pd.read_csv("../../data/Read Count Table/UnionTable/arti_SDR_union_count25.csv")
fepb = pd.read_csv("../../data/Read Count Table/UnionTable/fepB_SDR_union_count25.csv")
dmsc = pd.read_csv("../../data/Read Count Table/UnionTable/dmsC_SDR_union_count25.csv")
yfp = pd.read_csv("../../data/Read Count Table/yfp_SDR_rep1.csv") 
vienna = pd.read_csv("../../data/Vienna Package Prediction/Vienna_SeqExp.csv")



if __name__ == '__main__':
    
    # create library objects
    arti = library("arti",'union', arti)
    fepb = library("fepB",'union', fepb)
    dmsc = library("dmsC",'union', dmsc)
    vienna = library("vienna",'',vienna)
    
    yfp = pd.read_csv("../../data/Read Count Table/yfp_SDR_rep1.csv")
    yfp.loc[yfp.Counts < 25,'LogMean'] = np.nan
    yfp.loc[yfp.Counts < 25,'Counts'] = np.nan
    yfp = library("yfp",'rep1',yfp)


    # visualize data distribution  (Fig. 1B, S20A)
    '''
    fepb.distribution(normed = True)
    arti.distribution(normed = True)
    dmsc.distribution(normed = True)
    #vienna.distribution(normed = True, mini = min(vienna.dic.values())-0.5, maxi = max(vienna.dic.values())+0.5)
    yfp.distribution(normed = True, mini = 0, maxi = 2.7, ticks = 4)
    '''


    
    # compare the distribution of 3 libraries using pairwise scatterplots (Fig. 1C)
    '''
    print(">>>> 3 libraries pairwise scatterplots >>>>>")
    f = fepb.dataframe['LogMean']
    a = arti.dataframe['LogMean']
    d = dmsc.dataframe['LogMean']
    v = vienna.dataframe['LogMean']
    v = np.abs(v)
    ind_f = ~np.isnan(f)
    ind_a = ~np.isnan(a)
    ind_d = ~np.isnan(d)
    ind_v = ~np.isnan(v)


    ind = np.logical_and(ind_f, ind_a)
    print("Total data used = ", len(f[ind]))

    result = stats.pearsonr(f[ind],a[ind])
    print("Pearson-r:", result[0])  
    print("p_value:", result[1])
    
    result = stats.spearmanr(f[ind],a[ind])
    print("Spearman-r:", result[0])  
    print("p_value:", result[1])
    
    plt.figure(figsize=(6,6))
    plt.scatter(f,a, c = 'g',s = 9, alpha=1)
    plt.xlim(0,3.2)
    plt.ylim(0,3.2)
    #plt.savefig("../data/Result/pairwise_scatterplots/" + 'x_fepb_y_arti - Spearmanr_{}.tif'.format(result[0]))
    plt.show()
    
    ind = np.logical_and(ind_a, ind_d)
    # Method 1
    # r2 = r2_score(a[ind],d[ind])
    # print("r-squared:", r2) 
    # Method 2
    #slope, intercept, r_value, p_value, std_err = stats.linregress(a[ind],d[ind])
    print("Total data used = ", len(a[ind]))
    result = stats.pearsonr(a[ind],d[ind])
    print("Pearson-r:", result[0])  
    print("p_value:", result[1])
    result = stats.spearmanr(a[ind],d[ind])
    print("Spearman-r:", result[0])  
    print("p_value:", result[1])
    plt.figure(figsize=(6,6))
    plt.scatter(a,d, c = 'g',s = 9, alpha=1)
    plt.xlim(0,3.2)
    plt.ylim(0,3.2)
    #plt.savefig("../data/Result/pairwise_scatterplots/" + 'x_arti_y_dmsc - Spearmanr_{}.tif'.format(result[0]))
    plt.show()

    ind = np.logical_and(ind_d, ind_f)
    # r2 = r2_score(d[ind],f[ind])
    # print("r-squared:", r2)
    #slope, intercept, r_value, p_value, std_err = stats.linregress(d[ind],f[ind])
    print("Total data used = ", len(d[ind]))
    result = stats.pearsonr(d[ind],f[ind])    
    print("Pearson-r:", result[0])  
    print("p_value:", result[1])
    result = stats.spearmanr(d[ind],f[ind])
    print("Spearman-r:", result[0])  
    print("p_value:", result[1])
    plt.figure(figsize=(6,6))
    plt.scatter(d,f, c = 'g',s = 9, alpha=1)
    plt.xlim(0,3.2)
    plt.ylim(0,3.2)
    #plt.savefig("../data/Result/pairwise_scatterplots/" + 'x_dmsc_y_fepb - Spearmanr_{}.tif'.format(result[0]))
    plt.show()
    '''


    # fitness of Reference genotype vs mean fitness of its neighbors (Fig 1E)
    '''
    arti.out_neighbor_moment_trend(moment=1,s=0.5)
    fepb.out_neighbor_moment_trend(moment=1,s=0.5)
    dmsc.out_neighbor_moment_trend(moment=1,s=0.5)
    yfp.out_neighbor_moment_trend(mini=0, maxi = 2.7, ticks = 4, moment=1,s=0.5)
    '''
    

    # visualize trend in out-degree (by different backgrounds) (Fig 1F)
    '''
    a, c1 = fepb.out_degree_distribution_trend(error_bar=True)
    a, c2 = arti.out_degree_distribution_trend(error_bar=True)
    a, c3 = dmsc.out_degree_distribution_trend(error_bar=True)
    yfp.out_degree_distribution_trend(error_bar=True,threshold = 0.108,maxi = 2.7, ticks = 4)
    '''

    # nucleotide composition (Fig 2A)
    """
    fepb.nucleotide_composition()
    arti.nucleotide_composition()
    dmsc.nucleotide_composition()
    yfp.nucleotide_composition( mini = 0, maxi = 2.7, ticks = 4)  
    """

    # nucleotide trend (Fig 2B)
    """
    fepb.nucleotide_trend()
    arti.nucleotide_trend()
    dmsc.nucleotide_trend()   
    yfp.nucleotide_trend(ymax = 1.8, yticks = 4)
    """

    # Read the table and store into the library class    
    #if epistasis table does not exist --> create their epistasis tables
    #util.generate_epistasis_table(readcount_u_path,'arti_SDR_union_count25.csv',epi_path ,MaxOrder=9)
    arti.readIn_epistasis_table()
    fepb.readIn_epistasis_table()
    dmsc.readIn_epistasis_table()
    vienna.readIn_epistasis_table()
    yfp.readIn_epistasis_table()



    # visualize single nucleotide contribution (Fig 2C,D)
    '''
    fepb.order1_epistasis_viz()
    arti.order1_epistasis_viz()
    dmsc.order1_epistasis_viz()
    # Note that vienna --> sign of the values is the opposite --> use reverse argument
    vienna.order1_epistasis_viz(ylim_low = -1*(max(vienna.epistasis[0].values())+0.1) , ylim_high = -1*(min(vienna.epistasis[0].values())-0.1), upSidedown= True)
    yfp.order1_epistasis_viz(ylim_low = -0.16, ylim_high = 0.16)
    '''
    
    # visualize pairwise contribution (Fig 2E,F)
    '''
    fepb.order2_epistasis_viz()
    arti.order2_epistasis_viz()
    dmsc.order2_epistasis_viz()
    vienna.order2_epistasis_viz() # Note that the color is reversed in vienna case    
    yfp.order2_epistasis_viz()    
    '''
    


    print(">>>> Calculate Prediction value >>>>")
    arti.model()
    arti.model.readIn_predicted_npy_files()
    fepb.model()
    fepb.model.readIn_predicted_npy_files()
    dmsc.model()
    dmsc.model.readIn_predicted_npy_files()
    vienna.model()
    vienna.model.readIn_predicted_npy_files()
    yfp.model()
    yfp.model.readIn_predicted_npy_files()

    #  Explanatory power (R2) of single-nucleotide effects (N1) and pairwise epistasis (N2) on fitness (G) and the SD:aSD base-pairing energy (Fig 2G,H)
    # N1 vs N1-2 prediction scatter plots
    '''
    # arti.
    plt.figure(figsize=(6,6))
    x = np.array(list(arti.model.dic.values()))
    y = arti.model.fitted_expression_list[1-1]
    ind_x = ~np.isnan(x)
    ind_y = ~np.isnan(y)
    ind = np.logical_and(ind_x, ind_y)
    x = x[ind]
    y = y[ind]
    #print("Total data used = ", len(x))
    plt.scatter(y,x, s= 4, c = 'b', alpha=1)  #'lightcyan' 'skyblue'
    r = util.correlation_coef(y,x)
    print("r2 = :", r**2)
    print
    y2 = arti.model.fitted_expression_list[2-1]
    ind_y2 = ~np.isnan(y2)
    ind = np.logical_and(ind_x, ind_y2)
    y2 = y2[ind]
    plt.scatter(y2,x, s= 4, c = 'r', alpha=1) # 'lightpink' 'r'
    r = util.correlation_coef(y2,x)
    print("r2 = :", r**2)
    plt.xlim(0,3.2)
    plt.ylim(0,3.2)
    plt.xticks(np.arange(0, 3.2+0.8, 0.8))
    plt.yticks(np.arange(0, 3.2+0.8, 0.8))
    #plt.savefig("../data/Result/variance_explained/" + 'arti_N1vsN1-2.tif')
    plt.show()
    
    # fepb
    plt.figure(figsize=(6,6))
    x = np.array(list(fepb.model.dic.values()))
    y = fepb.model.fitted_expression_list[1-1]
    ind_x = ~np.isnan(x)
    ind_y = ~np.isnan(y)
    ind = np.logical_and(ind_x, ind_y)
    x = x[ind]
    y = y[ind]
    plt.scatter(y,x, s= 4, c = 'b', alpha=1)  #'lightcyan' 'skyblue'
    r = util.correlation_coef(y,x)
    print("r2 = :", r**2)

    y2 = fepb.model.fitted_expression_list[2-1]
    ind_y2 = ~np.isnan(y2)
    ind = np.logical_and(ind_x, ind_y2)
    y2 = y2[ind]    
    plt.scatter(y2,x, s= 4, c = 'r', alpha=1) # 'lightpink' 'r'
    r = util.correlation_coef(y2,x)
    print("r2 = :", r**2)
    plt.xlim(0,3.2)
    plt.ylim(0,3.2)
    plt.xticks(np.arange(0, 3.2+0.8, 0.8))
    plt.yticks(np.arange(0, 3.2+0.8, 0.8))
    #plt.savefig("../data/Result/variance_explained/" + 'fepb_N1vsN1-2.tif')
    plt.show()


    # dmsc
    plt.figure(figsize=(6,6))
    x = np.array(list(dmsc.model.dic.values()))
    y = dmsc.model.fitted_expression_list[1-1]
    ind_x = ~np.isnan(x)
    ind_y = ~np.isnan(y)
    ind = np.logical_and(ind_x, ind_y)
    x = x[ind]
    y = y[ind]
    plt.scatter(y,x, s= 4, c = 'b', alpha=1)  #'lightcyan' 'skyblue'
    r = util.correlation_coef(y,x)
    print("r2 = :", r**2)
    
    y2 = dmsc.model.fitted_expression_list[2-1]
    ind_y2 = ~np.isnan(y2)
    ind = np.logical_and(ind_x, ind_y2)
    y2 = y2[ind]
    plt.scatter(y2,x, s= 4, c = 'r', alpha=1) # 'lightpink' 'r'
    r = util.correlation_coef(y2,x)
    print("r2 = :", r**2)
    plt.xlim(0,3.2)
    plt.ylim(0,3.2)
    plt.xticks(np.arange(0, 3.2+0.8, 0.8))
    plt.yticks(np.arange(0, 3.2+0.8, 0.8))
    #plt.savefig("../data/Result/variance_explained/" + 'dmsc_N1vsN1-2.tif')
    plt.show()
    
    # vienna    
    plt.figure(figsize=(6,6))
    x = np.array(list(vienna.model.dic.values()))
    y = vienna.model.fitted_expression_list[1-1]
    ind_x = ~np.isnan(x)
    ind_y = ~np.isnan(y)
    ind = np.logical_and(ind_x, ind_y)
    x = x[ind]
    y = y[ind]
    plt.scatter(np.abs(y),np.abs(x), s= 4, c = 'b', alpha=1)  #'lightcyan' 'skyblue'
    r = util.correlation_coef(np.abs(y),np.abs(x))
    print("r2 = :", r**2)

    y2 = vienna.model.fitted_expression_list[2-1]
    ind_y2 = ~np.isnan(y2)
    ind = np.logical_and(ind_x, ind_y2)
    y2 = y2[ind]
    plt.scatter(np.abs(y2),np.abs(x), s= 4, c = 'r', alpha=1) # 'lightpink' 'r'
    r = util.correlation_coef(np.abs(y2),np.abs(x))
    print("r2 = :", r**2)
    plt.xlim(0,14)
    plt.ylim(0,14)
    #plt.savefig("../data/Result/variance_explained/" + 'vienna_N1vsN1-2.tif')
    plt.show()



    plt.figure(figsize=(6,6))
    x = np.array(list(yfp.model.dic.values()))
    y = yfp.model.fitted_expression_list[1-1]
    ind_x = ~np.isnan(x)
    ind_y = ~np.isnan(y)
    ind = np.logical_and(ind_x, ind_y)
    x = x[ind]
    y = y[ind]
    #print("Total data used = ", len(x))
    plt.scatter(y,x, s= 2, c = 'b', alpha=1)  #'lightcyan' 'skyblue'
    r = util.correlation_coef(y,x)
    print("r2 = :", r**2)
    
    y2 = yfp.model.fitted_expression_list[2-1]
    ind_y2 = ~np.isnan(y2)
    ind = np.logical_and(ind_x, ind_y2)
    y2 = y2[ind]    
    plt.scatter(y2,x, s= 2, c = 'r', alpha=1) # 'lightpink' 'r'
    r = util.correlation_coef(y2,x)
    print("r2 = :", r**2)
    plt.xlim(0,2.7)
    plt.ylim(0,2.7)
    plt.xticks(np.arange(0, 2.7+0.9, 0.9))
    plt.yticks(np.arange(0, 2.7+0.9, 0.9))
    #plt.savefig("../data/Result/variance_explained/" + 'yfp_N1vsN1-2.tif')
    plt.show()
    '''


    # Robustnessof sort-seq experiments. (Fig S8B)
    '''
    dmsc_rep1 = library("dmsC",'rep1',dmsc_rep1)
    dmsc_rep2 = library("dmsC",'rep2',dmsc_rep2)
    dmsc_rep3 = library("dmsC",'rep3',dmsc_rep3)
    dmsc_rep1.distribution()
    dmsc_rep2.distribution()
    dmsc_rep3.distribution()
    '''

    # common 7 scatterplot (Fig S9 A)
    # for BCD, change the argument of get_common_loci_dict by the according index 
    '''
    from scipy import stats
    fepb.get_common_loci_dict(2,8) 
    arti.get_common_loci_dict(1,7)
    dmsc.get_common_loci_dict(0,6)
    
    f = np.array(list(fepb.common_loci_dic.values()))
    a = np.array(list(arti.common_loci_dic.values()))
    d = np.array(list(dmsc.common_loci_dic.values()))
    
    ind_f = ~np.isnan(f)
    ind_a = ~np.isnan(a)
    ind_d = ~np.isnan(d)


    ind = np.logical_and(ind_f, ind_a)
    print("Total data used = ", len(f[ind]))

    result = stats.pearsonr(f[ind],a[ind])
    print("Pearson-r:", result[0])  
    print("p_value:", result[1])
    
    result = stats.spearmanr(f[ind],a[ind])
    print("Spearman-r:", result[0])  
    print("p_value:", result[1])
    
    plt.figure(figsize=(6,6))
    plt.scatter(f,a, c = 'g',s = 9, alpha=1)
    plt.xlim(0,3.2)
    plt.ylim(0,3.2)
    #plt.savefig("../data/Result/common7_scatter/" + 'x_fepb_y_arti - Spearmanr_{}.tif'.format(result[0]))
    plt.show()
    
    ind = np.logical_and(ind_a, ind_d)
    result = stats.pearsonr(a[ind],d[ind])
    print("Pearson-r:", result[0])  
    print("p_value:", result[1])
    result = stats.spearmanr(a[ind],d[ind])
    print("Spearman-r:", result[0])  
    print("p_value:", result[1])
    plt.figure(figsize=(6,6))
    plt.scatter(a,d, c = 'g',s = 9, alpha=1)
    plt.xlim(0,3.2)
    plt.ylim(0,3.2)
    #plt.savefig("../data/Result/common7_scatter/" + 'x_arti_y_dmsc - Spearmanr_{}.tif'.format(result[0]))
    plt.show()

    ind = np.logical_and(ind_d, ind_f)
    result = stats.pearsonr(d[ind],f[ind])    
    print("Pearson-r:", result[0])  
    print("p_value:", result[1])
    result = stats.spearmanr(d[ind],f[ind])
    print("Spearman-r:", result[0])  
    print("p_value:", result[1])
    plt.figure(figsize=(6,6))
    plt.scatter(d,f, c = 'g',s = 9, alpha=1)
    plt.xlim(0,3.2)
    plt.ylim(0,3.2)
    #plt.savefig("../data/Result/common7_scatter/" + 'x_dmsc_y_fepb - Spearmanr_{}.tif'.format(result[0]))
    plt.show()
    '''

 

    # visualize trend in out-degree with data shuffling (by different backgrounds) (Fig S11ACD)
    '''
    arti_2 = pd.read_csv("../data/Read Count Table/UnionTable/arti_SDR_union_count25.csv")
    fepb_2 = pd.read_csv("../data/Read Count Table/UnionTable/fepB_SDR_union_count25.csv")
    dmsc_2 = pd.read_csv("../data/Read Count Table/UnionTable/dmsC_SDR_union_count25.csv")
    arti_2['LogMean'] = np.random.RandomState(seed=42).permutation(arti_2['LogMean'])
    arti_rand = library("arti_rand",'union',arti)
    fepb_2['LogMean'] = np.random.RandomState(seed=42).permutation(fepb_2['LogMean'])
    fepb_rand = library("fepB_rand",'union',fepb)
    dmsc_2['LogMean'] = np.random.RandomState(seed=42).permutation(dmsc_2['LogMean'])
    dmsc_rand = library("dmsC_rand",'union',dmsc)
    fepb_rand.out_degree_distribution_trend(error_bar=True)
    arti_rand.out_degree_distribution_trend(error_bar=True)
    dmsc_rand.out_degree_distribution_trend(error_bar=True)
 
    arti.out_degree_distribution_trend(k=2, error_bar=True)
    fepb.out_degree_distribution_trend(k=2, error_bar=True)
    dmsc.out_degree_distribution_trend(k=2, error_bar=True)

    fepb_rand.out_degree_distribution_trend(k=2, error_bar=True)
    arti_rand.out_degree_distribution_trend(k=2, error_bar=True)
    dmsc_rand.out_degree_distribution_trend(k=2, error_bar=True)
    '''

    # Fig S12B
    """
    arti_rand.out_neighbor_moment_trend(moment=1,s=0.5)
    fepb_rand.out_neighbor_moment_trend(moment=1,s=0.5)
    dmsc_rand.out_neighbor_moment_trend(moment=1,s=0.5)
    """


    # Relationship between the nucleotide content and the fitness of SD variants (Fig S13)
    """
    plt.figure(figsize=(24,18))
    labelsize = 14
    length = 6
    color = 'k'
    ymax = 3.2
    yticks = 5
    s = 3 #5
    capsize = 6
    markersize = 20

    ax1 = plt.subplot(3,4,1)
    ax1.tick_params(axis = 'both', which = 'major', labelsize = labelsize ,length=length)
    ax2 = plt.subplot(3,4,2,sharex = ax1 ,sharey = ax1)
    ax2.tick_params(axis = 'both', which = 'major', labelsize = labelsize ,length=length)

    ax3 = plt.subplot(3,4,3,sharex = ax1 ,sharey = ax1)
    ax3.tick_params(axis = 'both', which = 'major', labelsize = labelsize ,length=length)
    ax4 = plt.subplot(3,4,4,sharex = ax1 ,sharey = ax1)
    ax4.tick_params(axis = 'both', which = 'major', labelsize = labelsize ,length=length)
    ax5 = plt.subplot(3,4,5,sharex = ax1 ,sharey = ax1)
    ax5.tick_params(axis = 'both', which = 'major', labelsize = labelsize ,length=length)
    ax6 = plt.subplot(3,4,6,sharex = ax1 ,sharey = ax1)
    ax6.tick_params(axis = 'both', which = 'major', labelsize = labelsize ,length=length)
    ax7 = plt.subplot(3,4,7,sharex = ax1 ,sharey = ax1)
    ax7.tick_params(axis = 'both', which = 'major', labelsize = labelsize ,length=length)
    ax8 = plt.subplot(3,4,8,sharex = ax1 ,sharey = ax1)
    ax8.tick_params(axis = 'both', which = 'major', labelsize = labelsize ,length=length)
    ax9 = plt.subplot(3,4,9,sharex = ax1 ,sharey = ax1)
    ax9.tick_params(axis = 'both', which = 'major', labelsize = labelsize ,length=length)
    ax10 = plt.subplot(3,4,10,sharex = ax1 ,sharey = ax1)
    ax10.tick_params(axis = 'both', which = 'major', labelsize = labelsize ,length=length)
    ax11 = plt.subplot(3,4,11,sharex = ax1 ,sharey = ax1)
    ax11.tick_params(axis = 'both', which = 'major', labelsize = labelsize ,length=length)
    ax12 = plt.subplot(3,4,12,sharex = ax1 ,sharey = ax1)
    ax12.tick_params(axis = 'both', which = 'major', labelsize = labelsize ,length=length)
    
    
    ### fepb ###
    mean = np.zeros((4,10))
    std = np.zeros((4,10))    
    vec_expression = np.array(list(fepb.dic.values()))
    vec_A = np.array(list(map(lambda x: x.count('A') , fepb.dic.keys())))
    vec_U = np.array(list(map(lambda x: x.count('U') , fepb.dic.keys())))
    vec_C = np.array(list(map(lambda x: x.count('C') , fepb.dic.keys())))
    vec_G = np.array(list(map(lambda x: x.count('G') , fepb.dic.keys())))
    
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
        
        x = i -0.5 + 0.8*np.random.rand(len(A))
        ax1.scatter(x , A, color = color, s = s)
        # add errorbar
        a = np.nanstd(A)
        #print(a,b)
        if not (i == 9):
            ax1.errorbar([i],[mean[0][i]], yerr=[[a],[a]],fmt='r.',markersize = markersize ,ecolor='red',elinewidth = 3, capsize =  capsize, capthick = 2)
        
        x = i -0.5 + 0.8*np.random.rand(len(U))
        ax2.scatter(x , U, color = color, s = s)
        # add errorbar
        a = np.nanstd(U)
        #print(a,b)
        if not (i == 9):
            ax2.errorbar([i],[mean[1][i]], yerr=[[a],[a]],fmt='r.',markersize = markersize, ecolor='red',elinewidth = 3, capsize =  capsize, capthick = 2)        
        
        x = i -0.5 + 0.8*np.random.rand(len(C))
        ax3.scatter(x , C, color = color, s = s)
        # add errorbar
        a = np.nanstd(C)
        #print(a,b)
        if not (i == 9):
            ax3.errorbar([i],[mean[2][i]], yerr=[[a],[a]],fmt='r.',markersize = markersize,ecolor='red',elinewidth = 3, capsize =  capsize, capthick = 2)        
        
        
        x = i-0.5 + 0.8*np.random.rand(len(G))
        ax4.scatter(x , G, color = color, s = s)
        # add errorbar
        a = np.nanstd(G)
        #print(a,b)
        if not (i == 9):
            ax4.errorbar([i],[mean[3][i]], yerr=[[a],[a]],fmt='r.',markersize = markersize,ecolor='red',elinewidth = 3, capsize =  capsize, capthick = 2)
    ### arti ###
    mean = np.zeros((4,10))
    std = np.zeros((4,10))    
    vec_expression = np.array(list(arti.dic.values()))
    vec_A = np.array(list(map(lambda x: x.count('A') , arti.dic.keys())))
    vec_U = np.array(list(map(lambda x: x.count('U') , arti.dic.keys())))
    vec_C = np.array(list(map(lambda x: x.count('C') , arti.dic.keys())))
    vec_G = np.array(list(map(lambda x: x.count('G') , arti.dic.keys())))
    
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
        
        x = i -0.5 + 0.8*np.random.rand(len(A))
        ax5.scatter(x , A, color = color, s = s)
        # add errorbar
        a = np.nanstd(A)
        #print(a,b)
        if not (i == 9):
            ax5.errorbar([i],[mean[0][i]], yerr=[[a],[a]],fmt='r.',markersize = markersize ,ecolor='red',elinewidth = 3, capsize =  capsize, capthick = 2)
        
        x = i -0.5 + 0.8*np.random.rand(len(U))
        ax6.scatter(x , U, color = color, s = s)
        # add errorbar
        a = np.nanstd(U)
        #print(a,b)
        if not (i == 9):
            ax6.errorbar([i],[mean[1][i]], yerr=[[a],[a]],fmt='r.',markersize = markersize, ecolor='red',elinewidth = 3, capsize =  capsize, capthick = 2)        
        
        x = i -0.5 + 0.8*np.random.rand(len(C))
        ax7.scatter(x , C, color = color, s = s)
        # add errorbar
        a = np.nanstd(C)
        #print(a,b)
        if not (i == 9):
            ax7.errorbar([i],[mean[2][i]], yerr=[[a],[a]],fmt='r.',markersize = markersize,ecolor='red',elinewidth = 3, capsize =  capsize, capthick = 2)        
        
        
        x = i-0.5 + 0.8*np.random.rand(len(G))
        ax8.scatter(x , G, color = color, s = s)
        # add errorbar
        a = np.nanstd(G)
        #print(a,b)
        if not (i == 9):
            ax8.errorbar([i],[mean[3][i]], yerr=[[a],[a]],fmt='r.',markersize = markersize,ecolor='red',elinewidth = 3, capsize =  capsize, capthick = 2)

    ### dmsc ###
    mean = np.zeros((4,10))
    std = np.zeros((4,10))    
    vec_expression = np.array(list(dmsc.dic.values()))
    vec_A = np.array(list(map(lambda x: x.count('A') , dmsc.dic.keys())))
    vec_U = np.array(list(map(lambda x: x.count('U') , dmsc.dic.keys())))
    vec_C = np.array(list(map(lambda x: x.count('C') , dmsc.dic.keys())))
    vec_G = np.array(list(map(lambda x: x.count('G') , dmsc.dic.keys())))
    
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
        
        x = i -0.5 + 0.8*np.random.rand(len(A))
        ax9.scatter(x , A, color = color, s = s)
        # add errorbar
        a = np.nanstd(A)
        #print(a,b)
        if not (i == 9):
            ax9.errorbar([i],[mean[0][i]], yerr=[[a],[a]],fmt='r.',markersize = markersize ,ecolor='red',elinewidth = 3, capsize =  capsize, capthick = 2)
        
        x = i -0.5 + 0.8*np.random.rand(len(U))
        ax10.scatter(x , U, color = color, s = s)
        # add errorbar
        a = np.nanstd(U)
        #print(a,b)
        if not (i == 9):
            ax10.errorbar([i],[mean[1][i]], yerr=[[a],[a]],fmt='r.',markersize = markersize, ecolor='red',elinewidth = 3, capsize =  capsize, capthick = 2)        
        
        x = i -0.5 + 0.8*np.random.rand(len(C))
        ax11.scatter(x , C, color = color, s = s)
        # add errorbar
        a = np.nanstd(C)
        #print(a,b)
        if not (i == 9):
            ax11.errorbar([i],[mean[2][i]], yerr=[[a],[a]],fmt='r.',markersize = markersize,ecolor='red',elinewidth = 3, capsize =  capsize, capthick = 2)        
        
        
        x = i-0.5 + 0.8*np.random.rand(len(G))
        ax12.scatter(x , G, color = color, s = s)
        # add errorbar
        a = np.nanstd(G)
        #print(a,b)
        if not (i == 9):
            ax12.errorbar([i],[mean[3][i]], yerr=[[a],[a]],fmt='r.',markersize = markersize,ecolor='red',elinewidth = 3, capsize =  capsize, capthick = 2)

    
    plt.ylim(0,ymax)
    plt.yticks(np.arange(0,ymax + (ymax-0)/(yticks-1), (ymax-0)/(yticks-1) ))
    plt.xticks([0,3,6,9])
    #plt.savefig("../data/Result/nucleotide_trend/" + 'nucleotide_trend_scatter_'+'all'+ '.tif')    
    plt.show()
    """





    # scatterplots between GFP and YFP (Fig S20 B)
    '''
    a = arti.dataframe['LogMean']
    y = yfp.dataframe['LogMean']
    ind_a = ~np.isnan(a)
    ind_y = ~np.isnan(y)


    ind = np.logical_and(ind_a, ind_y)
    print("Total data used = ", len(a[ind]))

    result = stats.pearsonr(a[ind],y[ind])
    print("Pearson-r:", result[0])  
    print("p_value:", result[1])
    0
    result = stats.spearmanr(a[ind],y[ind])
    print("Spearman-r:", result[0])  
    print("p_value:", result[1])
    
    plt.figure(figsize=(6,6))
    plt.scatter(a,y, c = 'forestgreen',s = 6, alpha = 1)
    plt.xlim(0,3.2)
    plt.ylim(0,2.7)
    
    plt.xticks(np.arange(0,3.2+3.2/4, 3.2/4 ))
    plt.yticks(np.arange(0,2.7+2.7/3, 2.7/3 ))
    #plt.savefig("../data/Result/pairwise_scatterplots/" + 'x_arti_y_yfp - Spearmanr_{}.tif'.format(result[0]))
    plt.show()
    '''

    













    '''
