# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 01:45:23 2018

@author: JeremyJ
"""
import numpy as np
import itertools
import matplotlib.pylab as plt 
from decimal import Decimal
import time
import collections
import pandas as pd
import networkx as nx
import scipy.stats
from sklearn.metrics import r2_score

from mpl_toolkits.mplot3d import Axes3D


def ceil(num, digit):
    return  np.ceil(num* 10**digit)/10**digit

def floor(num, digit):
    return  np.floor(num* 10**digit)/10**digit   

def negative_to_zero(array):
    return  np.array(list(map(lambda x: (abs(x)+x)/2 , array )))

def nth_root(value, n_root): 
    root_value = 1/float(n_root)
    return round(Decimal(value) ** Decimal(root_value),3)
 
def minkowski_distance(x,y,p_value):
    '''
    minkowski_distance can be used to calculate 
        Euclidean distance with p_value = 2
        Manhattan distance with p_value = 1
        Chebyshev distance with p_value = ∞
    '''
    return nth_root(sum(pow(abs(a-b),p_value) for a,b in zip(x, y)),p_value)

def square_rooted(x):
    return round(np.sqrt(np.nansum([a*a for a in x])),3)

def cosine_similarity(x,y):
    numerator = np.nansum([a*b for a,b in zip(x,y)])
    denominator = square_rooted(x)*square_rooted(y)
    return round(numerator/float(denominator),3)

def correlation_coef(x,y , cor_type = 'corrcoef', show_p = False):
    '''
    cor_type can only be 'corrcoef' or 'spearmanr'
    '''
    bad = ~np.logical_or(np.isnan(x), np.isnan(y))
    _x = np.compress(bad, x)
    _y = np.compress(bad, y)
    if cor_type == 'corrcoef':
        cor2, p = scipy.stats.pearsonr(_x,_y)
        if show_p == True:
            return cor2,p
        else:
            return cor2
    elif cor_type == 'spearmanr':
        if show_p == True:
            return scipy.stats.spearmanr(_x,_y)
        else:
            return scipy.stats.spearmanr(_x,_y)[0]
    else:
        print('Please specify correct cor_type!!')

def neighbors(seq):
    '''
    Generate a list of 1-mutation neighbors
    i.e. neighbors('AAA')--> 27 neighbors
    '''
    base = ['A','U','C','G']
    neighbors = [seq[:i]+ch+seq[i+1:] for i in range(len(seq)) if seq[i] != 'N' for ch in base if ch != seq[i]]
    return neighbors

def mutational_set(seq, letters = 'AUGC'):
    '''
    Generate all mutants at the non-N loci
    i.e. mutational_set('ANA') --> 4^2 variants
    '''
    locs = [i for i in range(len(seq)) if seq[i] != 'N']
    #print(locs)
    this_word = [[char] for char in seq]
    for loc in locs:
        #orig_char = seq[loc]
        this_word[loc] = [l for l in letters]
    for poss in itertools.product(*this_word):
        yield ''.join(poss)
    
        
def hm (matrix,xticklabel,yticklabel, cm = 'jet', title ="None",xlabel='-',ylabel = '+' ,cbarlabel = '', vmax= 4.450926202822742, vmin=0.0):
    from matplotlib import cm as CM
    fig, ax = plt.subplots()
    cmp = CM.get_cmap(cm)
    cmp.set_bad('w') 
    
    mask =  np.triu(np.ones(matrix.shape[0]),1)[::-1]
    mm = np.ma.array(matrix, mask=mask) 
    print("Max value ->", np.amax(mm))
    im = ax.imshow(mm, interpolation="nearest",cmap = cmp, vmin=vmin, vmax= vmax)

    # We want to show all ticks...
    ax.set_xticks(np.arange(0,len(xticklabel)))
    ax.set_yticks(np.arange(0,len(yticklabel)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(xticklabel)
    ax.set_yticklabels(yticklabel)
    ax.invert_yaxis()
    
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")
    
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    '''
    # Loop over data dimensions and create text annotations.
    for i in range(len(ylabel)):
        for j in range(len(xlabel)):
            text = ax.text(j, i, matrix[i, j],
                           ha="center", va="center", color="w")
    '''
    ax.set_title(title)
    for label in ax.get_xticklabels():
        label.set_visible(False)
    for label in ax.get_yticklabels():
        label.set_visible(False)    
    for label in ax.get_xticklabels()[::9]:
        label.set_visible(True)
    for label in ax.get_yticklabels()[::9]:
        label.set_visible(True)    
    
    fig.tight_layout()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

# some distance measuremetns
    
def max_separation_of_nt(seq):
    tmp =[index for index in range(len(seq)) if seq[index] != 'N']
    return (tmp[-1]-tmp[0]+1)

'''
def std_separation_of_nt(seq):
    tmp =[index for index in range(len(seq)) if seq[index] != 'N']
    return np.std(tmp)
'''

# functions to generate the epistasis tables

def SubSeq(seq, order_reduction):
    '''
    List is the position of non-N nt
    generate non-N positions and calculate its order
    '''
    order = 0
    index = 0
    array = []
    for i in seq:
        if i !='N':
            array.append(index)
            order += 1
        index += 1
    # print(array)
    List = [p for p in itertools.combinations(array, order-order_reduction)]
    # generate suborder sequences
    tmp = []    
    for p in List:
        reference = ['N']*len(seq)
        for index in p:
            reference[index] = seq[index]
        tmp.append(''.join(reference))
    return tmp

def SeqMean(df, order):
    print("\norder =", order)
    
    tmp = df['SeqID']
    if pd.isnull(tmp[0]) != True:
        SeqLen = len(tmp[0])
    else: 
        SeqLen = len(tmp[~pd.isnull(tmp)].iloc[0]) # use pd.isnull instead of np.isnan
    tStart = time.time()
    
    SeqList = []
    ExpList = []
    C = np.array(range(0,len(df.SeqID)))
    Ntable = list(itertools.product('ACGU', repeat=order))
    TotalCount = len(list(itertools.combinations(range(0,SeqLen), r=order)))*(4**order)//100 + 1
    
    for i in itertools.combinations(range(0,SeqLen), r=order):
        
        index = np.zeros(len(df.SeqID)).astype(int)
        for j in range(0,order):
            index += C//4**(SeqLen-1-i[j])%4*4**(order-1-j)
        
        tmp = np.array([['N']*SeqLen]*len(Ntable))
        tmp[:,i] = Ntable
        
        for j in range(0,4**order):
            SeqList.extend([''.join(tmp[j])])
            ExpList.append(np.mean(df.LogMean[np.where(index == j)[0]]))
            if len(ExpList)%TotalCount == 0: 
                print('=',end ='')
    tEnd = time.time()
    print("\ntime =", tEnd-tStart)
    return dict(zip(SeqList, ExpList))


#######  !!!!!  ####### The result is different ????
def generate_epistasis_table( input_path, file, output_path = "" , MaxOrder = 9, threshold =25):
    #output_path = "C:/Users/JeremyJ/Documents/Epistasis/data/Epistasis Table/" 
    print('######'+file+'######')
          
    df = pd.read_csv(input_path + file)
    #print(df)
    #print(df['SeqID'][0], type(df['SeqID'][0]))
    SeqLen = len(df['SeqID'][0])
    if 'rep' in file:
        df[df.Counts < threshold] = np.nan
        
    MultiSeqMean = collections.defaultdict(dict)
    MultiSeqMean[0] = {'N'*SeqLen: np.mean(df.LogMean)} # calculate order_0   # np.mean okay with nan value?
    
    for k in range(1,MaxOrder+1): # calculate order_1~9
        MultiSeqMean[k] = SeqMean(df, k)
    
    ### =======================================
    for order in range(1,MaxOrder+1):
        print("\norder = "+str(order))
        tStart = time.time()
        for subOrder in range(0,order):
            for key in MultiSeqMean[order]:
                MultiSeqMean[order][key] = MultiSeqMean[order][key] - np.nansum([MultiSeqMean[subOrder][v] for v in SubSeq(key,order-subOrder)])
        tEnd = time.time()
        print("time =", tEnd-tStart)
        pd.DataFrame(list(MultiSeqMean[order].items())).to_csv(output_path + file[:-4] +'_order'+str(order)+'.csv', header = ['SeqID','LogMean'], index = False) # save file

def hamming_distance(seq1,seq2):
    assert len(seq1) == len(seq2), 'Length of seq1 != seq2!!'
    dis = np.count_nonzero(   ~(np.array(list(seq1)) == np.array(list(seq2))) )
    return dis

def is_neighbor(seq1,seq2):
    '''
    Return a boolean depends on if hamming distance == 1
    '''
    
    dis = len(seq1) - np.count_nonzero( np.array(list(seq1)) == np.array(list(seq2)))
    if dis == 1:
        return True
    else:
        return False

def plot_network(list_of_seq, save =False):
    '''
    Given a list of sequences, plot their 2D network
    To Do
    '''
    g  = nx.Graph()
    for seq in list_of_seq:
        g.add_node(seq)
    # O(n^2) algorithm
    
    for seq1,seq2 in list(itertools.combinations(list_of_seq,2)):
        if is_neighbor(seq1,seq2):
            g.add_edge(seq1,seq2)
    print('Adding nodes/edges completed!!')
    pos = nx.spring_layout(g)
    #label_dict = {i: ''.join([i[loc] for loc in locs]) for i in g.nodes()}
    #nx.draw_networkx_labels(g, pos)#, label_dict, font_size=10)
    
    nx.draw_networkx_nodes(G = g, pos = pos, node_list = g.nodes(),node_color = 'r', alpha = 0.8, node_size = 15)
    nx.draw_networkx_edges(G = g, pos = pos, edge_color='g', alpha = 1.0, arrows = False)
    #nx.draw(g,pos, with_labels = False,node_size=10)
    
    print("Transitivity(the fraction of all possible triangles present in G)", nx.transitivity(g))
    print("Average clustering coefficient" , nx.average_clustering(g))

    #import graphistry
    #graphistry.register('myapikey')
    #graphistry.bind(source='src', destination='dst', node='nodeid').plot(g)
    
    plt.show()
    
    tmp  = list(nx.find_cliques(g))
    #print(tmp)
    print("# of cliques:", len(tmp))

    #return g, pos #, label_dict    


def normal_distribution_test(data, types = 'QQplot', save = False, figname = ''):
    '''
    (1) Shapiro-Wilk (S-W) test:
        SW test focuses on the tails, which is where we typically do care if the distributions are similar. 
        SW is usually preferred.
        
    (2) Kolmogorov-Smirnov test (K-S): 
        KS test looks at the quantile where your empirical CDF differs 
        maximally from the normal's theoretical CDF
        
    (3 )Anderson-Darling test:
    (4) QQ plot: essential
    '''
    data = np.array(list(data))
    ok = ~np.isnan(data)
    data = data[ok]

    if types == 'SW-test':
        test_stat, p_value = scipy.stats.shapiro(data)
        return test_stat, p_value
    elif types == 'KS-test':
        test_stat, p_value = scipy.stats.kstest(data,'norm',)
        return test_stat, p_value
    elif types == 'Andreson':
        # ‘norm’,’expon’,’logistic’,’gumbel’,’gumbel_l’, gumbel_r’, ‘extreme1’
        test_stat, p_value,significant_list = scipy.stats.anderson(data,'norm') 
        return test_stat, p_value, significant_list
    
    elif types == 'QQplot':
        scipy.stats.probplot(data, dist = 'norm',plot=plt)
        if save == True:
            plt.savefig("QQplot"+ figname +'.svg')
        # (osm, osr), (slope, intercept, r) = 
    else:
        print("Wrong type!!! Please enter: QQplot, SW-test, KS-test or Andreson ")
    

def scatter_plot(x_array,y_array, save = False, filename = ''):
    x = np.array(list(x_array))
    y = np.array(list(y_array))
    ind_x = ~np.isnan(x)
    ind_y = ~np.isnan(y)
    ind = np.logical_and(ind_x, ind_y)
    x = x[ind]
    y = y[ind]
    
    # calculate r2 score
    r2 = r2_score(x,y)
    print("r-squared:", r2)
    
    # plot scatterplot 
    plt.scatter(x,y, s= 0.05)
    if save == True:
        plt.savefig("scatterplot_fitted-exp_"+filename+".svg")
    plt.show()


def scatter_plot3D(x_array,y_array, z_array, axis_min =0, axis_max = 3.2, point_size = 3 ,elev = 30,azim = -60 ,ticks = 5 , labelsize = 14,save = False, filename = '',
                   label1 = 'Replicate 1', label2 = 'Replicate 2', label3 = 'Replicate 3'):
    x = np.array(list(x_array))
    y = np.array(list(y_array))
    z = np.array(list(z_array))
    ind_x = ~np.isnan(x)
    ind_y = ~np.isnan(y)
    ind_z = ~np.isnan(z)
    ind = np.logical_and(np.logical_and(ind_x, ind_y),ind_z)
    x = x[ind]
    y = y[ind]
    z = z[ind]
    
    # plot scatterplot 
    fig = plt.figure(figsize=(6,6))
    ax = plt.axes(projection='3d')
    
    ax = Axes3D(fig)  
    ax.scatter3D(x, y, z, c='k',  s= point_size) # cmap='Greens'

    ax.set_xlabel(label1, fontsize=labelsize)
    ax.set_ylabel(label2, fontsize=labelsize)
    ax.set_zlabel(label3, fontsize=labelsize)

    ax.set_xlim(axis_min, axis_max)
    ax.set_ylim(axis_min, axis_max)
    ax.set_zlim(axis_min, axis_max)
    
    ax.set_xticks(np.arange(axis_min, axis_max+(axis_max-axis_min)/(ticks-1), (axis_max-axis_min)/(ticks-1)  ))
    ax.set_yticks(np.arange(axis_min, axis_max+(axis_max-axis_min)/(ticks-1), (axis_max-axis_min)/(ticks-1)  ))
    ax.set_zticks(np.arange(axis_min, axis_max+(axis_max-axis_min)/(ticks-1), (axis_max-axis_min)/(ticks-1)  ))
    ax.grid(False)
    print(ax.azim)
    print(ax.elev)
    print(ax.dist)
    ax.tick_params(axis = 'both', which = 'major', labelsize = labelsize)
    ax.view_init(elev=elev, azim=azim)
    # ax = plt.axes(projection='3d')
    # ax.scatter(x, y, z, c=z, cmap='viridis', linewidth=0.5); 
    #
    
    if save == True:
        plt.savefig("../data/Result/scatterplot_3D/" +"scatterplot_3D_{}_".format(azim)+filename+".tif")
    plt.show()


def table_unify():
    pass

def sequence_logo(ax,list_of_variants, filename = "", save = False, show_tick = False):
    if len(list_of_variants) > 0:
        matrix = np.zeros((4, len(list_of_variants[0])))
        tmp = np.stack([ list(variant) for variant in list_of_variants ]).T
        for i,j in itertools.product(range(4),range(9)):
            nt_list = ['A','U','C','G']
            matrix[i,j] = tmp[j][tmp[j] == nt_list[i]].shape[0]
            
    
        print(matrix)
        matrix_normed = matrix*100 / matrix.sum(axis=0)
        matrix_normed = matrix_normed.T
        N = len(list_of_variants[0])
        ind = np.arange(N)    # the x locations for the groups
        width = 1       # the width of the bars: can also be len(x) sequence
        #x_interval = (maxi-mini)/(ticks-1)
        #fig = plt.figure(figsize=(9, 4))
        
        ax.bar(ind, matrix_normed[:,0], width, color = 'g', edgecolor='white') #A
        ax.bar(ind, matrix_normed[:,1], width, bottom=matrix_normed[:,0], color = 'r',  edgecolor='white') #U
        ax.bar(ind, matrix_normed[:,2], width, bottom=(matrix_normed[:,1]+ matrix_normed[:,0]), color ='deepskyblue', edgecolor='white') #C
        ax.bar(ind, matrix_normed[:,3], width, bottom=(matrix_normed[:,2]+ matrix_normed[:,1]+matrix_normed[:,0]), color = 'k', edgecolor='white') #G
        
        #plt.xticks(r, names, fontweight='bold')
        
        # This line is problematic
        #tmp = [N/(ticks-1)*i for i in range(ticks)]
        #plt.xticks(tmp,[str(round(i,2)) for i in np.arange(mini,maxi+x_interval,x_interval)])
        #plt.xticks(ind, ['R'+str(i) for i in range(1,21)])
        '''
        if show_tick == True:
            plt.yticks(np.arange(0, 100+10, 10))
            plt.xticks(np.arange(0,9+1,1))
        else:
            plt.yticks([])
            plt.xticks([])
            
        if save == True:
            plt.savefig( "../data/Result/sequence_logo/"+filename+".svg" )
        '''
        #plt.tight_layout()
        #ax = plt.gca().yaxis.grid(True)
        #plt.title(self.name)
        
        #return fig, 
        ax.axis('off')
    else:
        ax.axis('off')
        #plt.yticks([])
        #plt.xticks([])
        print("No variants in the set.")
    return ax

def all_shortest_paths(string1, string2, plot = False):
    if len(string1) != len(string2):
        raise Exception('The length of two strings are not the same !!')
        
    difference_index = [i for i in range(len(string1)) if string1[i] != string2[i]]
    max_steps = len(difference_index)
    paths = list(itertools.permutations(difference_index, max_steps))
    
    my_array = [ ([string1]+['' for j in range(max_steps)]) for i in range(len(paths))]#np.empty([len(paths),],dtype='<U9')

    def one_step(s1,s2,pos):
        template = list(s1)
        template[pos] = s2[pos]
        return ''.join(template)
    
    for i in range(len(paths)):
        for index,pos in enumerate(paths[i]):
            my_array[i][index+1] = one_step(my_array[i][index],string2,pos)
    
    return my_array
    


def permutation_test_on_correlation(array1,array2, n_samples):
    '''
    This function may take a long time 
    '''
    
    # Compute observed correlation: r_obs
    r_obs = correlation_coef(array1,array2)
    
    # Initialize permutation replicates: perm_replicates
    perm_replicates = np.empty(n_samples)
    # Draw replicates
    for i in range(len(perm_replicates)):
        # Permute illiteracy measurments: illiteracy_permuted
        array1_permuted = np.random.permutation(array1)
        # Compute Pearson correlation
        perm_replicates[i] = correlation_coef(array1_permuted, array2)
    # Compute p-value: p
    p = np.sum(perm_replicates>r_obs)/len(perm_replicates)
    print('p-val =', p)    
    
