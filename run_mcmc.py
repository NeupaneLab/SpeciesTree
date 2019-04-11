import simulation_settings
import re
import pandas as pd
import numpy as np
import scipy as sp
from scipy.stats import norm
from scipy.stats import multivariate_normal as mvn
#from random import randrange, uniform

# from statistics import variance
import itertools
import random
from random import shuffle, uniform, randint
from math import exp, log, lgamma, sqrt
import re, os, itertools, sys, glob
import matplotlib.pyplot as plt
# import seaborn as sns
import itertools
from scipy.stats import expon
from scipy.stats import lognorm
import pyper as pr
from collections import Counter
r = pr.R()

# import all the settings/prior from simulation_settings module
master_seed = simulation_settings.master_seed
random.seed(master_seed)
n_gen = simulation_settings.n_gen
save_every = simulation_settings.save_every
mean_prior_alpha = simulation_settings.mean_prior_alpha
mean_prior_sd = simulation_settings.mean_prior_sd
mean_sigma = simulation_settings.mean_sigma
sig_sd = simulation_settings.sig_sd
start_means = simulation_settings.start_means
start_sigmas = simulation_settings.start_sigmas
sim_root_brow = simulation_settings.sim_root_brow
start_root_brow = simulation_settings.start_root_brow
brown_mean_prior_alpha = simulation_settings.brown_mean_prior_alpha
brown_sig_sd = simulation_settings.brown_sig_sd
rate_multiplier = simulation_settings.rate_multiplier
run_on_prior = simulation_settings.run_on_prior
start_set =  simulation_settings.start_set
random_start_tree = simulation_settings.random_start_tree
data_file = simulation_settings.data_file

np.random.seed(7765219)
shuffle_nu = 4000

if simulate_morpho_data:
    data_size = simulation_settings.data_size
    define_sets = simulation_settings.define_sets
    morpho_sigma = simulation_settings.morpho_sigma



# converting edge-lengths to VCV matrix
def getvcv(tre):    
    r.assign("text.string", tre) ##assinging tree string to R

    matr = r("""
    library(ape)

    T<-read.tree(text=text.string)

    matrix = vcv(T) 
    matrix
    """)

    output = r.get("matrix")
    return output

#executing getvcv function
sim_cov = getvcv(random_start_tree)* rate_multiplier

#generate starting means under brownian motion with mean = sim_root_brow and VCV = sim_cov
def simulate_mean(rate_multiplier):    
    rand_simulated_data = mvn.rvs(mean = sim_root_brow, cov=sim_cov, size = 1)
    return rand_simulated_data

rand_simulated_means = simulate_mean(rate_multiplier) # rate_multiplier scales VCV matrix with rate_multiplier =1 as all edge-lengths are exactly converted  to VCV 
print 'rand_simulated_means=', rand_simulated_means

#simulating data for each tip with means generated under brownian motion (by executing function simulate_mean), parameter with size = data_size, and sigmas as morpho_sigma
def simulate_data(n):
    morpho_data = {}
    big_list= []
    

    for i, set in enumerate(define_sets):
        for k in set:
            sim_morpho = np.random.normal(rand_simulated_means[i], morpho_sigma[i], n)
            
            big_list.append(list(sim_morpho))
            morpho_data[k]=list(sim_morpho)
    return morpho_data
morpho_data = simulate_data(data_size)

#defining different updating class objects for mcmc
class node(object):
    def __init__(self, ndnum):              # initialization function
        self.rsib = None                    # right sibling
        self.lchild = None                  # left child
        self.par = None                     # parent node
        self.number = ndnum                 # node number (internals negative, tips 0 or positive)
        self.taxname = ''                   # leaf taxa
        self.edgelen = 0.0                  # branch length
        self.descendants = set([ndnum])     # set containing descendant leaf set
        self.partial_mean = None                 # will have length 4*npatterns
        self.partial_sigma = None                 # will have length 4*npatterns
        self.partial_morphos = None
        self.current_ln_like = None
        self.confset = None
        self.current_sigma_ln_prior_density = None
        self.current_mean_ln_prior_density = None
        self.current_brow_mean_ln_prior_density = None
        self.current_ln_posterior = None
        self.current_brow_ln_like = None
        self.current_set_model_prior = None
        self.types = None
        self.current_topology = ''
        self.collapse_topology = ''
        self.split_topology = ''
        self.vcv = None
        self.current_ln_posterior = None
        self.brow_mean = 0.0
        self.delta1 = None
        self.lam1 = None

    def __str__(self):
        # __str__ is a built-in function that is used by #print to show an object
        descendants_as_string = ','.join(['%d' % d for d in self.descendants])

        lchildstr = 'None'
        if self.lchild is not None:
            lchildstr = '%d' % self.lchild.number
            
        rsibstr = 'None'
        if self.rsib is not None:
            rsibstr = '%d' % self.rsib.number
            
        parstr = 'None'
        if self.par is not None:
            parstr = '%d' % self.par.number
            
        return 'node: number=%d edgelen=%g lchild=%s rsib=%s parent=%s descendants=[%s]' % (self.number, self.edgelen, lchildstr, rsibstr, parstr, descendants_as_string)


#####read newick tree and postorder tree traversal
def readnewick(tree):
    total_length = len(tree)
    internal_node_number = -1

    root = node(internal_node_number)
    nd = root
    i = 0
    pre = [root]
    while i < total_length:
        m = tree[i]

        if m =='(':
            internal_node_number -= 1

            child = node(internal_node_number)
            pre.append(child)
            nd.lchild=child

            child.par=nd
            nd=child
        elif m == ',':
            internal_node_number -= 1
            rsib = node(internal_node_number)
            pre.append(rsib)
            nd.rsib = rsib
            rsib.par=nd.par
            nd = rsib
        elif m == ')':
            nd = nd.par
    
        elif m == ':':
            edge_len_str = ''
            i+=1
            m = tree[i]
            assert m in ['0','1','2','3','4','5','6','7','8', '9','.']
            while m in ['0','1','2','3','4','5','6','7','8', '9','.']:
                edge_len_str += m
                i+=1
                m = tree[i]
            i -=1
            nd.edgelen = float(edge_len_str)


        else:
            if not m==';':
                internal_node_number += 1

                if True:
    #                 assert m in ['0','1','2','3','4','5','6','7','8','9', 's', 'a', string.ascii_lowercase], 'Error : expecting m to be a digit when in fact it was "%s"' % m
                    mm = ''
                    while m in ['0','1','2','3','4','5','6','7','8', '9'] or re.match('[a-zA-Z_ .;]+',m):
                        mm += m
                
                        i += 1
                        m = tree[i]
                    
    #                 nd.number = int(mm)
                    ext_node_number = range(1, 1000)
                
                    nd.number = ext_node_number[i]
                    nd.taxname = mm
                    i -= 1
        
            else:
                pass
        i += 1

    post = pre[:]
    post.reverse()
#     for i in pre:
#         #print i.taxname or i.number,
    return post

postorder = readnewick(random_start_tree)
morphovalues = []

#convert taxa names to digits for easier tree manipulation
manu_topology = random_start_tree[:]
for nd in postorder:
    if nd.number > 0:
        manu_topology = manu_topology.replace(nd.taxname, str(nd.number))


def taxa_list():
    taxa_list = {}
    for n,i in enumerate(postorder):

        if True:

            if i.taxname in morpho_data:
                taxa_list[i.number] = i.taxname
            else:
                pass
    return taxa_list
taxa_list = taxa_list()

preorder = postorder[:]
preorder.reverse()

def mega_list():
    mg_list = []

    for i,nd in enumerate(preorder):
        ##print 'starting new new node..............................................'   
    
        y = 1
        nodes = preorder[i:]
        length = len(nodes)
        newlist = []
        newlist.append(preorder[i].number)
        ##print 'root node.....', newlist
    
        while y < length:
            ##print 'possible children====', nodes[y].number,
            if nodes[y].par.number in newlist:
                newlist.append(nodes[y].number)
            y +=1
            ##print
        ##print 'newlist===', newlist
        newlist1 = [item for item in newlist if item >=0]
        ##print 'only tip newlist ==', newlist1
        mg_list.append(newlist1)
    return mg_list
mega_list = mega_list()

 
 

    
print '*************************************internal and external Nodes ********************************************************************'    
print mega_list
print '*************************************internal and external Nodes  ********************************************************************'    
print

#likelihood function
def ln_likelihood(data, me, si):
#     sigma = np.std(data)
    lh = norm(me, si).logpdf(data)
    return lh.sum()

# getting data values contained by each node (internal and external) during mcmc move between sets for likelihood calculation
def node_morphs():
    node_morphs = {}
    for i, each in enumerate(mega_list):
        morphs = []
        for j in each:
            taxa_name = taxa_list[j]
            morphs.extend(morpho_data[taxa_name])
            node_morphs[preorder[i].number] = morphs
    return node_morphs
node_morphs = node_morphs()

# getting tips contained by internal nodes during mcmc move between sets
def contained_tips():
    node_tips = {}
    for i, each in enumerate(mega_list):
        tip = []
        for j in each:
            taxa_name = taxa_list[j]
            tip.extend([j])
            node_tips[preorder[i].number] = tip
    return node_tips
contained_tips = contained_tips()

print 'contained_tips===', contained_tips


# Assigning data, means, and sigmas to each node during mcmc move
def node_values():
    externalnode_morphos = {}
    for key,value in taxa_list.items():
        externalnode_morphos[key] = morpho_data[value]
    node_morphs = {}
    for l,m in enumerate(postorder):
        node_morphs[m.number] = []

    node_morphs.update(externalnode_morphos) 
    
    node_means = {}
    h = 0
    for l,m in enumerate(preorder):
        if m.number >0:
            node_means[m.number] = start_means[h]
            h+=1 
        else:
            node_means[m.number] = 0.0 
       
    node_sigmas = {}
    h = 0
    for l,m in enumerate(preorder):
        if m.number >0:
            node_sigmas[m.number] = start_sigmas[h]
            h+=1 
        else:
             node_sigmas[m.number] = 0.0
  
    
    return node_morphs, node_means, node_sigmas, externalnode_morphos

# computing log-likelihood for each member of the set
def cal_like(each, mea,sig):
    loglikelihood_each_combo = 0
    morphs = []
    for j in each:
        taxa_name = taxa_list[j]
#         print 'taxa_name===========>', taxa_name
        morphs.extend(morpho_data[taxa_name])
#     #print 'each, mean, sig=======>', each, np.mean(morphs), np.std(morphs)
    loglikelihood_inner_set = ln_likelihood(morphs, mea, sig)
    
    loglikelihood_each_combo += loglikelihood_inner_set
    return loglikelihood_each_combo

#computing mean prior log-density for each set
def cal_prior(mean, value):
    current_mean_ln_prior_density = -value/mean - log(mean)
    return current_mean_ln_prior_density

# function to check the nodes that can be collapsed during mcmc move
def num_nodes_can_combine(search_set):
    nodes_in_set = []
    for nodenumber,search in contained_tips.iteritems():
        for i in search_set:
            if search == i:
                nodes_in_set.append(nodenumber)
    
    nodescontained = []
    for l in nodes_in_set:
        for t in postorder: 
            if t.number == l:
               nodescontained.append(t) 
    dupl = []
#     print 'len(nodescontained)=', len(nodescontained), nodescontained
#     print '****__'
    for j in nodescontained:
        dupl.append(j.par.number)
#     print '____****'
    
#     print 'dupl===', dupl
    num_of_combinablesets=[i for i in dupl if dupl.count(i)>1]
    return len(num_of_combinablesets)/2.
    
# getting assigned node ID (in digits) for each member of the set
def grabNode(dicto,search_set):
    for nodenumber,search in dicto.iteritems():
        if search == search_set:
            return nodenumber

# computing joint log-likelihood for the entire set
def jointLnLike(confset, means, sigmas):
    total_loglikelihood = 0
    for m,each in enumerate(confset):
        nod = grabNode(contained_tips,each)
        me =  means.get(nod)
        sig = sigmas.get(nod)
        current_ln_like = cal_like(each, me, sig)
#         print 'current lnlike===>', current_ln_like
        total_loglikelihood+=current_ln_like
#     print 'Total lnlike===>', total_loglikelihood
    return total_loglikelihood

# computing mean prior log-density for each set
def jointMeanLnprior(confset, means):
    total_mean_logprior = 0
    print 'confset=', confset
    for m,each in enumerate(confset):
        nod = grabNode(contained_tips,each)
        me =  means.get(nod)

        current_mean_logprior = norm.logpdf(me, loc = mean_prior_alpha, scale = mean_prior_sd)
        total_mean_logprior += current_mean_logprior
    return total_mean_logprior

# computing mean sigma (st.dev.) log-density for each set
def jointSigmaLnprior(confset, sigmas):
    total_sigma_logprior = 0
    for m,each in enumerate(confset):
        nod = grabNode(contained_tips,each)
        sig = sigmas.get(nod)
        current_sigma_logprior = lognorm.logpdf(sig,scale=mean_sigma, s=sig_sd)
        
        
        total_sigma_logprior += current_sigma_logprior
    return total_sigma_logprior

# computing log-density for brownian/root mean  
def joint_brow_mean(browian_mean):    
    current_brow_mean_logprior = norm.logpdf(browian_mean, loc = brown_mean_prior_alpha, scale = brown_sig_sd)
    return current_brow_mean_logprior
    
#collapse or merge topology during mcmc 
def mantree(tree2, taxal, taxar, par):
    m1 = re.search('\(%s\:([0-9]*[.]?[0-9]+)'%(taxal),tree2)
    m2 = re.search('\,%s\:([0-9]*[.]?[0-9]+)\)\:([0-9]*[.]?[0-9]+)'%(taxar),tree2)

    texttoremove =  m1.group(0)+m2.group(0)
    newedge = float(m2.group(1))+ float(m2.group(2))  
    findstring =  texttoremove
    replace_with = str(par)+':'+str(newedge)
    
    tree3 = tree2[:]
    tree4 = tree3.replace(findstring, replace_with)
#     print tree2
#     print tree4
    return tree4

# getting preorder sequence from node class
def getPreorder(nd, start = True): 
    blen= 0

    global _preorder
    if start:
        _preorder = [] 

    _preorder.append(nd)
    
    if nd.lchild:
        getPreorder(nd.lchild, False)  

    for ib in _preorder[1:]:
        blen+=ib.edgelen

    return blen


# 
def splitTree(t2,pick):
    addedge = getPreorder(pick, start = True)
    no = str(pick.number)
    m = re.search('%s\:[0-9]*[.]?[0-9]+'%(no),t2)
    if pick.number == -1:
        replace_with = '('+str(pick.lchild.number)+':'+ str(addedge)+','+ str(pick.lchild.rsib.number)+':'+ str(addedge)+')'
    else:
        replace_with = '('+str(pick.lchild.number)+':'+ str(addedge)+','+ str(pick.lchild.rsib.number)+':'+ str(addedge)+')'+':'+str(pick.edgelen)

    findstring = m.group(0)
    t3 = t2[:]
    t4 = t3.replace(findstring, replace_with)
    return t4


    
def collapseTree(tree2, pick):
    m1 = re.search('\(%s\:([0-9]*[.]?[0-9]+)'%(pick.par.lchild.number),tree2)
    if pick.par.number == -1:
        m2 = re.search('\,%s\:([0-9]*[.]?[0-9]+)\)'%(pick.par.lchild.rsib.number),tree2)
    
    else:
        m2 = re.search('\,%s\:([0-9]*[.]?[0-9]+)\)\:([0-9]*[.]?[0-9]+)'%(pick.par.lchild.rsib.number),tree2)

    texttoremove =  m1.group(0)+m2.group(0)
    if pick.par.number == -1:
        newedge = float(m2.group(1))  
    
    else:
        newedge = float(m2.group(1))+ float(m2.group(2))  
    findstring =  texttoremove
    replace_with = str(pick.par.number)+':'+str(newedge)

    tree3 = tree2[:]
    tree4 = tree3.replace(findstring, replace_with)
    return tree4



def brow_jointLnLike(confset, means, topology, browia_mean ):
    num_lineages = len(confset)
    means_for_brow = [browia_mean] * num_lineages
    total_loglikelihood = 0
    data_means = []
    data_sigmas = []
    for m,each in enumerate(confset):
        nod = grabNode(contained_tips,each)
        me =  means.get(nod)
        data_means.append(me)
    vcv = getvcv(topology)*rate_multiplier
    lgdensity_on_brow = mvn.logpdf(data_means, mean=means_for_brow, cov=vcv)
    return lgdensity_on_brow


def moves(pick, movetype):
    if movetype == 'mean':
        print '_MOVE MEAN_______________________________________________'

        print 'node.partial_mean before change', node.partial_means

        mean = node.partial_means.get(pick.number)
        sigma = node.partial_sigmas.get(pick.number)
        u = random.random()
        node.delta1 = uniform(0,10)
        node.lam1 = 0
        mean_new = mean-node.delta1/2.+node.delta1*u
        
        node.partial_means[pick.number]= mean_new
        print 'node.partial_mean after change', node.partial_means
        
        proposed_mean_ln_prior_density = jointMeanLnprior(node.confset, node.partial_means)
        
        if run_on_prior:
            proposed_ln_posterior =    proposed_mean_ln_prior_density + node.current_sigma_ln_prior_density + node.current_brow_mean_ln_prior_density
        else:
#             print '/////////////////////////START'
            proposed_ln_like = jointLnLike(node.confset, node.partial_means,node.partial_sigmas)

            if pick.number == -1:
                me = node.brow_mean
                si = sqrt(3.6*rate_multiplier)
                data = node.partial_means.get(pick.number)
            
                proposed_brow_ln_like = norm(me, si).logpdf(data) 
            
            else:
                node.vcv = getvcv(node.current_topology)
                print node.vcv
                proposed_brow_ln_like = brow_jointLnLike(node.confset, node.partial_means, node.current_topology, node.brow_mean)
#             print '/////////////////////////END'
        
            proposed_ln_posterior =  proposed_ln_like + proposed_brow_ln_like +  proposed_mean_ln_prior_density + node.current_sigma_ln_prior_density + node.current_brow_mean_ln_prior_density

#         hastings_ratio = log(m)
#         logR = proposed_ln_posterior - node.current_ln_posterior + hastings_ratio
        logR = proposed_ln_posterior - node.current_ln_posterior 
#         print 'ln_like before and after mean move =========......>>>', node.current_ln_like, proposed_ln_like
        
        print 'proposed_ln_posterior, current_ln_posterior', proposed_ln_posterior, node.current_ln_posterior
        print 'proposed_mean_ln_prior_density node.current_sigma_ln_prior_density + node.current_brow_mean_ln_prior_density =', proposed_mean_ln_prior_density, node.current_sigma_ln_prior_density + node.current_brow_mean_ln_prior_density

        u = random.random()
        print 'log(u) logR ==',  log(u), logR
        print '              **********************************'
        if log(u) < logR:
            node.types = 'Accept'
            node.partial_means[pick.number] = mean_new
            if run_on_prior is False:
                node.current_brow_ln_like = proposed_brow_ln_like
                node.current_ln_like = proposed_ln_like
            node.current_mean_ln_prior_density = proposed_mean_ln_prior_density
            node.current_ln_posterior = proposed_ln_posterior
            print 'ACCEPT, log(u) < logR, new proposal...mean =', mean_new

        else:
            node.types = 'Reject'
            node.partial_means[pick.number] = mean
            print 'REJECT, log(u) > logR, new proposal...mean =', mean
        print '              **********************************'
 
        print       

    elif movetype == 'brown_mean':
        print '_MOVE BROWN MEAN_______________________________________________'

        brow_mean = node.brow_mean
        u = random.random()
        node.delta1 = uniform(0,10)
        node.lam1 = 0
        brow_mean_new = brow_mean-node.delta1/2.+node.delta1*u
        proposed_brow_mean_ln_prior_density =  joint_brow_mean(brow_mean_new)
        
        print 'brow_mean, brow_mean_new==', brow_mean, brow_mean_new
        
        if run_on_prior:
            proposed_ln_posterior =  node.current_mean_ln_prior_density + node.current_sigma_ln_prior_density + proposed_brow_mean_ln_prior_density
        else:
            proposed_ln_like = jointLnLike(node.confset, node.partial_means,node.partial_sigmas)

#             print '/////////////////////////START'
            if pick.number == -1:
                me = node.brow_mean
                si = sqrt(3.6*rate_multiplier)
                data = node.partial_means.get(pick.number)
            
                proposed_brow_ln_like = norm(me, si).logpdf(data) 
            
#                 print 'me, si, data, proposed_brow_ln_like=', me, si, data, proposed_brow_ln_like 

            else:
                node.vcv = getvcv(node.current_topology)
                print node.vcv
                proposed_brow_ln_like = brow_jointLnLike(node.confset, node.partial_means, node.current_topology, brow_mean_new)
#                 print 'proposed_brow_ln_like ==', proposed_brow_ln_like
#             print '/////////////////////////END'

            proposed_ln_posterior = proposed_ln_like + proposed_brow_ln_like + node.current_mean_ln_prior_density + node.current_sigma_ln_prior_density + proposed_brow_mean_ln_prior_density
        
        logR = proposed_ln_posterior - node.current_ln_posterior
        print 'proposed_ln_posterior, current_ln_posterior', proposed_ln_posterior, node.current_ln_posterior
        print 'node.current_mean_ln_prior_density,  node.current_sigma_ln_prior_density, proposed_brow_mean_ln_prior_density =', node.current_mean_ln_prior_density,  node.current_sigma_ln_prior_density, proposed_brow_mean_ln_prior_density
        
        u = random.random()

        print 'log(u) logR ==',  log(u), logR
        print '              **********************************'
        if log(u) < logR:
            node.types = 'Accept'
            if run_on_prior is False:
                node.current_brow_ln_like = proposed_brow_ln_like
                node.current_ln_like = proposed_ln_like
            node.brow_mean = brow_mean_new
            node.current_ln_posterior = proposed_ln_posterior
            node.current_brow_mean_ln_prior_density = proposed_brow_mean_ln_prior_density
            print 'ACCEPT, log(u) < logR, new proposal...brow_mean =', brow_mean_new

        else:
            node.types = 'Reject'
            print 'REJECT, log(u) > logR, current...brow_mean =', node.brow_mean
        print '              **********************************'
 
    
    elif movetype == 'sigma':
        print '_MOVE SIGMA_______________________________________________'

        mean = node.partial_means.get(pick.number)
        sigma = node.partial_sigmas.get(pick.number)

        print 'node.partial_sigmas before change', node.partial_sigmas

        node.lam1 = uniform(0,10)
        node.delta1 = 0
        
        u = random.random()
        m = exp(node.lam1*(u-0.5))
        sigma_new = sigma * m
#         sigma_new = sigma-lam1/2.+lam1*u
#         if sigma_new < 0.0:
#             sigma_new = -sigma_new

        node.partial_sigmas[pick.number]= sigma_new
        print 'node.partial_sigmas after change', node.partial_sigmas

        proposed_sigma_ln_prior_density = jointSigmaLnprior(node.confset, node.partial_sigmas)

        if run_on_prior:
            proposed_ln_posterior =  node.current_mean_ln_prior_density + proposed_sigma_ln_prior_density + node.current_brow_mean_ln_prior_density
        else:
        
            proposed_ln_like = jointLnLike(node.confset, node.partial_means,node.partial_sigmas)

#             print '/////////////////////////START'
            if pick.number == -1:
                me = node.brow_mean
                si = sqrt(3.6*rate_multiplier)
                data = node.partial_means.get(pick.number)
            
                proposed_brow_ln_like = norm(me, si).logpdf(data) 
#                 print 'me, si, data, proposed_brow_ln_like=', me, si, data, proposed_brow_ln_like 
            else:
                node.vcv = getvcv(node.current_topology)
                print node.vcv
                proposed_brow_ln_like = brow_jointLnLike(node.confset, node.partial_means, node.current_topology, node.brow_mean)
#                 print 'proposed_brow_ln_like ==', proposed_brow_ln_like
#                 print '/////////////////////////END'
            proposed_ln_posterior = proposed_ln_like +  proposed_brow_ln_like + node.current_mean_ln_prior_density + proposed_sigma_ln_prior_density + node.current_brow_mean_ln_prior_density

        hastings_ratio = log(m)
        logR = proposed_ln_posterior - node.current_ln_posterior + hastings_ratio

        print 'proposed_ln_posterior, current_ln_posterior', proposed_ln_posterior, node.current_ln_posterior
        print 'node.current_mean_ln_prior_density , proposed_sigma_ln_prior_density , node.current_brow_mean_ln_prior_density=', node.current_mean_ln_prior_density , proposed_sigma_ln_prior_density , node.current_brow_mean_ln_prior_density
        u = random.random()
        print 'log(u), logR,  ==', log(u), logR
        print '              **********************************'
        if log(u) < logR:
            node.types = 'Accept'
            node.partial_sigmas[pick.number] = sigma_new
            if run_on_prior is False:
                node.current_brow_ln_like = proposed_brow_ln_like
                node.current_ln_like = proposed_ln_like
            node.current_sigma_ln_prior_density = proposed_sigma_ln_prior_density
            node.current_ln_posterior = proposed_ln_posterior
            print 'ACCEPT, log(u) < logR, new proposal...sigma =', sigma_new

        else:
            node.types = 'Reject'
            node.partial_sigmas[pick.number] = sigma
            print 'REJECT, log(u) > logR, new proposal...sigma =', sigma
        print '***************************************************'
        print

    elif movetype == 'movedown':
        print '_MOVE DOWN_______________________________________________'
        node.delta1 = uniform(0,10)
        node.lam1 = uniform(0,10)
        
        
        topology = node.current_topology[:]
        node.collapse_topology = collapseTree(topology, pick)
        


        prob_movedown = 0
        prob_moveup = 0

        
        if pick.number > 1:
            prob_movedown += 1
            
        else:
            prob_movedown += 1/2.
        
        if pick.par.number == -1:
            prob_moveup += 1
        
        else:
            prob_moveup += 1/2.
# 
        print' node.partial_means before collapsed =', node.partial_means
        print' node.partial_sigmas before collapsed=', node.partial_sigmas


        new_confset = node.confset[:]
        nodeadd = pick.par.number
        nodel = contained_tips[pick.par.lchild.number]
        nodelrisb = contained_tips[pick.par.lchild.rsib.number]
        print '>>>>>>>>>>>>>>>>>>>>>>>movedown'
    
        print 'confset====', node.confset
        new_confset[new_confset.index(nodel)] = contained_tips[nodeadd]
        new_confset.remove(nodelrisb)
        print 'new_confset =',  new_confset
        
        mean_old = node.partial_means[pick.par.number]
        sigma_old = node.partial_sigmas[pick.par.number]
        
        mean_cherries = [node.partial_means.get(pick.par.lchild.number), node.partial_means.get(pick.par.lchild.rsib.number)]
        sigma_cherries = [node.partial_sigmas.get(pick.par.lchild.number), node.partial_sigmas.get(pick.par.lchild.rsib.number)]
        
        mean_new_chosen = random.choice([node.partial_means.get(pick.par.lchild.number), node.partial_means.get(pick.par.lchild.rsib.number)])
        sigma_new_chosen = random.choice([node.partial_sigmas.get(pick.par.lchild.number), node.partial_sigmas.get(pick.par.lchild.rsib.number)])



        v = random.random()
        mean_new1 = mean_new_chosen - node.delta1/2.+ node.delta1*v


        w = random.random()
        vv2 = sigma_new_chosen - node.lam1/2.+ node.lam1*w
        if vv2 < 0.0:
            vv2 = -vv2
        sigma_new1 = vv2

        
        node.partial_means[pick.par.number]= mean_new1
        node.partial_sigmas[pick.par.number]= sigma_new1
        print' node.partial_means after collapsed =', node.partial_means
        print' node.partial_sigmas after collapsed=', node.partial_sigmas

        n_c = num_nodes_can_combine(node.confset)
        n_s = 0.0
        for m in node.confset:
            if len(m)>1:
                n_s+=1

        if pick.par.number ==-1:
            n_p_c = 0
        else:
            n_p_c = num_nodes_can_combine(new_confset)
        
        n_p_s = 0.0
        for m in new_confset:
            if len(m)>1:
                n_p_s+=1

 
        print 'topology=',  topology
        print 'node.collapse_topology=',  node.collapse_topology

        proposed_mean_ln_prior_density = jointMeanLnprior(new_confset, node.partial_means)
        proposed_sigma_ln_prior_density = jointSigmaLnprior(new_confset, node.partial_sigmas)
 
        log_proposalRatio = (log(n_c+n_s))-(log(n_p_c + n_p_s)+log(node.delta1)+log(node.lam1))
        log_Jacobian = -log(node.delta1)-log(node.lam1)
#         print 'log_Jacobian==', log_Jacobian

        if run_on_prior:
            proposed_ln_posterior =   proposed_mean_ln_prior_density + proposed_sigma_ln_prior_density + node.current_brow_mean_ln_prior_density
        else:
            proposed_ln_like = jointLnLike(new_confset, node.partial_means,node.partial_sigmas)
#             print 'proposed_ln_like==', proposed_ln_like
#             print '/////////////////////////START'
            if pick.par.number == -1:
                me = node.brow_mean
                si = sqrt(3.6*rate_multiplier)
                data = node.partial_means.get(pick.par.number)
            
                proposed_brow_ln_like = norm(me, si).logpdf(data) 
#                 print 'me, si, data, proposed_brow_ln_like=', me, si, data, proposed_brow_ln_like 
            else:
                node.vcv = getvcv(node.collapse_topology)
                print node.vcv
                proposed_brow_ln_like = brow_jointLnLike(new_confset, node.partial_means, node.collapse_topology, node.brow_mean)
#                 print 'proposed_brow_ln_like ==', proposed_brow_ln_like
#                 print '/////////////////////////END'

        
            proposed_ln_posterior =  proposed_ln_like + proposed_brow_ln_like + proposed_mean_ln_prior_density + proposed_sigma_ln_prior_density + node.current_brow_mean_ln_prior_density

        logR = proposed_ln_posterior - node.current_ln_posterior + log_proposalRatio + log_Jacobian
        print 'proposed_ln_posterior, node.current_ln_posterior=', proposed_ln_posterior, node.current_ln_posterior

        print 'proposed_mean_ln_prior_density, proposed_sigma_ln_prior_density, node.current_brow_mean_ln_prior_density=', proposed_mean_ln_prior_density, proposed_sigma_ln_prior_density, node.current_brow_mean_ln_prior_density
        print 'log_proposalRatio, log_Jacobian, logR=', log_proposalRatio, log_Jacobian, logR

        u = random.random()
        print 'log(u) logR ==',  log(u), logR
        print '              **********************************'
        if log(u) < logR:
            node.types = 'Accept'
            node.confset = new_confset
            node.partial_means[pick.par.number] = mean_new1
            node.partial_sigmas[pick.par.number] = sigma_new1
            if run_on_prior is False:
                node.current_brow_ln_like = proposed_brow_ln_like
                node.current_ln_like = proposed_ln_like
            node.current_ln_posterior = proposed_ln_posterior
            node.current_topology = node.collapse_topology
            print 'ACCEPT, log(u) < logR, new proposal...mean =', new_confset, mean_new1

        else:
            node.types = 'Reject'
            node.partial_means[pick.par.number]= mean_old
            node.partial_sigmas[pick.par.number]= sigma_old
#             node.partial_means[pick.par.number] = mean
#             node.partial_sigmas[pick.par.number] = sigma
            print 'REJECT, log(u) > logR, new proposal...mean =', node.confset
        print '              **********************************'

        print

 
    elif movetype == 'moveup':
    
        print '_MOVE UP_______________________________________________', node.confset

        node.delta1 = uniform(0,10)
        node.lam1 = uniform(0,10)
    
        topology = node.current_topology[:]
        node.split_topology = splitTree(topology, pick)
        prob_movedown = 0
        prob_moveup = 0
        
        if pick.number == -1:
            prob_moveup +=1
            
        else:
            prob_moveup += 1/2.
        
        if pick.lchild.number > 0:
            prob_movedown += 1
        
        else:
            prob_movedown += 1/2.
    
    
        lchild = contained_tips.get(pick.lchild.number)
        rsib = contained_tips.get(pick.lchild.rsib.number)
        current_node = contained_tips.get(pick.number)
        tmp_confset = node.confset[:]
        new_confset =  tmp_confset[:tmp_confset.index(current_node)]+[lchild]+[rsib]+node.confset[node.confset.index(current_node):]
        new_confset.remove(current_node)
        print 'new_confset==', new_confset
 
        print 'node.partial_means before split=', node.partial_means
        print 'node.partial_sigmas before split=', node.partial_sigmas
 
 
       
        if pick.number == -1:
            n_c = 0
        else:
            n_c = num_nodes_can_combine(node.confset)
        
        n_s = 0.0
        for m in node.confset:
            if len(m)>1:
                n_s+=1


        n_p_c = num_nodes_can_combine(new_confset)
        
        n_p_s = 0.0
        for m in new_confset:
            if len(m)>1:
                n_p_s+=1
        
        current_mean = node.partial_means.get(pick.number)
        current_sigma = node.partial_sigmas.get(pick.number)
        v = random.random()
        child_means = [current_mean, current_mean - node.delta1/2.+ node.delta1*v]
        print 'child_means=', child_means


        w = random.random()
        vv = current_sigma - node.lam1/2.+ node.lam1*w
#         print 'vv==', vv
        if vv < 0.0:
            vv = -vv
        child_sigmas = [current_sigma, vv]
        print 'child_sigma=', child_sigmas

#         print '_____MOVE UP_____num_nodes_can_split, num_node_sets_can_combine ', num_nodes_can_split, num_node_sets_can_combine 
        
        left_right_means = [node.partial_means.get(pick.lchild.number), node.partial_means.get(pick.lchild.rsib.number)]
        left_right_sigmas = [node.partial_sigmas.get(pick.lchild.number), node.partial_sigmas.get(pick.lchild.rsib.number)]
        
        lrnode_list = [pick.lchild.number, pick.lchild.rsib.number]
        ranchoice = random.choice(lrnode_list)
        
        node.partial_means[ranchoice] = child_means[0]
        node.partial_sigmas[ranchoice] = child_sigmas[0]
        
        lrnode_list.remove(ranchoice)

        node.partial_means[lrnode_list[0]] = child_means[1]
        node.partial_sigmas[lrnode_list[0]] = child_sigmas[1]
        
        print 'node.partial_means after split=', node.partial_means
        print 'node.partial_sigmas split=', node.partial_sigmas


        print topology
        print node.split_topology
 
        proposed_mean_ln_prior_density = jointMeanLnprior(new_confset, node.partial_means)

        proposed_sigma_ln_prior_density = jointSigmaLnprior(new_confset, node.partial_sigmas)
        
        log_proposalRatio = (log(n_c+ n_s) + log(node.delta1)+log(node.lam1)) - (log(n_p_c+n_p_s))
        log_Jacobian = log(node.delta1)+log(node.lam1)

        if run_on_prior:
            proposed_ln_posterior =  proposed_mean_ln_prior_density + proposed_sigma_ln_prior_density  + node.current_brow_mean_ln_prior_density
        else:
            ######work below this######
            proposed_ln_like = jointLnLike(new_confset, node.partial_means,node.partial_sigmas)
        
#             print '/////////////////////////START'
            node.vcv = getvcv(node.collapse_topology)
            print node.vcv
            proposed_brow_ln_like = brow_jointLnLike(new_confset, node.partial_means, node.split_topology, node.brow_mean)
#             print 'proposed_brow_ln_like==', proposed_brow_ln_like
#             print '/////////////////////////END'
       
            proposed_ln_posterior =  proposed_ln_like + proposed_brow_ln_like + proposed_mean_ln_prior_density + proposed_sigma_ln_prior_density  + node.current_brow_mean_ln_prior_density

        logR = proposed_ln_posterior - node.current_ln_posterior + log_proposalRatio + log_Jacobian
        print 'proposed_ln_posterior, current_ln_posterior', proposed_ln_posterior, node.current_ln_posterior
        print 'proposed_mean_ln_prior_density, proposed_sigma_ln_prior_density, node.current_brow_mean_ln_prior_density=', proposed_mean_ln_prior_density, proposed_sigma_ln_prior_density, node.current_brow_mean_ln_prior_density
        print 'log_proposalRatio, log_Jacobian, logR=', log_proposalRatio, log_Jacobian, logR
        u = random.random()
        print 'log(u) logR ==',  log(u), logR
        print '              **********************************'
        if log(u) < logR:
            node.types = 'Accept'
            node.confset = new_confset
            if run_on_prior is False:
                node.current_brow_ln_like = proposed_brow_ln_like
                node.current_ln_like = proposed_ln_like
            node.current_ln_posterior = proposed_ln_posterior
            node.current_topology = node.split_topology
            print 'ACCEPT, log(u) < logR, new proposal...mean =', new_confset

        else:
            node.types = 'Reject'
            node.partial_means[pick.lchild.number]= left_right_means[0]
            node.partial_means[pick.lchild.rsib.number]= left_right_means[1]
            node.partial_sigmas[pick.lchild.number]= left_right_sigmas[0]
            node.partial_sigmas[pick.lchild.rsib.number]= left_right_sigmas[1]            
            print 'REJECT, log(u) > logR, new proposal...mean =', node.confset
        print '              **********************************'

        print

def exploreNodes(mega_list, n_gen):
#    old_stdout = sys.stdout
#    log_file = open("message.log.txt","w")
#    sys.stdout = log_file

    start_morphs = []
    node.partial_morphos = node_values()[0]
    #print 'node.partial_morphos==', node.partial_morphos
    node.partial_means = node_values()[1]
    node.partial_sigmas = node_values()[2]
    print' node.partial_means=', node.partial_means
    print' node.partial_sigmas=', node.partial_sigmas
    node.confset = start_set
#     node.confset =[[19, 43], [74, 91], [120, 140]]

    node.current_ln_like = jointLnLike(node.confset, node.partial_means, node.partial_sigmas)
    print 'curren_ln_lik===', node.current_ln_like

#     node.current_set_model_prior = log(2)-log(len(node.confset))
    #print 'current_set_model_prior========', node.current_set_model_prior

    node.current_mean_ln_prior_density = jointMeanLnprior(node.confset, node.partial_means)
    print 'current_mean_ln_prior_density========', node.current_mean_ln_prior_density

    node.current_sigma_ln_prior_density = jointSigmaLnprior(node.confset, node.partial_sigmas)
    print 'current_sigma_ln_prior_density========', node.current_sigma_ln_prior_density
    

    print '*********%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print '*********%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

    node.brow_mean = start_root_brow[0]
    node.current_brow_mean_ln_prior_density =  joint_brow_mean(node.brow_mean)
    print 'node.brow_mean ==', node.brow_mean
    print 'node.current_brow_mean_ln_prior_density==', node.current_brow_mean_ln_prior_density
    node.current_topology = manu_topology[:]
    node.vcv = getvcv(node.current_topology)
    print node.vcv
    node.current_brow_ln_like = brow_jointLnLike(node.confset, node.partial_means, node.current_topology, node.brow_mean)
    print 'node.brow_jointLnLike==', node.current_brow_ln_like

    print '*********%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    print '*********%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'


    if run_on_prior:
        node.current_ln_posterior  = node.current_mean_ln_prior_density + node.current_sigma_ln_prior_density + node.current_brow_mean_ln_prior_density
    else:
        node.current_ln_posterior  = node.current_ln_like + node.current_mean_ln_prior_density + node.current_sigma_ln_prior_density + node.current_brow_mean_ln_prior_density
    print ' node.current_ln_posterior=',  node.current_ln_posterior
    print ' node.current_mean_ln_prior_density node.current_sigma_ln_prior_density =',  node.current_mean_ln_prior_density, node.current_sigma_ln_prior_density

    node.types = 'Start'
    #print node.types

    mcmc = 0
    output = os.path.join('mcmc.txt')
    newf = open(output, 'w')
    newf.write('%s\t'%('n_gen'))

    for each, h in node.partial_means.iteritems():
        newf.write('%d_mu\t'%(each))
    for each, h in node.partial_sigmas.iteritems():
        newf.write('%d_sig\t'%(each))
    newf.write('node.brow_mean\t')
    newf.write('node.current_brow_mean_ln_prior_density\t')
    newf.write('ln_posterior\t')
#     newf.write('ln_likelihood\t')
#     newf.write('node.current_brow_ln_like\t')
    newf.write('\n')

    output = os.path.join('model.txt')

#     newf2 = open(output, 'w')
#     newf2.write('%s\t%s\n\n'%('data-generating means from brownian motion=', str(rand_simulated_means)))
#     newf2.write('%s\t'%('n_gen'))
#     for i in (node.confset):
#         newf2.write('%s\t'%(i))
#     newf2.write('\n')

    output3 = os.path.join('model.txt')
    newf2 = open(output3, 'w')
    newf2.write('%s\t'%('n_gen'))
    for i in (node.confset):
        newf2.write('%s\t'%(i))
    newf2.write('\n')


#     output2 = os.path.join('posterior_samples.txt')
#     newf3 = open(output2, 'w')
#     newf3.write(' **********starting values ******************************************************\n')
#     newf3.write('Starting_topology=%s\t\n\n'%(node.current_topology))
#     newf3.write('%s\t\n\n'%(node.vcv))
#     newf3.write('#############################\n')
#     newf3.write('#####log-prior-densities#####\n')
#     newf3.write('node.current_mean_ln_prior_density = %s\t\n'%(node.current_mean_ln_prior_density))
#     newf3.write('node.current_sigma_ln_prior_density = %s\t\n'%(node.current_sigma_ln_prior_density))
#     newf3.write('node.current_brow_mean_ln_prior_density = %s\t\n'%(node.current_brow_mean_ln_prior_density))
#     newf3.write('\n')

    #print
    #print
    list_of_sets = []

    for i in range(n_gen):
        pick = random.choice(postorder)
        print 'iter and picked node#########################################################################################################################################################>', i+1, pick.number
        print 'pick.number=>', pick.number
        
        if (pick.lchild == None) and (contained_tips.get(pick.number) in node.confset):
            movetype = random.choice(('movedown', 'mean_sigma'))
#             print 'movetype====', movetype
            if movetype == 'mean_sigma':
                movetype = random.choice(('mean', 'sigma', 'brown_mean'))
                moves(pick, movetype)
            else:
                moves(pick, movetype)
            set_string = " ".join(str(x) for x in node.confset)
            list_of_sets.append(set_string)
        
        elif (pick.number is not -1) and (pick.lchild is not None) and (contained_tips.get(pick.number) in node.confset):
            if (contained_tips.get(pick.par.lchild.number) in node.confset) and (contained_tips.get(pick.par.lchild.rsib.number) in node.confset):
                movetype = random.choice(('mov',  'mean_sigma' ))
                
                if movetype == 'mean_sigma':
                    movetype = random.choice(('mean', 'sigma', 'brown_mean' ))
                    moves(pick, movetype)
                else:
                    movetype = random.choice(('moveup', 'movedown' ))
                    moves(pick, movetype)
 
            else:
                movetype = random.choice(('moveup', 'mean_sigma'))
                if movetype == 'mean_sigma':
                    movetype = random.choice(('mean', 'sigma', 'brown_mean' ))
                    moves(pick, movetype)
                else:
                    moves(pick, movetype)
            set_string = " ".join(str(x) for x in node.confset)
            list_of_sets.append(set_string)

        elif (pick.number is -1):
            #print 'This is the root=====>', pick.number

            current_node = contained_tips.get(pick.number)
            if current_node in node.confset:
                movetype = random.choice(('moveup', 'mean_sigma'))
                if movetype == 'mean_sigma':
                    movetype = random.choice(('mean', 'sigma', 'brown_mean'))
                    moves(pick, movetype)
                else:
                    moves(pick, movetype)
                set_string = " ".join(str(x) for x in node.confset)
                list_of_sets.append(set_string)

            else:
                node.types = 'Pass'
                print 'passing node -1 does not qualify=', pick.number

        else:
            node.types = 'Pass'
            print 'passing does not qualify=', pick.number

        if (i+1) % save_every == 0:
            newf.write('%s\t'%(mcmc+1))

            for y, h in node.partial_means.iteritems():
                newf.write('%.6f\t'%(h))
            for m, n in node.partial_sigmas.iteritems():
                newf.write('%.6f\t'%(n))
            newf.write('%.6f\t'%(node.brow_mean))
            newf.write('%.6f\t'%(node.current_brow_mean_ln_prior_density))
            newf.write('%.6f\t'%(node.current_ln_posterior))
#             newf.write('%.6f\t'%(node.current_ln_like))
#             newf.write('%.6f\t'%(node.current_brow_ln_like))

            newf.write('\n')


        if node.types == ('Accept') or node.types ==  ('Reject'):
#             print 'yeeeeeeeeeeeeeeeeeee', node.types

            if (i+1) % save_every == 0:
                newf2.write('%s\t'%(mcmc+1))
            
                collt = []
                for i,list in enumerate(node.confset):
                    collt.append([i+1]*(len(list)))
            
            
 
                for i in collt:
                    newf2.write(('%s\t\t') %("       ".join(str(x) for x in i)))
                
                newf2.write('%.6f\t%.6f\t%s\t%.6f\t%.6f\t%s\t'%(node.current_ln_posterior, node.current_ln_like, node.types, node.delta1, node.lam1,movetype))
				
                newf2.write('\n')
        else:
            pass
        print
#             set_string = " ".join(str(x) for x in node.confset)
#             list_of_sets.append(set_string)

#             newf3.write('\n')
#             newf3.write('iter and picked node ###################################################################################>%s\t%s\t\n'%(mcmc+1,pick.number))
#             newf3.write('%s\n'%(movetype))
#             newf3.write('ACCEPT or REJECT or PASS ==> %s\t\n\n'%(node.types))
#             newf3.write('Accepted_topology=%s\t\n\n'%(node.current_topology))
#             newf3.write('%s\t\n\n'%(node.vcv))
#             newf3.write('#############################\n')
#             newf3.write('#####log-prior-densities#####\n')
#             newf3.write('node.current_mean_ln_prior_density = %s\t\n'%(node.current_mean_ln_prior_density))
#             newf3.write('node.current_sigma_ln_prior_density = %s\t\n'%(node.current_sigma_ln_prior_density))
#             newf3.write('node.current_brow_mean_ln_prior_density = %s\t\n'%(node.current_brow_mean_ln_prior_density))




        mcmc+=1
    newf2.write('\n\n')
#     newf3.write('\n\n')
    myDict=Counter(list_of_sets)
    print myDict
    numadd = 0
    for ele in myDict:
        numadd+=myDict[ele]
    for ele in myDict:
        print (ele, myDict[ele])
        newf2.write('%s\t%.6f\t\t\n'%(ele, float(myDict[ele])/float(numadd)))
#         print float(myDict[ele])/float(numadd)

#    sys.stdout = old_stdout
#    log_file.close()


if __name__ == '__main__':
    exploreNodes(mega_list, n_gen)
