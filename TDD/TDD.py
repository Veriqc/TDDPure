import numpy as np
import copy
import time
import random
from TDD.ComplexTable import *

"""Define global variables"""
computed_table = dict()
unique_table = dict()
global_index_order = dict()
global_node_idx=0
add_find_time=0
add_hit_time=0
cont_find_time=0
cont_hit_time=0
add_hit_time2=0
add_times = 0
cont_times = 0
# epi=1e-5

terminal_node = None

class Index:
    """The index in the form (name,hyper_idx), here name is expected to be a string, hyper_idx is used to cope with hyper index"""
    def __init__(self, name ,hyper_idx=0):
        self.name = name
        self.hyper_idx = hyper_idx
        
    def __eq__(self,other):
        if self.name == other.name and self.hyper_idx == other.hyper_idx:
            return True
        else:
            return False
        
    def __lt__(self,other):
        if global_index_order[self.name] < global_index_order[other.name]:
            return True
        elif self.name == other.name and self.hyper_idx<other.hyper_idx:
            return True
        
        return False
    
    def __str__(self):
        return str((self.name,self.hyper_idx))

    
class Node:
    """To define the node of TDD"""
    def __init__(self,key,num=2):
        self.key = key
        self.succ_num = num
        self.id = None
        self.out_weight = None #suppose to be a list of length num
        self.succ = None #suppose to be a list of length num
        self.ref_num = None


class TDD:
    def __init__(self,node,weight=cn1):
        """TDD"""
        self.weight = weight
        self.node = node
        self.index_set = None #suppose to be a list or set
        self.key_2_index = None #suppose to be a dict
            
            
    def node_number(self):
        node_set=set()
        node_set=get_node_set(self.node,node_set)
        return len(node_set)
    
    def self_copy(self):
        temp = TDD(self.node)
        temp.weight = self.weight
        temp.index_set = copy.copy(self.index_set)
        temp.key_2_index=copy.copy(self.key_2_index)
        return temp
    
    def show(self,real_label=True):
        from graphviz import Digraph
        from IPython.display import Image
        edge=[]              
        dot=Digraph(name='reduced_tree')
        dot=layout(self.node,self.key_2_index,dot,edge,real_label)
        dot.node('-0','',shape='none')
        dot.edge('-0',str(self.node.id),color="blue",label=str(complex(round(self.weight.r.val,2),round(self.weight.i.val,2))))
        dot.format = 'png'
        return Image(dot.render('output'))
    
    def to_array(self,var=[]):
        split_pos=0
        key_repeat_num=dict()
        var_idx=dict()       
        if var:
            for idx in var:
                if not idx.key in var_idx:
                    var_idx[idx.key]=1
                else:
                    var_idx[idx.key]+=1
        elif self.index_set:
            for idx in self.index_set:
                if not idx.key in var_idx:
                    var_idx[idx.key]=1
                else:
                    var_idx[idx.key]+=1
        if var:
            split_pos=len(var_idx)-1
        elif self.key_2_index:
            split_pos=max(self.key_2_index)
        else:
            split_pos=self.node.key
        orig_order=[]
        for k in range(split_pos+1):
            if k in self.key_2_index:
                if self.key_2_index[k] in var_idx:
                    key_repeat_num[k] = var_idx[self.key_2_index[k]]
            else:
                key_repeat_num[k]=1
            if k in self.key_2_index:
                for k1 in range(key_repeat_num[k]):
                    orig_order.append(self.key_2_index[k])
                     

        res = tdd_2_np(self,split_pos,key_repeat_num)

        return res

    def measure(self,split_pos=None):
        res=[]
        get_measure_prob(self)
        if split_pos==None:
            if self.key_2_index:
                split_pos=max(self.key_2_index)
            else:
                split_pos=self.node.key

        if split_pos==-1:
            return ''
        else:
            if split_pos!=self.node.key:
                l=random.randint(0,1)
                temp_res=self.measure(split_pos-1)
                res=str(l)+temp_res
                return res
            l=random.uniform(0,sum(self.node.meas_prob))
            if l<self.node.meas_prob[0]:
                temp_tdd=Slicing(self,self.node.key,0)
                temp_res=temp_tdd.measure(split_pos-1)
                res='0'+temp_res
            else:
                temp_tdd=Slicing(self,self.node.key,1)
                temp_res=temp_tdd.measure(split_pos-1)
                res='1'+temp_res
#         print(res)
        return res
    def get_amplitude(self,b):
        """b is the term for calculating the amplitude"""
        if len(b)==0:
            return self.weight.r.val+1j*self.weight.i.val
        
        if len(b)!=self.node.key+1:
            b.pop(0)
            return self.get_amplitude(b)
        else:
            temp_tdd=Slicing(self,self.node.key,b[0])
            b.pop(0)
            res=temp_tdd.get_amplitude(b)
            w=res*self.weight
            return w.r.val+1j*w.i.val
            
    def sampling(self,k):
        res=[]
        for k1 in range(k):
            temp_res=self.measure()
            res.append(temp_res )
#         print(res)
        return res
        
        
    def __eq__(self,other):
        if self.node==other.node and self.weight == other.weight:# and self.key_2_index==other.key_2_index
            return True
        else:
            return False
        
def layout(node,key_2_idx,dot=None,succ=[],real_label=True):
    col=['red','blue','black','green']
    if real_label and node.key in key_2_idx:
        if node.key==-1:
            dot.node(str(node.id), str(1), fontname="helvetica",shape="circle",color="red")
        else:
            dot.node(str(node.id), key_2_idx[node.key], fontname="helvetica",shape="circle",color="red")
    else:
        dot.node(str(node.id), str(node.key), fontname="helvetica",shape="circle",color="red")
    for k in range(node.succ_num):
        if node.succ[k]:
            label1=str(complex(round(node.out_weight[k].r.val,2),round(node.out_weight[k].i.val,2)))
            if not node.succ[k] in succ:
                dot=layout(node.succ[k],key_2_idx,dot,succ,real_label)
                dot.edge(str(node.id),str(node.succ[k].id),color=col[k%4],label=label1)
                succ.append(node.succ[k])
            else:
                dot.edge(str(node.id),str(node.succ[k].id),color=col[k%4],label=label1)
    return dot        

        
def Ini_TDD(index_order=[],max_rank=200):
    """To initialize the unique_table,computed_table and set up a global index order"""
    global unique_table,computed_table,terminal_node,global_node_idx
    global add_find_time,add_hit_time,cont_find_time,cont_hit_time,add_times,cont_times

    global_node_idx=0
    add_find_time=0
    add_hit_time=0
    cont_find_time=0
    cont_hit_time=0
    add_times = 0
    cont_times = 0

    unique_table = {k:dict() for k in range(max_rank)}
    computed_table = dict()
    computed_table['+'] = {k:{k1: dict() for k1 in range(-1,max_rank)} for k in range(-1,max_rank)}
    computed_table['*'] = {k:{k1: dict() for k1 in range(-1,max_rank)} for k in range(-1,max_rank)}
    set_index_order(index_order)
    ini_complex(max_rank)
    terminal_node = Find_Or_Add_Unique_table(-1)
    return get_identity_tdd()

def Clear_TDD():
    """To initialize the unique_table,computed_table and set up a global index order"""
    global unique_table,computed_table,terminal_node,global_node_idx
    global add_find_time,add_hit_time,cont_find_time,cont_hit_time
    global_node_idx=0
    unique_table.clear()
    computed_table['+'].clear()
    computed_table['*'].clear()
    add_find_time=0
    add_hit_time=0
    cont_find_time=0
    cont_hit_time=0
    global_node_idx=0


def get_identity_tdd():
    tdd = TDD(terminal_node,cn1)
    tdd.index_set = []
    tdd.key_2_index = {}
    return tdd

def get_unique_table():
    return unique_table

def get_unique_table_num():
    return len(unique_table)

def set_index_order(var_order):
    """the index with a bigger value should appear on the top"""
    global global_index_order
    global_index_order=dict()
    if isinstance(var_order,list):
        for k in range(len(var_order)):
            global_index_order[var_order[k]] = k
    if isinstance(var_order,dict):
        global_index_order = copy.copy(var_order)
    global_index_order[-1] = -1
    
def get_index_order():
    global global_index_order
    return global_index_order
    
def get_int_key(weight):
    """To transform a complex number to a tuple with int values"""
    return (int(weight.r.val*epi_inv) ,int(weight.i.val*epi_inv))

def get_node_set(node,node_set=set()):
    """Only been used when counting the node number of a TDD"""
    if not node in node_set:
        node_set.add(node)
        for k in range(node.succ_num):
            if node.succ[k]:
                node_set = get_node_set(node.succ[k],node_set)
    return node_set

def Find_Or_Add_Unique_table(x,weigs=[],succ_nodes=[]):
    """To return a node if it already exist, creates a new node otherwise"""
    global global_node_idx,unique_table
    
    if x==-1:
        if unique_table.__contains__(x):
            return unique_table[x]
        else:
            res=Node(x,0)
            res.id=0
            unique_table[x]=res
        return res
    temp_key = [(w.r, w.i) for w in weigs] + succ_nodes

    temp_key=tuple(temp_key)
    if temp_key in unique_table[x]:
        return unique_table[x][temp_key]
    else:
        res=Node(x,len(succ_nodes))
        global_node_idx+=1
        res.id = global_node_idx
        res.out_weight = weigs
        res.succ = succ_nodes
        unique_table[x][temp_key]=res
    return res
    

def normalize(x,the_successors,cached = False):
    """The normalize and reduce procedure"""

    m = len(the_successors)
    is_zero=[equalsZero(the_successors[k].weight) for k in range(m)]
    if cached:
        for k in range(0,m):
            if is_zero[k] and the_successors[k].weight!=cn0:
                releaseCached(the_successors[k].weight)
                the_successors[k].weight = cn0
            
    
#     weigs_abs=[int(round(succ.weight.norm()/epi)) for succ in the_successors]
#     max_pos = weigs_abs.index(max(weigs_abs))
#     weig_max = the_successors[max_pos].weight
    
    max_pos = -1
    weig_max = cn1
    for k in range(m):
        if is_zero[k]: continue
        if max_pos == -1:
            max_pos = k
            weig_max = the_successors[k].weight
            weig_max_norm = weig_max.norm()
        else:
            mag = the_successors[k].weight.norm()
            if mag - weig_max_norm>epi:
                max_pos=k
                weig_max = the_successors[k].weight
                weig_max_norm = mag
    
    if max_pos == -1:
        return TDD(terminal_node,cn0)
            
    weigs = []
    for k in range(m):
        if k==max_pos:
            weigs.append(cn1)
        else:
            if the_successors[k].weight==cn0:
                weigs.append(cn0)
            elif (not cached) and equalsOne(weig_max):
                weigs.append(the_successors[k].weight)
            else:
                weigs.append(Find_Or_Add_Complex_table(the_successors[k].weight/weig_max))
    
    succ_nodes=[succ.node for succ in the_successors]
    
    all_equal=True
    for k in range(1,m):
        if succ_nodes[k]!=succ_nodes[0]:
            all_equal=False
            break
        if weigs[k]!=weigs[0]:
            all_equal=False
            break            
    if all_equal:
        if not cached:
            return the_successors[0]
        if cached:
            for k in range(1,len(the_successors)):
                if the_successors[k].weight!=cn0 and the_successors[k].weight!=cn1:
                    releaseCached(the_successors[k].weight)
            return the_successors[0]
    
    
    node=Find_Or_Add_Unique_table(x,weigs,succ_nodes)
    
    if cached:
        res=TDD(node,weig_max)
    else:
        res=TDD(node,Find_Or_Add_Complex_table(weig_max))
        
    if cached:
        for k in range(m):
            if k ==max_pos:
                continue
            if the_successors[k].weight!=cn0 and the_successors[k].weight!=cn1:
                releaseCached(the_successors[k].weight)

    return res
              
              

def get_count():
    global add_find_time,add_hit_time,cont_find_time,cont_hit_time,add_times,cont_times
    print("add:",add_hit_time,'/',add_find_time,'/',add_hit_time/add_find_time)
    print("cont:",cont_hit_time,"/",cont_find_time,"/",cont_hit_time/cont_find_time)
    print('add:',add_times,'cont:',cont_times)

ttt_table = dict()

ttt_table2 = dict()

ttt_n = 0

ttt_n2 = 0    
    
def find_computed_table(item):
    """To return the results that already exist"""
    global computed_table,add_find_time,add_hit_time,cont_find_time,cont_hit_time,add_hit_time2
    if item[0] == '+':
            add_find_time+=1
#         if item[1].node.id > item[2].node.id:
#             the_key=(get_int_key(item[2].weight),item[2].node,get_int_key(item[1].weight),item[1].node)
#             k1 = item[2].node.key
#             k2 = item[1].node.key
#             if the_key in computed_table['+'][item[2].node.key][item[1].node.key]:
#                 res = computed_table['+'][item[2].node.key][item[1].node.key][the_key]
#                 add_hit_time+=1
#                 if abs(res[0])<epi and abs(res[1])<epi:
#                     return TDD(terminal_node,cn0)
#                 else:
#                     tdd = TDD(res[2],getCachedComplex2(res[0],res[1]))
#                     return tdd            
#         else:
            the_key=(get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node)
            if the_key in computed_table['+'][item[1].node.key][item[2].node.key]:
                res = computed_table['+'][item[1].node.key][item[2].node.key][the_key]
                add_hit_time+=1
                if abs(res[0])<epi and abs(res[1])<epi:
                    return TDD(terminal_node,cn0)
                else:
                    tdd = TDD(res[2],getCachedComplex2(res[0],res[1]))
                    return tdd
    else:
        cont_find_time+=1
        
        if item[1].node.id > item[2].node.id:
            the_key=(item[2].node,item[1].node,item[4],item[3],item[5])
            
            if not the_key in ttt_table:
                ttt_table[the_key]=0
            
            if computed_table['*'][item[2].node.key][item[1].node.key].__contains__(the_key):
                res = computed_table['*'][item[2].node.key][item[1].node.key][the_key]
                ttt_table[the_key]+=1
                cont_hit_time+=1
                if abs(res[0])<epi and abs(res[1])<epi:
                    return TDD(terminal_node,cn0)
                else:
                    tdd = TDD(res[2],getCachedComplex2(res[0],res[1]))
                    return tdd
        else: 
            the_key=(item[1].node,item[2].node,item[3],item[4],item[5])
            
            if not the_key in ttt_table:
                ttt_table[the_key]=0            
            
            if computed_table['*'][item[1].node.key][item[2].node.key].__contains__(the_key):
                res = computed_table['*'][item[1].node.key][item[2].node.key][the_key]
                ttt_table[the_key]+=1
                cont_hit_time+=1
                if abs(res[0])<epi and abs(res[1])<epi:
                    return TDD(terminal_node,cn0)
                else:
                    tdd = TDD(res[2],getCachedComplex2(res[0],res[1]))
                    return tdd
    return None

def insert_2_computed_table(item,res):
    """To insert an item to the computed table"""
    global computed_table,cont_time,find_time,hit_time

    if item[0] == '+':
        
#         global ttt_table2,ttt_n2
#         h = hash((get_int_key(item[2].weight),item[2].node,get_int_key(item[1].weight),item[1].node))
#         if not h in ttt_table2:
#             ttt_table2[h]=(get_int_key(item[2].weight),item[2].node,get_int_key(item[1].weight),item[1].node)
#         else:
#             if ttt_table2[h] != (get_int_key(item[2].weight),item[2].node,get_int_key(item[1].weight),item[1].node):
#                 ttt_n2+=1
#                 ttt_table2[h] = (get_int_key(item[2].weight),item[2].node,get_int_key(item[1].weight),item[1].node)
                
#         h = hash((get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node))
#         if not h in ttt_table2:
#             ttt_table2[h]=(get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node)
#         else:
#             if ttt_table2[h] != (get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node):
#                 ttt_n2+=1
#                 ttt_table2[h] = (get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight),item[2].node)                
        
#         if item[1].node.id>item[2].node.id:
#             computed_table['+'][item[2].node.key][item[1].node.key][(get_int_key(item[2].weight) ,item[2].node,get_int_key(item[1].weight),item[1].node)] = (res.weight.r.val,res.weight.i.val,res.node)
#         else:
            computed_table['+'][item[1].node.key][item[2].node.key][(get_int_key(item[1].weight),item[1].node,get_int_key(item[2].weight) ,item[2].node)] = (res.weight.r.val,res.weight.i.val,res.node)
    else:
#         global ttt_table,ttt_n
#         h = hash((item[1].node,item[2].node,item[3],item[4],item[5]))
#         if not h in ttt_table:
#             ttt_table[h]=(item[1].node,item[2].node,item[3],item[4],item[5])
#         else:
#             if ttt_table[h] != (item[1].node,item[2].node,item[3],item[4],item[5]):
#                 ttt_n+=1
#                 ttt_table[h] = (item[1].node,item[2].node,item[3],item[4],item[5])
                
#         h = hash((item[2].node,item[1].node,item[4],item[3],item[5]))
#         if not h in ttt_table:
#             ttt_table[h]=(item[2].node,item[1].node,item[4],item[3],item[5])
#         else:
#             if ttt_table[h] != (item[2].node,item[1].node,item[4],item[3],item[5]):
#                 ttt_n+=1
#                 ttt_table[h] = (item[2].node,item[1].node,item[4],item[3],item[5])               
        
        if item[1].node.id>item[2].node.id:
            computed_table['*'][item[2].node.key][item[1].node.key][(item[2].node,item[1].node,item[4],item[3],item[5])] = (res.weight.r.val,res.weight.i.val,res.node)
        else:
            computed_table['*'][item[1].node.key][item[2].node.key][(item[1].node,item[2].node,item[3],item[4],item[5])] = (res.weight.r.val,res.weight.i.val,res.node)
            
        
        
        
    
def get_index_2_key(var):
    var_sort = copy.copy(var)
    var_sort.sort()
    n=0
    idx_2_key={}
    key_2_idx={}
    for idx in var_sort:
        if not idx.name in idx_2_key:
            idx_2_key[idx.name] = n
            key_2_idx[n] = idx.name
            n+=1
    return key_2_idx,idx_2_key
    
def get_tdd(U,var=[]):
 
    key_2_idx,idx_2_key = get_index_2_key(var)
    
    order=[]
    
    for idx in var:
        order.append(idx_2_key[idx.name])
        
    tdd = np_2_tdd(U,order)
    tdd.index_set = var
    tdd.key_2_index = key_2_idx
    return tdd

def np_2_tdd(U,order=[]):
    #index is the index_set as the axis order of the matrix
    U_dim=U.ndim
    U_shape=U.shape
    if sum(U_shape)==U_dim:
        for k in range(U_dim):
            U=U[0]
        res=TDD(terminal_node,Find_Or_Add_Complex_table(getTempCachedComplex2(U)))
        return res
    
    if not order:
        order=list(range(U_dim))
            
    x=max(order)
    split_pos=order.index(x)
    order[split_pos]=-1
    split_U=np.split(U,U_shape[split_pos],split_pos)
    
    while x in order:
        split_pos=order.index(x)
        for k in range(len(split_U)):
            split_U[k]=np.split(split_U[k],U_shape[split_pos],split_pos)[k]
        order[split_pos]=-1
    
    the_successors=[]
    for k in range(U_shape[split_pos]):
        res=np_2_tdd(split_U[k],copy.copy(order))
        the_successors.append(res)
    tdd = normalize(x,the_successors)
    
    return tdd
    
    
def np_2_tdd2(U,split_pos=None):
    #index is the index_set as the axis order of the matrix
    U_dim=U.ndim
    U_shape=U.shape
    if sum(U_shape)==U_dim:
        for k in range(U_dim):
            U=U[0]
        res=TDD(terminal_node,Find_Or_Add_Complex_table(getTempCachedComplex2(U)))
        return res
    if split_pos==None:
        split_pos=U_dim-1
        
    split_U=np.split(U,U_shape[split_pos],split_pos)
    the_successors=[]
    for k in range(U_shape[split_pos]):
        res=np_2_tdd(split_U[k],split_pos-1)
        the_successors.append(res)
    tdd = normalize(split_pos,the_successors)
    for k in range(len(U_shape)):
        tdd.key_width[k]=U_shape[k]
    return tdd
    
def tdd_2_np(tdd,split_pos=None,key_repeat_num=dict()):
#     print(split_pos,key_repeat_num)
    if split_pos==None:
        split_pos=tdd.node.key
            
    if split_pos==-1:
        return tdd.weight.r.val+1j*tdd.weight.i.val
    else:
        the_succs=[]
        for k in range(tdd.key_width[split_pos]):
            succ=Slicing2(tdd,split_pos,k)
            succ.key_width=tdd.key_width
            temp_res=tdd_2_np(succ,split_pos-1,key_repeat_num)
            the_succs.append(temp_res)
        if not split_pos in key_repeat_num:
            r = 1
        else:
            r = key_repeat_num[split_pos]
            
        if r==1:
            res=np.stack(tuple(the_succs), axis=the_succs[0].ndim)
        else:
            new_shape=list(the_succs[0].shape)
            for k in range(r):
                new_shape.append(tdd.key_width[split_pos])
            res=np.zeros(new_shape)
            for k1 in range(tdd.key_width[split_pos]):
                f='res['
#                 print(the_succs[0].ndim,r-1)
                for k2 in range(the_succs[0].ndim):
                    f+=':,'
                for k3 in range(r-1):
                    f+=str(k1)+','
                f=f[:-1]+']'
                eval(f)[k1]=the_succs[k1]
        return res
    
    
"""need modify"""                
def get_measure_prob(tdd):
    if tdd.node.meas_prob:
        return tdd
    if tdd.node.key==-1:
        tdd.node.meas_prob=[0.5,0.5]
        return tdd
    if not tdd.node.succ_num==2:
        print("Only can be used for binary quantum state")
        return tdd
    get_measure_prob(Slicing(tdd,tdd.node.key,0))
    get_measure_prob(Slicing(tdd,tdd.node.key,1))
    tdd.node.meas_prob=[0]*2
    p0=tdd.node.out_weight[0].r.val*tdd.node.out_weight[0].r.val+tdd.node.out_weight[0].i.val*tdd.node.out_weight[0].i.val
    p1=tdd.node.out_weight[1].r.val*tdd.node.out_weight[1].r.val+tdd.node.out_weight[1].i.val*tdd.node.out_weight[1].i.val
                
    tdd.node.meas_prob[0]=p0*sum(tdd.node.succ[0].meas_prob)*2**(tdd.node.key-tdd.node.succ[0].key-1)
    tdd.node.meas_prob[1]=p1*sum(tdd.node.succ[1].meas_prob)*2**(tdd.node.key-tdd.node.succ[1].key-1)
    return tdd



class new_key_node:
    def __init__(self,level,new_key):
        self.level = level
        self.new_key = new_key
        self.next = dict()
        self.father = None
        
    def __eq__(self,other):
        if self.level==other.level and self.new_key==other.new_key:
            return true
        else:
            return False
        
    def __hash__(self):
        return id(self)
    
    def __str__(self):
        temp = self
        s=[]
        while temp.level!=-1:
            s.append(temp.new_key)
            temp = temp.father
        return str(s)
        
    def append_new_key(self,new_key):
        if new_key in self.next:
            return self.next[new_key]
        else:
            temp = new_key_node(self.level+1,new_key)
            self.next[new_key] = temp
            temp.father = self
            return temp


key_2_new_key_tree_header = new_key_node(-1,-1)

    
def cont(tdd1,tdd2):
    
    var_cont = [] #to record the indices to be contracted
    var_out = [] #ro record the indices to be remained
    var_out_name = []
    var_cont_name = []
    
    for var1 in tdd1.index_set:
        if not var1 in tdd2.index_set:
            var_out.append(var1)
            var_out_name.append(var1.name)
        else:
            var_cont.append(var1)
            
    for var2 in tdd2.index_set:
        if not var2 in tdd1.index_set:
            var_out.append(var2)
            var_out_name.append(var2.name)
            
    for var in var_cont:
        if var.name in var_out_name:
            var_cont.remove(var)
            
    var_cont_name = [idx.name for idx in var_cont]
            
    k1 = 0
    k2 = 0
    new_key = 0
    m1 = len(tdd1.key_2_index)
    m2 = len(tdd2.key_2_index)
    key_2_new_key1 = key_2_new_key_tree_header
    key_2_new_key2 = key_2_new_key_tree_header
    new_key_2_index = dict()

    while k1 < m1 or k2 < m2 :
        if k1 == m1: 
            for k2 in range(k2,m2):
                key_2_new_key2 = key_2_new_key2.append_new_key(new_key)
                new_key_2_index[new_key]=tdd2.key_2_index[k2]
                new_key+=1
            break
                
        if k2 == m2: 
            for k1 in range(k1,m1):
                key_2_new_key1 = key_2_new_key1.append_new_key(new_key)
                new_key_2_index[new_key]=tdd1.key_2_index[k1]
                new_key+=1
            break
        if global_index_order[tdd1.key_2_index[k1]] < global_index_order[tdd2.key_2_index[k2]]:
            key_2_new_key1 = key_2_new_key1.append_new_key(new_key)
            new_key_2_index[new_key]=tdd1.key_2_index[k1]
            new_key+=1
            k1+=1
        elif global_index_order[tdd1.key_2_index[k1]] > global_index_order[tdd2.key_2_index[k2]]:
            key_2_new_key2 = key_2_new_key2.append_new_key(new_key)
            new_key_2_index[new_key]=tdd2.key_2_index[k2]
            new_key+=1
            k2+=1
        elif tdd1.key_2_index[k1] in var_out_name:
            key_2_new_key1 = key_2_new_key1.append_new_key(new_key)
            key_2_new_key2 = key_2_new_key2.append_new_key(new_key)
            new_key_2_index[new_key]=tdd1.key_2_index[k1]
            new_key+=1
            k1+=1
            k2+=1
        else:
            key_2_new_key1 = key_2_new_key1.append_new_key(new_key-0.5)
            key_2_new_key2 = key_2_new_key2.append_new_key(new_key-0.5)
            k1+=1
            k2+=1
    
#     print(key_2_new_key1)
#     print(key_2_new_key2)
    
    cacheCount_in=cacheCount
    tdd=contract(tdd1,tdd2,key_2_new_key1,key_2_new_key2,len(set(var_cont_name)))
    if not tdd.weight==cn0 and not tdd.weight==cn1:
        releaseCached(tdd.weight)
    tdd.weight=Find_Or_Add_Complex_table(tdd.weight)

    if not cacheCount==cacheCount_in:
        print('Something went wrong, cacheCount not match')
    
    tdd.index_set=var_out
    tdd.key_2_index = new_key_2_index

    return tdd
    

def contract(tdd1,tdd2,key_2_new_key1,key_2_new_key2,cont_num):
    """The contraction of two TDDs, var_cont is in the form [[4,1],[3,2]]"""
    global cont_times
    cont_times+=1
    
    k1=tdd1.node.key
    k2=tdd2.node.key
    w1=tdd1.weight
    w2=tdd2.weight

    if w1== cn0 or w2==cn0:
        return TDD(terminal_node,cn0)    
    
    if k1==-1 and k2==-1:
        tdd=TDD(terminal_node,cn_mulCached(w1,w2))
        if cont_num>0:
            temp = getTempCachedComplex2(2**cont_num)
            cn_mul(tdd.weight,tdd.weight,temp)
        return tdd
    

        
    while key_2_new_key2.level > k2:
        key_2_new_key2=key_2_new_key2.father        

    if k1==-1:
        if cont_num ==0 and key_2_new_key2.new_key==k2:
            tdd=TDD(tdd2.node,cn_mulCached(w1,w2))
            return tdd
            
    while key_2_new_key1.level > k1:
        key_2_new_key1=key_2_new_key1.father
        
    if k2==-1:      
        if cont_num ==0 and key_2_new_key1.new_key==k1:
            tdd=TDD(tdd1.node,cn_mulCached(w1,w2))
            return tdd
    
    tdd1.weight = cn1
    tdd2.weight = cn1

    tdd = find_computed_table(['*',tdd1,tdd2,key_2_new_key1,key_2_new_key2,cont_num])
    if tdd:
        tdd1.weight=w1
        tdd2.weight=w2
        if not tdd.weight==cn0:
            cn_mul(tdd.weight,tdd.weight,w1)
            cn_mul(tdd.weight,tdd.weight,w2)
            if equalsZero(tdd.weight):
                releaseCached(tdd.weight)
                return TDD(terminal_node,cn0)
        return tdd
                
    new_key1 = key_2_new_key1.new_key
    new_key2 = key_2_new_key2.new_key
        
    if new_key1 > new_key2: #the bigger new_key is supposed to be coped first
        temp_key_2_new_key2 = key_2_new_key2
        if (2*new_key1)%2==0:
            the_successors=[]
            for k in range(tdd1.node.succ_num):
                e1 = TDD(tdd1.node.succ[k],tdd1.node.out_weight[k])
                res=contract(e1,tdd2,key_2_new_key1.father,temp_key_2_new_key2,cont_num)
                the_successors.append(res)
            tdd=normalize(new_key1,the_successors,True)
        else:
            tdd=TDD(terminal_node,cn0)
            for k in range(tdd11.node.succ_num):
                e1 = TDD(tdd1.node.succ[k],tdd1.node.out_weight[k])
                res=contract(e1,tdd2,key_2_new_key1.father,temp_key_2_new_key2,cont_num-1)
                if tdd.weight==cn0:
                    tdd=res
                elif res.weight != cn0:
                    old_w = tdd.weight
                    tdd=add(tdd,res)
                    releaseCached(old_w)
                    releaseCached(res.weight)
    elif new_key1 == new_key2:
        if (2*new_key1)%2==0:
            the_successors=[]
            for k in range(tdd1.node.succ_num):
                e1 = TDD(tdd1.node.succ[k],tdd1.node.out_weight[k])
                e2 = TDD(tdd2.node.succ[k],tdd2.node.out_weight[k])
                res=contract(e1,e2,key_2_new_key1.father,key_2_new_key2.father,cont_num)
                the_successors.append(res)
            tdd=normalize(new_key1,the_successors,True)
        else:
            tdd=TDD(terminal_node,cn0)
            for k in range(tdd1.node.succ_num):
                e1 = TDD(tdd1.node.succ[k],tdd1.node.out_weight[k])
                e2 = TDD(tdd2.node.succ[k],tdd2.node.out_weight[k])                
                res=contract(e1,e2,key_2_new_key1.father,key_2_new_key2.father,cont_num-1)           
                if tdd.weight==cn0:
                    tdd=res
                elif res.weight != cn0:
                    old_w = tdd.weight
                    tdd=add(tdd,res)
                    releaseCached(old_w)
                    releaseCached(res.weight)
    else:
        temp_key_2_new_key1 = key_2_new_key1
        if (2*new_key2)%2==0:
            the_successors=[]
            for k in range(tdd2.node.succ_num):
                e2 = TDD(tdd2.node.succ[k],tdd2.node.out_weight[k])
                res=contract(tdd1,e2,temp_key_2_new_key1,key_2_new_key2.father,cont_num)
                the_successors.append(res)
            tdd=normalize(new_key2,the_successors,True)
        else:
            tdd=TDD(terminal_node,cn0)
            for k in range(tdd2.node.succ_num):
                e2 = TDD(tdd2.node.succ[k],tdd2.node.out_weight[k])
                res=contract(tdd1,e2,temp_key_2_new_key1,key_2_new_key2.father,cont_num-1)           
                if tdd.weight==cn0:
                    tdd=res
                elif res.weight != cn0:
                    old_w = tdd.weight
                    tdd=add(tdd,res)
                    releaseCached(old_w)
                    releaseCached(res.weight)
     
    insert_2_computed_table(['*',tdd1,tdd2,key_2_new_key1,key_2_new_key2,cont_num],tdd)
    tdd1.weight=w1
    tdd2.weight=w2
    if not tdd.weight==cn0 and (w1!=cn1 or w2!=cn1):
        if tdd.weight==cn1:
            tdd.weight = cn_mulCached(w1,w2)
        else:
            cn_mul(tdd.weight,tdd.weight,w1)
            cn_mul(tdd.weight,tdd.weight,w2)
        if equalsZero(tdd.weight):
            releaseCached(tdd.weight)
            return TDD(terminal_node,cn0)
    return tdd
    
def Slicing(tdd,x,c):
    """Slice a TDD with respect to x=c"""

    k=tdd.node.key
    
    if k==-1:
        res = TDD(tdd.node,tdd.weight)
        return res
    
    if k<x:
        res = TDD(tdd.node,tdd.weight)
        return res
    
    if k==x:
        res=TDD(tdd.node.succ[c],tdd.node.out_weight[c])
        return res
    else:
        print("Not supported yet!!!")
        
        
def Slicing2(tdd,x,c):
    """Slice a TDD with respect to x=c"""

    k=tdd.node.key
    
    if k==x:
        res=TDD(tdd.node.succ[c])
        if tdd.node.out_weight[c]==cn0:
            res.weight=cn0
            return res
        res.weight=cn_mulCached(tdd.node.out_weight[c],tdd.weight)
        return res    
    
    if k==-1:
        res = TDD(tdd.node,tdd.weight)
        return res
    
    if k<x:
        res = TDD(tdd.node,tdd.weight)
        return res
    else:
        print("Not supported yet!!!")        
        
        

def add(tdd1,tdd2):
    global add_times
    add_times+=1
    k1=tdd1.node.key
    k2=tdd2.node.key
#     print('add 1078',k1,k2,tdd1.weight,tdd2.weight)
    if tdd1.weight==cn0:
        if tdd2.weight==cn0:
            return TDD(terminal_node,cn0)
        else:
            res = TDD(tdd2.node,getCachedComplex2(tdd2.weight.r.val,tdd2.weight.i.val))
            return res        
    
    if tdd2.weight == cn0:
        res = TDD(tdd1.node,getCachedComplex2(tdd1.weight.r.val,tdd1.weight.i.val))
        return res
    
    if tdd1.node==tdd2.node:
        weig=cn_addCached(tdd1.weight,tdd2.weight)
        if equalsZero(weig):
            releaseCached(weig)
            return TDD(terminal_node,cn0)
        else:
            res=TDD(tdd1.node,weig)
            return res
        
    if tdd1.node.id > tdd2.node.id:
        return add(tdd2,tdd1)
    
    w1 = tdd1.weight
    w2 = tdd2.weight
    
    tdd1.weight = cn1
    
    tdd2.weight = cn_divCached(tdd2.weight,w1)
        
    res = find_computed_table(['+',tdd1,tdd2])
    if res:
        tdd1.weight = w1
        releaseCached(tdd2.weight)
        tdd2.weight = w2
        if not res.weight==cn0:
            cn_mul(res.weight,res.weight,w1)
            if equalsZero(res.weight):
                releaseCached(res.weight)
                return TDD(terminal_node,cn0)
        return res
    
    the_successors=[]
    if k1>k2:
        x=k1
        for k in range(tdd1.node.succ_num):
            e1 = Slicing2(tdd1,x,k)
            e2 = tdd2
            res=add(e1,e2)
            the_successors.append(res)
            if not e1.weight==cn0:
                releaseCached(e1.weight)
    elif k1==k2:
        x=k1
        for k in range(tdd1.node.succ_num):
            e1=Slicing2(tdd1,x,k)
            e2=Slicing2(tdd2,x,k)
            res=add(e1,e2)
            the_successors.append(res)
            if not e1.weight==cn0:
                releaseCached(e1.weight)
            if not e2.weight==cn0:
                releaseCached(e2.weight)              
    else:
        x=k2
        for k in range(tdd2.node.succ_num):
            e1=tdd1
            e2=Slicing2(tdd2,x,k)
            res=add(e1,e2)
            the_successors.append(res)
            if not e2.weight==cn0:
                releaseCached(e2.weight)
                
    res = normalize(x,the_successors,True)
    insert_2_computed_table(['+',tdd1,tdd2],res)
    tdd1.weight = w1
    releaseCached(tdd2.weight)
    tdd2.weight = w2
    if not res.weight == cn0:
        cn_mul(res.weight,res.weight,w1)
    return res


def incRef(tdd):
    if tdd.node.key==-1:
        return
    
    tdd.node.ref_num+=1
    
    if tdd.node.ref_num==1:
        for k in range(tdd1.node.succ_num):
            incRef(Slicing(tdd,tdd.node.key,k))
            
def decRef(tdd):
    if tdd.node.key==-1:
        return
    
    if tdd.node.ref_num==1:
        print('Error In defRef')
    
    tdd.node.ref_num-=1
    
    if tdd.node.ref_num==0:
        for k in range(tdd1.node.succ_num):
            decRef(Slicing(tdd,tdd.node.key,k))            
    
    

def garbageCollect():
    global computed_table
    global unique_table
    temp_unique_table = dict()
    for item in unique_table:
        if not unique_table[item].ref_num==0:
            temp_unique_table[item] = unique_table[item]
    
    unique_table.clear()
    
    unique_table = temp_unique_table              
