import numpy as np


class Node:
    """
    contain all the information about one bus(node)
    power is in MW/MVAR/MVA
    """
    def __init__(self,num):
        
        self.seq_num=num
        self.label=1   #its label in the IEEE 
        
        self.phase=np.array([0,0,0],dtype=int)
        self.control_info=np.array([0,0,0],dtype=int) # 0: no control, 1: real power, 2:reactive power control, 3: both real and reactive
        self.fixed_load=np.array([0,0,0],dtype=complex)
        self.lower_control=np.array([0,0,0],dtype=complex)
        self.upper_control=np.array([0,0,0],dtype=complex)
        
        self.neighbor=Neighbor()
        self.impedance=None # line impedance connected from this node to ancestor
        
        #Variables below are for computation only
        
        #primal variables
        self.x_self_vol=np.matrix(np.zeros([3,3],dtype=complex))
        self.x_self_current=np.matrix(np.zeros([3,3],dtype=complex))
        self.x_self_bpower=np.matrix(np.zeros([3,3],dtype=complex))
        self.x_self_npower=np.matrix(np.zeros([3,1],dtype=complex))
        
        self.x_parent_vol=np.matrix(np.zeros([3,3],dtype=complex))
        self.x_parent_current=np.matrix(np.zeros([3,3],dtype=complex))
        self.x_parent_bpower=np.matrix(np.zeros([3,3],dtype=complex))
        
        self.z_vol=np.matrix(np.zeros([3,3],dtype=complex))
        self.z_current=np.matrix(np.zeros([3,3],dtype=complex))
        self.z_bpower=np.matrix(np.zeros([3,3],dtype=complex))
        self.z_npower=np.matrix(np.zeros([3,1],dtype=complex))
        
        #multipliers
        self.lambda_vol=np.matrix(np.zeros([3,3],dtype=complex))
        self.lambda_current=np.matrix(np.zeros([3,3],dtype=complex))
        self.lambda_bpower=np.matrix(np.zeros([3,3],dtype=complex))
        self.lambda_npower=np.matrix(np.zeros([3,1],dtype=complex))
        
        self.gamma_vol=np.matrix(np.zeros([3,3],dtype=complex))
        self.mu_current=np.matrix(np.zeros([3,3],dtype=complex))
        self.mu_bpower=np.matrix(np.zeros([3,3],dtype=complex))
        
        
        #backward,forward sweep variables
        self.node_vol=np.matrix(np.zeros((3,1))*(1+0j))
        self.node_cur=np.matrix(np.zeros((3,1))*(1+0j))
        self.branch_cur=np.matrix(np.zeros((3,1))*(1+0j))
        
        
        
        
        
        
        self.load_control=np.array([0,0,0],dtype=complex)
        
        
    def __str__(self):
        s='****************************\n'+\
        'Label: '+str(self.label)+'\n'+\
        'Sequence number: '+str(self.seq_num)+'\n'+\
        'Phase: '+str(self.phase)+'\n'+\
        'Fixed load: '+str(self.fixed_load)+'\n'+\
        'Control upper: '+str(self.upper_control)+'\n'+\
        'x voltage: '+str(self.x_self_vol)+'\n'+\
        'x current: '+str(self.x_self_current)+'\n'+\
        'x branch power: '+str(self.x_self_bpower)+'\n'+\
        'x node power: '+str(self.x_self_npower)+'\n'+\
        'x parent voltage: '+str(self.x_parent_vol)+'\n'+\
        'x parent current: '+str(self.x_parent_current)+'\n'+\
        'x parent branch power: '+str(self.x_parent_bpower)+'\n'+\
        'z voltage: '+str(self.z_vol)+'\n'+\
        'z current: '+str(self.z_current)+'\n'+\
        'z branch power: '+str(self.z_bpower)+'\n'+\
        'z node power: '+str(self.z_npower)+'\n'+\
        '****************************\n'
        return s
        
   
   
   
   
   
        
class Neighbor:
    def __init__(self):
        self.parent=None
        self.children=set([])
        self.neighbor=set([])
        
    def __str__(self):
        s="Parent: "+repr(self.parent)+'.\n Children: '+repr(self.children)+'.\n Neighbor: '+repr(self.neighbor) 
        return s
        
if __name__=='__main__':
    n=Node()