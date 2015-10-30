import cmath
import math
import ForwardBackSweep
import numpy as np
from copy import deepcopy

import scipy.linalg as LA

class Initialization:
    """
    Initialize all the primal and dual variables
    """
    
    def __init__(self,node_set):
        
        FBS=ForwardBackSweep.ForwardBackSweep()
        FBS.solver(node_set)
        
        #Intilization for ADMM
        self.__DFS_init(node_set,0)

            
    
    def __DFS_init(self,node_set,nodeID):
        """
        Initialize using depth first search
        """
        
        n=node_set[nodeID]
        n.x_self_vol=n.node_vol*n.node_vol.getH()
        n.x_self_current=n.branch_cur*n.branch_cur.getH()
        n.x_self_bpower=n.node_vol*n.branch_cur.getH()
        if(nodeID==0):
            n.x_self_npower=np.matrix([[-n.node_vol[i,0]*n.branch_cur[i,0].conjugate()] for i in range(3)])
            #n.x_self_npower=np.matrix(np.zeros([3,1],dtype=complex))
        else:
            n.x_self_npower=np.matrix([[x] for x in n.fixed_load],dtype=complex)
            
        
        
        

        
        """
        W=np.matrix(np.zeros([6,6],dtype=complex))
        W[0:3,0:3]=node_set[nodeID].x_self_vol
        W[0:3,3:6]=node_set[nodeID].x_self_bpower
        W[3:6,0:3]=node_set[nodeID].x_self_bpower.getH()
        W[3:6,3:6]=node_set[nodeID].x_self_current
        w,v=LA.eigh(W)
        print w
        """
        
        
        n.z_vol=deepcopy(n.x_self_vol)
        n.z_current=deepcopy(n.x_self_current)
        n.z_bpower=deepcopy(n.x_self_bpower)
        n.z_npower=deepcopy(n.x_self_npower)
        
        n.x_parent_current=deepcopy(n.x_self_current)
        n.x_parent_bpower=deepcopy(n.x_self_bpower)
        if(nodeID!=0):
            n.x_parent_vol=deepcopy(node_set[n.neighbor.parent].x_self_vol)
            
        
        for i in n.neighbor.children:
            self.__DFS_init(node_set,i)
        
        

        
        
        
        
        

        
        
        
        
        
        
        
        
        



            
            
            
            
            
            
            
            
            