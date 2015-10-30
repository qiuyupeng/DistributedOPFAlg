import numpy as np

class VarDict:
    
    def __init__(self,n):
        """
        -1 means that entry is not a variable
        """
        self.parent=dict([])
        self.parent['node_vol']=-1*np.ones([3,3],dtype=int)
        
        self.me=dict([])
        self.me['node_vol']=-1*np.ones([3,3],dtype=int)
        self.me['node_power']=-1*np.ones([3,1],dtype=int)
        self.me['branch_cur']=-1*np.ones([3,3],dtype=int)
        self.me['branch_power']=-1*np.ones([3,3],dtype=int)
        
        self.children=dict([])
        for childID in n.neighbor.children:
            self.children[childID]=dict([])
            self.children[childID]['branch_cur']=-1*np.ones([3,3],dtype=int)
            self.children[childID]['branch_power']=-1*np.ones([3,3],dtype=int)
        
