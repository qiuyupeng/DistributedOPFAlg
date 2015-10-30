import numpy as np
import Constant as ct
import time
#import matplotlib.pyplot as plt


from Xupdate import *
from Zupdate import *
from Mupdate import *

class Solver:
    """
    Optimization solver for both x and z update step
    """
    
    def __init__(self,node_set):
        
        self.rho=ct.Constant().rho
        self.node_set=node_set
        self.xupdate=Xupdate(self.node_set,self.rho)
        self.zupdate=Zupdate(self.node_set,self.rho)
        self.mupdate=Mupdate(self.node_set,self.rho)
        
    
    def start(self,num_iteration=None):
        """
        Do xupdate, zupdate, multipliers update sequentially
        if num_iteration=None, using stopping criteria
        Else, using fixed num_iteration
        """
        primalRes=[]
        dualRes=[]
        exe_time=[]
        if(num_iteration!=None):
            for i in range(num_iteration):
                print "Iteration: ", i,
                start_time=time.time()
                self.xupdate.update(self.node_set)
                dualRes.append(self.zupdate.update(self.node_set))
                primalRes.append(self.mupdate.update(self.node_set))
                exe_time.append(time.time()-start_time)
                print " Primal: ", primalRes[i], " Dual: ",dualRes[i], " Time: ", exe_time[i]
        
        
        #plot the primal and dual residual v.s. iteration
        #plt.semilogy(range(num_iteration),primalRes)
        #plt.semilogy(range(num_iteration),dualRes)
        return (primalRes,dualRes,exe_time)
        
        
        
        
    

            
            
                
            
            

        
        

        
        
    
        
        
        
        
    