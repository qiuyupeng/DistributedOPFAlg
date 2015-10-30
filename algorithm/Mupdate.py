import numpy as np
import scipy.linalg as LA
import math



class Mupdate:
    """
    update all the lagrangian multipliers
    """
    
    def __init__(self,node_set,rho):
        
        self.rho=rho
    
    
    
    def update(self,node_set):
        """
        lambda=lamda+rho(x-z)
        
        return primal residual
        """
        primalRes=0
        for n in node_set:
            
            numChildren=len(n.neighbor.children)
            
            n.lambda_bpower=n.lambda_bpower+self.rho*(n.x_self_bpower-n.z_bpower)*(2*numChildren+3)
            n.lambda_npower=n.lambda_npower+self.rho*(n.x_self_npower-n.z_npower)
            
            for i in range(3):
                n.lambda_vol[i,i]=n.lambda_vol[i,i]+self.rho*(n.x_self_vol[i,i]-n.z_vol[i,i])*2
                n.lambda_current[i,i]=n.lambda_current[i,i]+self.rho*(n.x_self_current[i,i]-n.z_current[i,i])*(numChildren+1)
                for j in range(i):
                    n.lambda_vol[i,j]=n.lambda_vol[i,j]+self.rho*(n.x_self_vol[i,j]-n.z_vol[i,j])*2*2
                    n.lambda_vol[j,i]=n.lambda_vol[i,j].conjugate()
                    n.lambda_current[i,j]=n.lambda_current[i,j]+self.rho*(n.x_self_current[i,j]-n.z_current[i,j])*(numChildren+1)*2
                    n.lambda_current[j,i]=n.lambda_current[i,j].conjugate()
                    
            #print n.seq_num
            primalRes+=LA.norm(n.x_self_vol-n.z_vol)**2
            #print primalRes
            primalRes+=LA.norm(n.x_self_current-n.z_current)**2
            #print primalRes
            primalRes+=LA.norm(n.x_self_bpower-n.z_bpower)**2
            #print primalRes
            primalRes+=LA.norm(n.x_self_npower-n.z_npower)**2
            #print primalRes
            
            #if n.seq_num==0:
            #    print n.x_self_npower,'\n',n.z_npower
            
            
            if n.seq_num!=0:
                p=node_set[n.neighbor.parent]
                n.mu_bpower=n.mu_bpower+self.rho*(n.x_parent_bpower-n.z_bpower)
                for i in range(3):
                    n.mu_current[i,i]=n.mu_current[i,i]+self.rho*(n.x_parent_current[i,i]-n.z_current[i,i])
                    n.gamma_vol[i,i]=n.gamma_vol[i,i]+self.rho*(n.x_parent_vol[i,i]-p.z_vol[i,i])
                    for j in range(i):
                        n.mu_current[i,j]=n.mu_current[i,j]+self.rho*(n.x_parent_current[i,j]-n.z_current[i,j])*2
                        n.mu_current[j,i]=n.mu_current[i,j].conjugate()
                        n.gamma_vol[i,j]=n.gamma_vol[i,j]+self.rho*(n.x_parent_vol[i,j]-p.z_vol[i,j])*2
                        n.gamma_vol[j,i]=n.gamma_vol[i,j].conjugate()

                
                
                primalRes+=LA.norm(n.x_parent_current-n.z_current)**2
                #print primalRes
                primalRes+=LA.norm(n.x_parent_bpower-n.z_bpower)**2
                #print primalRes
                primalRes+=LA.norm(n.x_parent_vol-p.z_vol)**2
                #print primalRes
        
        return math.sqrt(primalRes)
            
            