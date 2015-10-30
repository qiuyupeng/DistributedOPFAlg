import numpy as np

import PreCache as PC

class Xupdate:
    """
    Optimization solver for both x update step
    """
    
    def __init__(self,node_set,rho):
        
        self.pCache=PC.PreCache(node_set)
        self.rho=rho
        
        
        
    def update(self,node_set):
        """
        Solve the following for each node 
        min  0.5x^H*diag(a)x+c^Hx    s.t. Bx=0
        return x
        """
        #print self.pCache.d[1]
        for i in range(len(node_set)):
            n=node_set[i]
            c=self.__construct_c(node_set,i)
            #calculate the results
            x=self.pCache.parameter_c[i]*c+self.pCache.offset[i]
            self.__writeResultBack(x,node_set,i)


            
            
        
            
    
    
    
    def __writeResultBack(self,x,node_set,nodeID):
        """
        Write the results in the x_update back to the physical variables in node_set
        """
        #print nodeID,'\n'
        for i in range(3):
            for j in range(3):
                pos_real=self.pCache.varDictReal[nodeID].me['branch_power'][i,j]
                pos_imag=self.pCache.varDictImg[nodeID].me['branch_power'][i,j]
                #oldval=node_set[nodeID].x_self_bpower[i,j]
                node_set[nodeID].x_self_bpower[i,j]=self.__writeSingleEntry(node_set[nodeID].z_bpower[i,j],pos_real,pos_imag,x)
                #print abs(oldval-node_set[nodeID].x_self_bpower[i,j])
                
                for childID in node_set[nodeID].neighbor.children:
                    pos_real=self.pCache.varDictReal[nodeID].children[childID]['branch_power'][i,j]
                    pos_imag=self.pCache.varDictImg[nodeID].children[childID]['branch_power'][i,j]
                    node_set[childID].x_parent_bpower[i,j]=self.__writeSingleEntry(node_set[childID].z_bpower[i,j],pos_real,pos_imag,x)
                
        for i in range(3):
            pos_real=self.pCache.varDictReal[nodeID].me['node_power'][i,0]
            pos_imag=self.pCache.varDictImg[nodeID].me['node_power'][i,0]
            #print pos_real, pos_imag
            node_set[nodeID].x_self_npower[i,0]=self.__writeSingleEntry(node_set[nodeID].z_npower[i,0],pos_real,pos_imag,x)

                
            
            
        parentID=node_set[nodeID].neighbor.parent
        for i in range(3):
            
            pos_real=self.pCache.varDictReal[nodeID].me['node_vol'][i,i]
            node_set[nodeID].x_self_vol[i,i]=self.__writeSingleEntry(node_set[nodeID].z_vol[i,i],pos_real,-1,x)
            pos_real=self.pCache.varDictReal[nodeID].me['branch_cur'][i,i]
            node_set[nodeID].x_self_current[i,i]=self.__writeSingleEntry(node_set[nodeID].z_current[i,i],pos_real,-1,x)
            
            if(nodeID!=0):
                pos_real=self.pCache.varDictReal[nodeID].parent['node_vol'][i,i]
                node_set[nodeID].x_parent_vol[i,i]=self.__writeSingleEntry(node_set[parentID].z_vol[i,i],pos_real,-1,x)
                #if(nodeID==2):
                #    print pos_real, node_set[nodeID].x_parent_vol[i,i],node_set[parentID].z_vol[i,i]
                
            for j in range(i):
                pos_real=self.pCache.varDictReal[nodeID].me['node_vol'][i,j]
                pos_imag=self.pCache.varDictImg[nodeID].me['node_vol'][i,j]
                node_set[nodeID].x_self_vol[i,j]=self.__writeSingleEntry(node_set[nodeID].z_vol[i,j],pos_real,pos_imag,x)
                node_set[nodeID].x_self_vol[j,i]=node_set[nodeID].x_self_vol[i,j].conjugate()
                
                pos_real=self.pCache.varDictReal[nodeID].me['branch_cur'][i,j]
                pos_imag=self.pCache.varDictImg[nodeID].me['branch_cur'][i,j]
                node_set[nodeID].x_self_current[i,j]=self.__writeSingleEntry(node_set[nodeID].z_current[i,j],pos_real,pos_imag,x)
                node_set[nodeID].x_self_current[j,i]=node_set[nodeID].x_self_current[i,j].conjugate()
                
                
                if(nodeID!=0):
                    pos_real=self.pCache.varDictReal[nodeID].parent['node_vol'][i,j]
                    pos_imag=self.pCache.varDictImg[nodeID].parent['node_vol'][i,j]
                    node_set[nodeID].x_parent_vol[i,j]=self.__writeSingleEntry(node_set[parentID].z_vol[i,j],pos_real,pos_imag,x)
                    node_set[nodeID].x_parent_vol[j,i]=node_set[nodeID].x_parent_vol[i,j].conjugate()
            
            for j in range(i+1):
                for childID in node_set[nodeID].neighbor.children:
                    if(i==j):
                        pos_real=self.pCache.varDictReal[nodeID].children[childID]['branch_cur'][i,i]
                        node_set[childID].x_parent_current[i,i]=self.__writeSingleEntry(node_set[childID].z_current[i,i],pos_real,-1,x)
                    else:
                        pos_real=self.pCache.varDictReal[nodeID].children[childID]['branch_cur'][i,j]
                        pos_imag=self.pCache.varDictImg[nodeID].children[childID]['branch_cur'][i,j]
                        node_set[childID].x_parent_current[i,j]=self.__writeSingleEntry(node_set[childID].z_current[i,j],pos_real,pos_imag,x)
                        node_set[childID].x_parent_current[j,i]=node_set[childID].x_parent_current[i,j].conjugate()

                
                        
                        
                   
                
                
                
            
            
    
    def __writeSingleEntry(self,originalVal,pos_real,pos_imag,x):
        """
        Return a new update if either pos_real or pos_imag are not -1
        """
        if(pos_real>-0.5 and pos_imag>-0.5):
            return complex(x[pos_real,0],x[pos_imag,0])
        elif(pos_real>-0.5 and pos_imag<-0.5):
            return complex(x[pos_real,0],originalVal.imag)
        elif(pos_real<-0.5 and pos_imag>-0.5):
            return complex(originalVal.real,x[pos_imag,0])
        else:
            return originalVal
        
    
    
    def __construct_c(self,node_set,nodeID):
        """
        construct the coefficient c in the following optimzation problem:
        min 0.5x^TAx+c^Tx   s.t. Bx=d
        """
        
        n=node_set[nodeID]
        num_child=len(n.neighbor.children)
        c=np.matrix(np.zeros([self.pCache.numVar[nodeID],1]))
        
        #branch power S_i and S_j j \in C_i
        for i in range(3):
            for j in range(3):
                pos=self.pCache.varDictReal[nodeID].me['branch_power'][i,j]
                if(pos>-0.5):
                    c[pos,0]=-(2*num_child+3)*n.z_bpower[i,j].real+n.lambda_bpower[i,j].real/self.rho
                pos=self.pCache.varDictImg[nodeID].me['branch_power'][i,j]
                if(pos>-0.5):
                    c[pos,0]=-(2*num_child+3)*n.z_bpower[i,j].imag+n.lambda_bpower[i,j].imag/self.rho
                
                for childID in n.neighbor.children:
                    child=node_set[childID]
                    pos=self.pCache.varDictReal[nodeID].children[childID]['branch_power'][i,j]
                    if(pos>-0.5):
                        c[pos,0]=-child.z_bpower[i,j].real+child.mu_bpower[i,j].real/self.rho
                    pos=self.pCache.varDictImg[nodeID].children[childID]['branch_power'][i,j]
                    if(pos>-0.5):
                        c[pos,0]=-child.z_bpower[i,j].imag+child.mu_bpower[i,j].imag/self.rho
        
            
        
        #power injection
        for i in range(3):
            pos=self.pCache.varDictReal[nodeID].me['node_power'][i,0]
            if(pos>-0.5):
                c[pos,0]=(-n.z_npower[i,0].real+n.lambda_npower[i,0].real/self.rho)
            pos=self.pCache.varDictImg[nodeID].me['node_power'][i,0]
            if(pos>-0.5):
                c[pos,0]=(-n.z_npower[i,0].imag+n.lambda_npower[i,0].imag/self.rho)
        
        #voltage, current (hermittian and only lower triangular is required)
        if(n.seq_num!=0):
            p=node_set[n.neighbor.parent]
        for i in range(3):
            for j in range(i+1):
                
                isdiag=1 if (i==j) else 2
                
                pos=self.pCache.varDictReal[nodeID].me['node_vol'][i,j]
                if(pos>-0.5):
                    c[pos,0]=-2*n.z_vol[i,j].real*isdiag+n.lambda_vol[i,j].real/self.rho
                pos=self.pCache.varDictImg[nodeID].me['node_vol'][i,j]
                if(pos>-0.5):
                    c[pos,0]=-2*n.z_vol[i,j].imag*isdiag+n.lambda_vol[i,j].imag/self.rho
                
                pos=self.pCache.varDictReal[nodeID].me['branch_cur'][i,j]
                if(pos>-0.5):
                    c[pos,0]=-(num_child+1)*n.z_current[i,j].real*isdiag+n.lambda_current[i,j].real/self.rho
                pos=self.pCache.varDictImg[nodeID].me['branch_cur'][i,j]
                if(pos>-0.5):
                    c[pos,0]=-(num_child+1)*n.z_current[i,j].imag*isdiag+n.lambda_current[i,j].imag/self.rho
                    
                pos=self.pCache.varDictReal[nodeID].parent['node_vol'][i,j]
                if(pos>-0.5):
                    c[pos,0]=-p.z_vol[i,j].real*isdiag+n.gamma_vol[i,j].real/self.rho
                    #if(n.seq_num==2):
                    #    print pos, p.z_vol[i,j], c[pos,0],self.pCache.A[2][pos,pos]
                pos=self.pCache.varDictImg[nodeID].parent['node_vol'][i,j]
                if(pos>-0.5):
                    c[pos,0]=-p.z_vol[i,j].imag*isdiag+n.gamma_vol[i,j].imag/self.rho
                    
                for childID in n.neighbor.children:
                    child=node_set[childID]
                    pos=self.pCache.varDictReal[nodeID].children[childID]['branch_cur'][i,j]
                    if(pos>-0.5):
                        c[pos,0]=-child.z_current[i,j].real*isdiag+child.mu_current[i,j].real/self.rho
                    pos=self.pCache.varDictImg[nodeID].children[childID]['branch_cur'][i,j]
                    if(pos>-0.5):
                        c[pos,0]=-child.z_current[i,j].imag*isdiag+child.mu_current[i,j].imag/self.rho
        
        return c
                
                
            
        
        
        


            
            
                
            
            

        
        

        
        
    
        
        
        
        
    