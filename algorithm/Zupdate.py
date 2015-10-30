import numpy as np
import scipy.linalg as LA
import math
from copy import deepcopy

class Zupdate:
    """
    Optimization solver for z update step
    
    member variables:
    self.rho
    self.phaseDict
    self.numPhase
    self.numChildren
    """
    
    def __init__(self,node_set,rho):
        
        self.rho=rho
        self.phaseDict=[]
        self.numPhase=[]
        self.numChildren=[]
        self.numBus=len(node_set)
        self.Xvar=[None]*self.numBus
        
        #the parameter in front of p in the objective
        self.objPara=np.matrix(np.ones([self.numBus,3],dtype=float))*0.1
        #self.objPara=np.matrix(np.zeros([self.numBus,3],dtype=float))*0.1
        
        
        self.__construct_phaseDict(node_set)
        

        
        
    def update(self,node_set):
        """
        update the z variables
        
        return the dual residual
        """
        
        
        dualRes=0
        dualRes+=self.__update_s(node_set)
        dualRes+=self.__update_v_ell_S(node_set)
        
        #print node_set[1].z_vol
        #print oldval.z_current-node_set[1].z_current
        #print LA.norm(oldval.z_bpower-node_set[1].z_bpower)
        #print LA.norm(oldval.z_npower-node_set[1].z_npower)
        
        return math.sqrt(dualRes)
        
        
        
    def __update_s(self,node_set):
        """
        Update the node power injection
        """
        dualRes_square=0
        for n in node_set:
            for j in range(3):
                oldVal=n.z_npower[j,0]
                if(n.phase[j]==1):
                    if(n.control_info[j]==1 or n.control_info[j]==3):
                        tmp_upper=n.fixed_load[j].real+n.upper_control[j].real
                        tmp_lower=n.fixed_load[j].real+n.lower_control[j].real
                        tmp_update=n.x_self_npower[j,0].real+(n.lambda_npower[j,0].real-self.objPara[n.seq_num,j])/self.rho
                        n.z_npower[j,0]=complex(self.__band(tmp_upper,tmp_lower,tmp_update),n.z_npower[j].imag)
                    if(n.control_info[j]==2 or n.control_info[j]==3):
                        tmp_upper=n.fixed_load[j].imag+n.upper_control[j].imag
                        tmp_lower=n.fixed_load[j].imag+n.lower_control[j].imag
                        tmp_update=n.x_self_npower[j,0].imag+(n.lambda_npower[j,0].imag)/self.rho

                        n.z_npower[j,0]=complex(n.z_npower[j,0].real,self.__band(tmp_upper,tmp_lower,tmp_update))
                        

                        
                dualRes_square+=LA.norm(oldVal-n.z_npower[j,0])**2
        return dualRes_square
                        
                    
    def __band(self,upper,lower,a):
        return max(lower,min(a,upper))
        
    def __update_v_ell_S(self,node_set):
        """
        update v_z, ell_z, S_z by solving 
        min ||W-Xvar||_2^2 s.t. W>=0
        """
        self.__construct_Xvar(node_set)
        
        
        
        dualRes_square=0
        for i in range(self.numBus):
            n=node_set[i]
            if n.seq_num!=0:
                w,v=LA.eigh(self.Xvar[i])
                v=np.matrix(v)
                W=np.matrix(np.zeros([self.numPhase[i]*2,self.numPhase[i]*2],dtype=complex))
                
                #print w
                if(w.min()>=0):
                    W=self.Xvar[i]
                else:
                    for j in range(self.numPhase[i]*2):
                        w[j]=max(w[j],0)
                        W+=w[j]*v[:,j]*(v[:,j].getH())
                
                dualRes_square+=self.__write_v_ell_S_back(W,node_set,i)                

        return dualRes_square
        
        
                    
        
        
    def __write_v_ell_S_back(self,W,node_set,nodeID):
        """
        Write the results in the x_update back to the physical variables (v,ell,S only) in node_set
        """
        dualRes_square=0
        n=node_set[nodeID]
        
        for i in range(3):
            for j in range(3):
                pos_1=self.phaseDict[nodeID][i,0]
                pos_2=self.phaseDict[nodeID][j,0]
                if(pos_1>-0.5 and pos_2>-0.5):
                    
                    dualRes_square+=LA.norm(n.z_vol[i,j]-W[pos_1,pos_2])**2
                    #print LA.norm(n.z_vol[i,j]-W[pos_1,pos_2])**2
                    dualRes_square+=LA.norm(n.z_current[i,j]-W[pos_1+self.numPhase[nodeID],pos_2+self.numPhase[nodeID]])**2
                    #print LA.norm(n.z_current[i,j]-W[pos_1+self.numPhase[nodeID],pos_2+self.numPhase[nodeID]])**2

                    dualRes_square+=(LA.norm(n.z_bpower[i,j]-W[pos_1,pos_2+self.numPhase[nodeID]])**2)*2
                    #print (LA.norm(n.z_bpower[i,j]-W[pos_1,pos_2+self.numPhase[nodeID]])**2)*2
                    
                    #update v_z
                    n.z_vol[i,j]=W[pos_1,pos_2]
                    #update current
                    n.z_current[i,j]=W[pos_1+self.numPhase[nodeID],pos_2+self.numPhase[nodeID]]
                    #update branch power
                    n.z_bpower[i,j]=W[pos_1,pos_2+self.numPhase[nodeID]]          
        """
        WW=np.matrix(np.zeros([6,6],dtype=complex))
        WW[0:3,0:3]=n.z_vol
        WW[3:6,3:6]=n.z_current
        WW[0:3,3:6]=n.z_bpower
        WW[3:6,0:3]=n.z_bpower.getH()
        vv,ww=LA.eigh(WW)
        print vv
        
        print n.z_vol[0,0].real*n.z_current[0,0].real-abs(n.z_bpower[0,0])**2
        """
        return dualRes_square
                
        
    
    def __construct_Xvar(self,node_set):
        """
        min ||Z-X||_2^2  s.t. Z>=0
        Construct the X
        """
        
        for i in range(1,self.numBus):
            n=node_set[i]
            self.Xvar[i]=np.matrix(np.zeros([self.numPhase[i]*2,self.numPhase[i]*2],dtype=complex))
            for j in range(3):     
                for k in range(3):
                    pos_1=self.phaseDict[i][j,0]
                    pos_2=self.phaseDict[i][k,0]
                    if(pos_1>-0.5 and pos_2>-0.5):
                        
                        if j==k:
                            #voltage itself
                            self.Xvar[i][pos_1,pos_2]=n.x_self_vol[j,k]*2+n.lambda_vol[j,k]/self.rho
                            
                            #current itself and parent
                            self.Xvar[i][pos_1+self.numPhase[i],pos_2+self.numPhase[i]]=n.x_self_current[j,k]*(self.numChildren[i]+1)+n.lambda_current[j,k]/self.rho
                            self.Xvar[i][pos_1+self.numPhase[i],pos_2+self.numPhase[i]]+=(n.x_parent_current[j,k]+n.mu_current[j,k]/self.rho)
                        else:
                            #voltage itself
                            self.Xvar[i][pos_1,pos_2]=n.x_self_vol[j,k]*2+n.lambda_vol[j,k]/(2*self.rho)
                            
                            #current itself and parent
                            self.Xvar[i][pos_1+self.numPhase[i],pos_2+self.numPhase[i]]=n.x_self_current[j,k]*(self.numChildren[i]+1)+n.lambda_current[j,k]/(2*self.rho)
                            self.Xvar[i][pos_1+self.numPhase[i],pos_2+self.numPhase[i]]+=(n.x_parent_current[j,k]+n.mu_current[j,k]/(2*self.rho))
                            
                        
                        #branch power itself and parent
                        self.Xvar[i][pos_1,pos_2+self.numPhase[i]]=n.x_self_bpower[j,k]*(self.numChildren[i]+1.5)+n.lambda_bpower[j,k]/(2*self.rho)
                        self.Xvar[i][pos_1,pos_2+self.numPhase[i]]+=(n.x_parent_bpower[j,k]*0.5+n.mu_bpower[j,k]/(2*self.rho))
                        #branch power symmetric part
                        self.Xvar[i][pos_2+self.numPhase[i],pos_1]=self.Xvar[i][pos_1,pos_2+self.numPhase[i]].conjugate()
                    #voltage duplicates from children
                    if(n.seq_num!=0):
                        for childID in n.neighbor.children:
                            c=node_set[childID]
                            if(pos_1>-0.5 and pos_2>-0.5):
                                if j==k:
                                    self.Xvar[i][pos_1,pos_2]+=(c.x_parent_vol[j,k]+c.gamma_vol[j,k]/self.rho)
                                else:
                                    self.Xvar[i][pos_1,pos_2]+=(c.x_parent_vol[j,k]+c.gamma_vol[j,k]/(2*self.rho))
                   
            self.Xvar[i]/=(self.numChildren[i]+2)
            
            #if(LA.norm(self.Xvar[i]-self.Xvar[i].getH())>1e-10):
            #    print "The hermittian matrix at bus", i, "is not correct!"
        
        
                
    
    def __construct_phaseDict(self,node_set):
        """
        Construct the phase dictionary
        """
        
        for n in node_set:
            n_dict=-np.matrix(np.ones([3,1],dtype=int))
            num_phase=0
            for i in range(3):
                if(n.phase[i]==1):
                    n_dict[i]=num_phase
                    num_phase+=1
            self.phaseDict.append(n_dict)
            self.numPhase.append(num_phase)
            self.numChildren.append(len(n.neighbor.children))
            
            
            
            
            
            
            
            
            
            