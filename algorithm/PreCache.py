import numpy as np
import numpy.linalg as LA
from VarDict import *

class PreCache:
    """
    Calculate the iteration independent variables.

    Member variables:
    self.varDictReal
    self.varDictImg
    self.numVar : array of number of variables associated with each bus

    self.numBus
    
    
    min 0.5x^TAx+c^Tx   s.t. Bx=d
    
    self.A
    self.B
    self.d
    
    Solution: x=parameter_c*c+offset
    
    self.parameter_c
    self.offset
    """

    def __init__(self,node_set):

        self.varDictReal=[]
        self.varDictImg=[]
        self.numVar=[]
        self.numBus=len(node_set)
        
        self.A=[None]*self.numBus
        self.B=[None]*self.numBus
        self.d=[None]*self.numBus
        
        self.parameter_c=[None]*self.numBus
        self.offset=[None]*self.numBus

        self.__hash_pos(node_set)
        self.__construct_Matrix(node_set)
        self.__construct_parameter(node_set)
        
        
        
    def __construct_parameter(self,node_set):
        
        
        for i in range(self.numBus):
            tmp0=LA.inv(self.A[i])
            
            tmp1=LA.inv(self.B[i]*tmp0*self.B[i].getT())
            
            self.parameter_c[i]=tmp0*self.B[i].getT()*tmp1*self.B[i]*tmp0-tmp0
            self.offset[i]=tmp0*self.B[i].getT()*tmp1*self.d[i]
        
    
    def __construct_Matrix(self,node_set):
        """
        Construct A,B,d
        """
        self.__construct_A(node_set)
        self.__construct_B_d(node_set)
        

    def __construct_B_d(self,node_set):
        """
        Construct the B and d in the constraint
        """
        
        for i in range(self.numBus):
            n=node_set[i]
            B=np.matrix(np.zeros([15,self.numVar[i]]))
            d=np.matrix(np.zeros([15,1]))
            
            
            #diag(sum_{j in C_i}(S_j-z_jl_j)-S_i)+s_i=0
            for j in range(3):
                
                pos=self.varDictReal[i].me['branch_power'][j,j]
                if(pos>-0.5):
                    B[9+j,pos]=-1

                pos=self.varDictImg[i].me['branch_power'][j,j]
                if(pos>-0.5):
                    B[12+j,pos]=-1
                
                pos=self.varDictReal[i].me['node_power'][j,0]
                
                if(pos>-0.5):
                    B[9+j,pos]=1
                else:
                    d[9+j,0]=-n.fixed_load[j].real
                
                pos=self.varDictImg[i].me['node_power'][j,0]
                if(pos>-0.5):
                    B[12+j,pos]=1
                else:
                    d[12+j,0]=-n.fixed_load[j].imag
                
                for childID in n.neighbor.children:
                    c=node_set[childID]
                    pos=self.varDictReal[i].children[childID]['branch_power'][j,j]
                    if(pos>-0.5):
                        B[9+j,pos]=1
                    pos=self.varDictImg[i].children[childID]['branch_power'][j,j]
                    if(pos>-0.5):
                        B[12+j,pos]=1
                        
                    pos=self.varDictReal[i].children[childID]['branch_cur'][0,0]
                    #B[9+j,pos]+=(-c.impedance[0,0].real)
                    #B[12+j,pos]+=(-c.impedance[0,0].imag)
                    
                    for k in range(3):
                        if(k<j):
                            pos=self.varDictReal[i].children[childID]['branch_cur'][j,k]
                            if(pos>-0.5):
                                B[9+j,pos]+=(-c.impedance[j,k].real)
                                B[12+j,pos]+=-(c.impedance[j,k].imag)
                            pos=self.varDictImg[i].children[childID]['branch_cur'][j,k]
                            if(pos>-0.5):
                                B[9+j,pos]+=(-c.impedance[j,k].imag)
                                B[12+j,pos]+=(c.impedance[j,k].real)
                        elif(k==j):
                            pos=self.varDictReal[i].children[childID]['branch_cur'][j,j]
                            if(pos>-0.5):
                                B[9+j,pos]+=(-c.impedance[j,j].real)
                                B[12+j,pos]+=(-c.impedance[j,j].imag)
                        else:
                            pos=self.varDictReal[i].children[childID]['branch_cur'][k,j]
                            if(pos>-0.5):
                                B[9+j,pos]+=(-c.impedance[j,k].real)
                                B[12+j,pos]+=(-c.impedance[j,k].imag)
                            pos=self.varDictImg[i].children[childID]['branch_cur'][k,j]
                            if(pos>-0.5):
                                B[9+j,pos]+=(c.impedance[j,k].imag)
                                B[12+j,pos]+=(-c.impedance[j,k].real)
                    
                    
            
            #v_{A_i}-v_i+zS^*+Sz^*-zlz=0
            if(n.seq_num!=0):
                Z=np.ndarray(shape=(3,3,3,3),dtype=complex)
                for j1 in range(3):
                    for j2 in range(3):
                        for j3 in range(3):
                            for j4 in range(3):
                                Z[j1,j2,j3,j4]=n.impedance[j1,j2]*(n.impedance[j3,j4].conjugate())
            
            
            for j in range(3):
                for k in range(j+1):
                    pos=self.varDictReal[i].parent['node_vol'][j,k]
                    if(n.phase[j]==1 and n.phase[k]==1 and pos>-0.5):
                        B[(j*(j+1))/2+k,pos]=1
                    pos=self.varDictImg[i].parent['node_vol'][j,k]
                    if(n.phase[j]==1 and n.phase[k]==1 and pos>-0.5):
                        B[6+(j*(j-1))/2+k,pos]=1
                    
                    pos=self.varDictReal[i].me['node_vol'][j,k]
                    if(pos>-0.5):
                        B[(j*(j+1))/2+k,pos]=-1
                    pos=self.varDictImg[i].me['node_vol'][j,k]
                    if(pos>-0.5):
                        B[6+(j*(j-1))/2+k,pos]=-1
                    
                    #zS^H+Sz^H
                    for l in range(3):
                        pos=self.varDictReal[i].me['branch_power'][k,l]
                        if(pos>-0.5):
                            B[(j*(j+1))/2+k,pos]+=n.impedance[j,l].real
                        pos=self.varDictImg[i].me['branch_power'][k,l]
                        if(pos>-0.5):
                            B[(j*(j+1))/2+k,pos]+=n.impedance[j,l].imag
                            
                        pos=self.varDictReal[i].me['branch_power'][j,l]
                        if(pos>-0.5):
                            B[(j*(j+1))/2+k,pos]+=n.impedance[k,l].real
                        pos=self.varDictImg[i].me['branch_power'][j,l]
                        if(pos>-0.5):
                            B[(j*(j+1))/2+k,pos]+=n.impedance[k,l].imag
                        
                        if(j!=k):
                            pos=self.varDictReal[i].me['branch_power'][k,l]
                            if(pos>-0.5):
                                B[6+(j*(j-1))/2+k,pos]+=n.impedance[j,l].imag
                            pos=self.varDictImg[i].me['branch_power'][k,l]
                            if(pos>-0.5):
                                B[6+(j*(j-1))/2+k,pos]+=(-n.impedance[j,l].real)
                            
                            pos=self.varDictReal[i].me['branch_power'][j,l]
                            if(pos>-0.5):
                                B[6+(j*(j-1))/2+k,pos]+=(-n.impedance[k,l].imag)
                            pos=self.varDictImg[i].me['branch_power'][j,l]
                            if(pos>-0.5):
                                B[6+(j*(j-1))/2+k,pos]+=(n.impedance[k,l].real)
                    
                    #zlz
                    
                    for l1 in range(3):
                        pos=self.varDictReal[i].me['branch_cur'][l1,l1]
                        if(pos>-0.5):
                            B[(j*(j+1))/2+k,pos]+=(-Z[j,l1,k,l1].real)
                            if(j!=k):
                                B[6+(j*(j-1))/2+k,pos]+=(-Z[j,l1,k,l1].imag)

                        for l2 in range(l1):
                            pos=self.varDictReal[i].me['branch_cur'][l1,l2]
                            if(pos>-0.5):
                                B[(j*(j+1))/2+k,pos]+=(-Z[j,l1,k,l2].real)
                                B[(j*(j+1))/2+k,pos]+=(-Z[j,l2,k,l1].real)
                                if(j!=k):
                                    B[6+(j*(j-1))/2+k,pos]+=(-Z[j,l1,k,l2].imag)
                                    B[6+(j*(j-1))/2+k,pos]+=(-Z[j,l2,k,l1].imag)

                            pos=self.varDictImg[i].me['branch_cur'][l1,l2]
                            if(pos>-0.5):
                                B[(j*(j+1))/2+k,pos]-=(-Z[j,l1,k,l2].imag)
                                B[(j*(j+1))/2+k,pos]+=(-Z[j,l2,k,l1].imag)
                                if(j!=k):
                                    B[6+(j*(j-1))/2+k,pos]+=(-Z[j,l1,k,l2].real)
                                    B[6+(j*(j-1))/2+k,pos]-=(-Z[j,l2,k,l1].real)
                        
                                


                            
                            
                    
            
            #remove the all zero rows

            self.B[i]=np.matrix(np.zeros([0,self.numVar[i]],dtype=complex))
            self.d[i]=np.matrix(np.zeros([0,1],dtype=complex))
            for j in range(15):
                if(np.sum(abs(B[j]))>1e-10):
                    self.B[i]=np.vstack((self.B[i],B[j]))
                    self.d[i]=np.vstack((self.d[i],d[j]))                  
    
                    
            
            
            
    
                    
    
    def __construct_A(self,node_set):
        """
        Construct the A in the objective
        """
        
        for i in range(self.numBus):
            n=node_set[i]
            A=np.matrix(np.zeros([self.numVar[i],self.numVar[i]]))
            
            
            num_children=len(n.neighbor.children)
            
            #update the coefficient in front of S_i, S_j j\in C_i
            for j in range(3):
                for k in range(3):
                    pos=self.varDictReal[i].me['branch_power'][j,k]
                    if(pos>-0.5):
                        A[pos,pos]=2*num_children+3.0
                        
                    pos=self.varDictImg[i].me['branch_power'][j,k]
                    if(pos>-0.5):
                        A[pos,pos]=2*num_children+3.0
                    
                    for childID in n.neighbor.children:
                        pos=self.varDictReal[i].children[childID]['branch_power'][j,k]
                        if(pos>-0.5):
                            A[pos,pos]=1
                        pos=self.varDictImg[i].children[childID]['branch_power'][j,k]
                        if(pos>-0.5):
                            A[pos,pos]=1

            
            #update the coefficient in front of s_i(power injection)
            for j in range(3):
                pos=self.varDictReal[i].me['node_power'][j,0]
                if(pos>-0.5):
                    A[pos,pos]=1
                pos=self.varDictImg[i].me['node_power'][j,0]
                if(pos>-0.5):
                    A[pos,pos]=1
                    
                
            
            # update the coefficients in front of v, ell,
            for j in range(3):
                for k in range(j+1):
                    
                    isdiag=1 if (j==k) else 2
                    
                    pos=self.varDictReal[i].parent['node_vol'][j,k]
                    if(pos>-0.5):
                        A[pos,pos]=isdiag
                    pos=self.varDictImg[i].parent['node_vol'][j,k]
                    if(pos>-0.5):
                        A[pos,pos]=isdiag
                        
                    pos=self.varDictReal[i].me['node_vol'][j,k]
                    if(pos>-0.5):
                        A[pos,pos]=isdiag*2
                    pos=self.varDictImg[i].me['node_vol'][j,k]
                    if(pos>-0.5):
                        A[pos,pos]=isdiag*2
                        
                    pos=self.varDictReal[i].me['branch_cur'][j,k]
                    if(pos>-0.5):
                        A[pos,pos]=isdiag*(num_children+1)
                    pos=self.varDictImg[i].me['branch_cur'][j,k]
                    if(pos>-0.5):
                        A[pos,pos]=isdiag*(num_children+1)
                        
                    for childID in n.neighbor.children:
                        pos=self.varDictReal[i].children[childID]['branch_cur'][j,k]
                        if(pos>-0.5):
                            A[pos,pos]=isdiag
                        pos=self.varDictImg[i].children[childID]['branch_cur'][j,k]
                        if(pos>-0.5):
                            A[pos,pos]=isdiag   
            self.A[i]=A
    
    def __hash_pos(self,node_set):
        """
        Hash the real/imaginary part of each variable in a one dimensional
        real arrays
        """

        for i in range(self.numBus):
            n=node_set[i]
            self.varDictReal.append(VarDict(n))
            self.varDictImg.append(VarDict(n))
            tmp_numVar=0
            
            
            #update the branch power variables
            #S=VI^* is not a hermittian matrix and needs to check all 9 entries
            for j in range(3):
                for k in range(3):
                    if(n.seq_num!=0 and n.phase[j]==1 and n.phase[k]==1):
                        self.varDictReal[i].me['branch_power'][j,k]=tmp_numVar
                        tmp_numVar+=1
                        self.varDictImg[i].me['branch_power'][j,k]=tmp_numVar
                        tmp_numVar+=1
                    for childID in n.neighbor.children: 
                        c=node_set[childID]
                        if(c.phase[j]==1 and c.phase[k]==1):
                            self.varDictReal[i].children[childID]['branch_power'][j,k]=tmp_numVar
                            tmp_numVar+=1
                            self.varDictImg[i].children[childID]['branch_power'][j,k]=tmp_numVar
                            tmp_numVar+=1   
            
 
            

            #update power injection s, which is a 3x1 vector
            for j in range(3):
                if(n.phase[j]==1):
                    if(n.control_info[j]==1 or n.control_info[j]==3):
                        self.varDictReal[i].me['node_power'][j,0]=tmp_numVar
                        tmp_numVar+=1
                    if(n.control_info[j]==2 or n.control_info[j]==3):
                        self.varDictImg[i].me['node_power'][j,0]=tmp_numVar
                        tmp_numVar+=1
                        


            
            #Update v, ell, which is hermitian and lower diagonal is enough
            #root does not have parents
            if n.seq_num!=0:
                p=node_set[n.neighbor.parent]
                for j in range(3):
                    for k in range(j+1):
                        if(p.phase[j]==1 and p.phase[k]==1):
                            self.varDictReal[i].parent['node_vol'][j,k]=tmp_numVar
                            tmp_numVar+=1
                            #Diagonal element is real
                            if k!=j :
                                self.varDictImg[i].parent['node_vol'][j,k]=tmp_numVar
                                tmp_numVar+=1 
                                
                        if(n.phase[j]==1 and n.phase[k]==1):
                            self.varDictReal[i].me['node_vol'][j,k]=tmp_numVar
                            tmp_numVar+=1
                            self.varDictReal[i].me['branch_cur'][j,k]=tmp_numVar
                            tmp_numVar+=1
                            #Diagonal element is real
                            if(k!=j):
                                self.varDictImg[i].me['node_vol'][j,k]=tmp_numVar
                                tmp_numVar+=1
                                self.varDictImg[i].me['branch_cur'][j,k]=tmp_numVar
                                tmp_numVar+=1
                            

            for childID in n.neighbor.children:
                c=node_set[childID]
                for j in range(3):
                    for k in range(j+1):
                        if(c.phase[j]==1 and c.phase[k]==1):
                            self.varDictReal[i].children[childID]['branch_cur'][j,k]=tmp_numVar
                            tmp_numVar+=1
                            if(k!=j):
                                self.varDictImg[i].children[childID]['branch_cur'][j,k]=tmp_numVar
                                tmp_numVar+=1

            self.numVar.append(tmp_numVar)
            
            
            
            
        
