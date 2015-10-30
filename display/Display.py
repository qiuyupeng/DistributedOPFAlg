import matplotlib.pyplot as plt
import os
import math

import numpy as np

class Display:
    """
    Display the useful results regarding all the raw data:
    
    1.Primal residual v.s. Dual residual
    2.Mean executing time.
    3.Objective value (loss)
    """
    
    def __init__(self,pRes,dRes,exeTime,var,obj,feederName):
        
        self.pRes=np.array(pRes)
        self.dRes=np.array(dRes)
        self.exeTime=np.array(exeTime)
        self.node_set=var
        self.obj=obj
        self.dir=os.getcwd()+'/result/'+feederName
        self.feederName=feederName
        
        if(not os.path.exists(self.dir)):
            os.makedirs(self.dir)
        
        
        
    def plot_primalDualRes(self):
        """
        Plot the primal and dual residual 
        """
        plt.semilogy(self.pRes,linewidth=2.0)
        plt.semilogy(self.dRes,linewidth=2.0,ls='--')
        plt.legend(['Primal Residual','Dual Residual'])
        plt.xlabel('Number of Iterations')
        plt.ylabel('Residual error')
        plt.title('Primal and Dual Residual')
        
        plt.savefig(self.dir+'/residualerror.png')
        plt.clf()
        
    def plot_normal_primalDualRes(self):
        """
        Plot the primal and dual residual (divide by number of agents )
        """
        numBus=math.sqrt(len(self.node_set))
        plt.semilogy(self.pRes/numBus,linewidth=2.0)
        plt.semilogy(self.dRes/numBus,linewidth=2.0,ls='--')
        plt.legend(['Primal Residual','Dual Residual'])
        plt.xlabel('Number of Iterations')
        plt.ylabel('Residual error')
        plt.title('Primal and Dual Residual')
        
        plt.savefig(self.dir+'/residualerror.png')
        plt.clf()
    
    def show_simulationRst(self):
        
        for n in self.node_set:
            print n
    
    def write_log(self):
        
        with open(self.dir+'/log.txt','w') as f:
            f.write('The network thermal loss is '+str(self.thermalLoss())+'MW\n\n')
            x=self.statTime()
            f.write('The total computation time is '+str(x[0])[0:5]+'s\n')
            f.write('The mean computation time is '+str(x[1])[0:5]+'s\n')
            f.write('The std of computation time is '+str(x[2])[0:5]+'s\n\n')
            
            
    
    def statTime(self):
        
        totalTime=np.sum(self.exeTime)
        avgTime=np.mean(self.exeTime)
        stdTime=np.std(self.exeTime)
        
        return (totalTime,avgTime,stdTime)
    
    def thermalLoss(self):
        
        loss=0
        for n in self.node_set:
            loss+=n.z_npower.sum().real
        
        return loss
        
        
        
        
