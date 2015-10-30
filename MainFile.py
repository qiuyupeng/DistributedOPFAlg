import os
import sys
sys.path.insert(0,os.getcwd()+'/translator')
sys.path.insert(0,os.getcwd()+'/algorithm')
sys.path.insert(0,os.getcwd()+'/display')
import time

from Line import *
from SpotLoad import *
from Cap import *
import Algorithm as Alg
import Display as Disp


class MainFile:
    """
    define some basics about the position of the existing IEEE files.
    """
    def __init__(self,feedername='feeder2'):
        """
        string of the file names
        """
        self.directory=os.getcwd()
        self.feedername=feedername
        self.cap_data='Cap Data.xls'
        self.distributed_load='Distriuted Load Data.xls'
        self.line_data='Line Data.xls'
        self.spot_load='Spot Load Data.xls'
        
        
        tmp=Line(self.directory+'/'+self.feedername+'/'+self.line_data,feedername)  #the node set is inferred from Line Data.xls
        self.node_set=tmp.node_set
        self.node_dict=tmp.node_dict
        
        tmp=Spot_Load(self.directory+'/'+self.feedername+'/'+self.spot_load,self.node_set,self.node_dict)
        tmp=Cap(self.directory+'/'+self.feedername+'/'+self.cap_data,self.node_set,self.node_dict)
        
        #for n in self.node_set:
        #    print n
        #self.node_set[1].lower_control=np.array([-10j,-10j,-10j],dtype=complex)

        algorithm=Alg.Algorithm(self.node_set)
        pRes,dRes,exeTime,var,obj=algorithm.getResult()

        print "************Simulation Result*************"
        self.disp=Disp.Display(pRes,dRes,exeTime,var,obj,self.feedername)
        self.disp.plot_primalDualRes()
        self.disp.plot_normal_primalDualRes()
        self.disp.write_log()
        #self.disp.show_simulationRst()

        
        
        
        

        
        
        
        
        




if __name__=="__main__":
    x=MainFile()
    #print x.node_set[1].impedance
    #for y in x.node_set:
    #    print y


    
    
   
    

    
    
    
    