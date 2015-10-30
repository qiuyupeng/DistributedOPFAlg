import cmath
import math
import numpy as np
import numpy.linalg as LA

class ForwardBackSweep:
    """
    The backward/forward sweep algorithm to solve unbalanced 3-phase radial network
    """
    
    def __init__(self):
        pass
        
    
    def solver(self,node_set):
        
        
        self.__initilization(node_set)
        
        
        num_iter=0
        error=1000
        while (error>1e-10 and num_iter<100):
            self.__Forward(node_set,0)
            self.__Backward(node_set,0)
            error=self.__error(node_set)
            num_iter+=1
            print "Iteration "+str(num_iter)+": error residual: "+str(error)
        
        self.showResult(node_set)
        
        print "Forward/Backforward Algorithm Solution Done!"
        print "*********************************"
        
        
        
    def showResult(self,node_set):
        
        print "Results:"
        for n in node_set:
            s='Label: '+str(n.label)+'\n'+\
            'Sequence number: '+str(n.seq_num)+'\n'+\
            'Node voltage: a: '+str(cmath.polar(n.node_vol[0,0]))+'b: '+str(cmath.polar(n.node_vol[1,0]))+'c: '+str(cmath.polar(n.node_vol[2,0]))+'\n'+\
            'Branch current: '+str(n.branch_cur)+'\n'\
            'Load: '+str(n.fixed_load)+'\n\n'
            print s
        

        
        
        
    def __initilization(self,node_set):
        """
        Intialize the self.node_vol variables
        """
        
        print "*********************************"
        
        for x in node_set:
            x.node_vol=np.transpose(np.matrix([cmath.exp(0), cmath.exp(complex(0,math.pi*2/3)), cmath.exp(complex(0,-math.pi*2/3))]))
        
        print "Forward/Backward Algorithm Initialization Done!"
        
        
        
        
    def __error(self,node_set):
        """
        Check the infeasibility of the solution
        """
        error=0
        for n in node_set:
            if(n.seq_num!=0):
                error+=LA.norm(n.node_vol-node_set[n.neighbor.parent].node_vol-n.impedance*n.branch_cur)
                #print n.node_vol, '\n', node_set[n.neighbor.parent].node_vol
                
        return error
                
       
    def __Forward(self,node_set,cur_nodeID):
        """
        Forward Sweep
        """
        
        n=node_set[cur_nodeID]
        if(cur_nodeID!=0):
            n.node_vol=node_set[n.neighbor.parent].node_vol+n.impedance*n.branch_cur
        for i in n.neighbor.children:
            self.__Forward(node_set,i)
     
    def __Backward(self,node_set,cur_nodeID):
        """
        Backward sweep
        """
        
        n=node_set[cur_nodeID]
        for i in range(3):
            n.node_cur[i,0]=(n.fixed_load[i]/n.node_vol[i,0]).conjugate()
            #n.node_cur[i,0]=n.fixed_load[i]/n.node_vol[i,0]
            n.branch_cur[i,0]=n.node_cur[i,0]
        
        for i in n.neighbor.children:
            self.__Backward(node_set,i)
            n.branch_cur=n.branch_cur+node_set[i].branch_cur
            


            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        
        
        
