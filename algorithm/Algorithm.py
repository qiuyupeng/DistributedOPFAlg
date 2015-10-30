import Initialization
from Solver import Solver
import Constant as ct

class Algorithm:
    """
    Interface from mainfile
    
    Input: node_set
    """
    
    def __init__(self,node_set):
        self.node_set=node_set
        
        #Initialize all the variables
        Initialization.Initialization(self.node_set)
        
        self.solver=Solver(self.node_set)
    
    def getResult(self):
        
        primalRes,dualRes,exeTime=self.solver.start(ct.Constant().numIter)
        var=self.node_set
        
        obj=0
        for n in var:
            obj+=sum(n.z_npower)[0,0].real
        return (primalRes,dualRes,exeTime,var,obj)
        
        