import xlrd
import os
import numbers
from numpy import linalg as LA

from Config import *
from Exception import *
from Node import *

class Line:
    """
    Extract data from Line Data.xls and construct the graph 
    represented adjancey matrix
    
    Member variables:
    
    self.node_set
    self.node_dict
    
    
    self.rawData
    self.Substation
    
    """
    
    def __init__(self,FileName,FeederName):
        
        #Consant
        self.IMPEDANCE_PER_UNIT=1.0e6/(7200**2)
        
        self.FeederName=FeederName
        
        self.Substation=Config().feeder_sub[self.FeederName]
        self.node_dict=dict([])
        
        self.node_dict[self.Substation]=0 #substation bus as 0
        
        
        self.__readDataFromFile(FileName)
        self.__HandleRawData()
        self.__updateImpedance()
        

        
    def __updateImpedance(self):
        
        config_file=Config().feeder_dict[self.FeederName]()
        
        fr=[x[0] for x in self.rawData]
        to=[x[1] for x in self.rawData]
        z_length=[x[2]/5280. for x in self.rawData]
        z_config=[x[3] for x in self.rawData]
        
        try:
            for i in range(len(z_length)):
                if self.node_set[self.node_dict[fr[i]]].neighbor.parent==self.node_dict[to[i]] :
                    if(self.node_set[self.node_dict[fr[i]]].impedance!=None):
                        raise MyException('Impedance at Node'+str(self.node_dict[fr[i]])+' already set!')
                    self.node_set[self.node_dict[fr[i]]].impedance=config_file[z_config[i]]*z_length[i]*self.IMPEDANCE_PER_UNIT
                    self.__definePhase(config_file[z_config[i]],self.node_set[self.node_dict[fr[i]]])
                        
                elif self.node_set[self.node_dict[to[i]]].neighbor.parent==self.node_dict[fr[i]] :
                    if(self.node_set[self.node_dict[to[i]]].impedance!=None):
                        raise MyException('Impedance at Node'+str(self.node_dict[to[i]])+' already set!')
                    self.node_set[self.node_dict[to[i]]].impedance=config_file[z_config[i]]*z_length[i]*self.IMPEDANCE_PER_UNIT
                    self.__definePhase(config_file[z_config[i]],self.node_set[self.node_dict[to[i]]])
        except MyException as e:
            print e             
        
        
    def __definePhase(self,z_config,node):
        """
        Determine the phase of a node
        """
        #switch and xfm are assumed to be three phase
        if(LA.norm(z_config)<1e-5):
            node.phase=np.ones(3,dtype=np.int)
        else:
            for i in range(3):
                if(z_config[i,i]<1e-10):
                    node.phase[i]=0
                else:
                    node.phase[i]=1
        
            
    
    def __HandleRawData(self):
        
        seq_num=1  #name the current node 
        
        self.node_set=[Node(0)]
        self.node_set[0].label=self.Substation
        self.node_set[0].phase=np.array([1,1,1])
        self.node_set[0].control_info=np.array([3,3,3])
        self.node_set[0].upper_control=1e4*np.array([1+1j,1+1j,1+1j],dtype=complex)
        
        for x in self.rawData:
            if(x[0] not in self.node_dict.keys()):
                self.node_dict[x[0]]=seq_num
                self.node_set.append(Node(seq_num))
                self.node_set[seq_num].label=x[0]
                seq_num=seq_num+1
            if(x[1] not in self.node_dict.keys()):
                self.node_dict[x[1]]=seq_num
                self.node_set.append(Node(seq_num))
                self.node_set[seq_num].label=x[1]
                seq_num=seq_num+1
                
            fr=self.node_dict[x[0]]
            to=self.node_dict[x[1]]
            self.node_set[fr].neighbor.neighbor.add(to)
            self.node_set[to].neighbor.neighbor.add(fr)
            
        self.__BFS()  #populate the parent and children entries in neighbor_set 
    
    
    def __BFS(self):
        
        node_queue=[0]
        node_set=set([])
        while len(node_queue)>0 :
            n=node_queue.pop(0)
            node_set.add(n)
            for i in self.node_set[n].neighbor.neighbor:
                if i not in node_set :
                    node_queue.append(i)
                    self.node_set[i].neighbor.parent=n
                    self.node_set[n].neighbor.children.add(i)
            
                
        
            
            
                
            
        
    def __readDataFromFile(self,File):
        
        try:
            if(os.path.exists(File)==False):
                raise MyException(File)
        except MyException as e:
            print "File", e, "does not exists!"
        raw_data=xlrd.open_workbook(File).sheets()[0] #only sheet 1 has data
        
        nrows=raw_data.nrows
        ncols=raw_data.ncols
        
        data=[[None]*ncols for i in range(nrows-3)]
        
        try:
            for row in range(nrows-3):
                for col in range(ncols):
                    if(raw_data.cell_type(row+3,col)==0):
                        raise MyException("row: "+str(row)+" col:"+str(col))
                    if(isinstance(raw_data.cell(row+3,col).value,numbers.Number)):
                        data[row][col]=int(raw_data.cell(row+3,col).value)
                    else:
                        data[row][col]=raw_data.cell(row+3,col).value
        except MyException as e:
            print "Data in", e, "does not exists!"
            
        self.rawData=data




        



if __name__=="__main__":
    x=Line()
    x.readDataFromFile('cao.xls')




        