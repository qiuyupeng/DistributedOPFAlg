import xlrd
import os
import numbers

from Exception import *
from Node import *



class Cap:
    """
    Put all the load infos in the Node object. 
    
    Member variables:
    
    self.rawData
    """
    
    
    
    def __init__(self,FileName,node_set,node_dict):
        
        self.__readDataFromFile(FileName)
        self.__updateNode(node_set,node_dict)
        
        
        
        
    def __updateNode(self,node_set,node_dict):
        """
        Update fixed_load,phase,control info in the node
        
        """
        try:
            for x in self.rawData:
                if(x[0] not in node_dict.keys()):
                    raise MyException('Bus '+str(x[0])+' does not exist')
                node_id=node_dict[x[0]]
                node_set[node_id].upper_control=node_set[node_id].upper_control+np.array([complex(0,x[1]),complex(0,x[2]),complex(0,x[3])])/1000 #convert it to MW.
                if(x[1]!=0):
                    node_set[node_id].phase[0]=1
                    if(node_set[node_id].control_info[0]==0 or node_set[node_id].control_info[0]==1):
                        node_set[node_id].control_info[0]+=2     
                if(x[2]!=0):
                    node_set[node_id].phase[1]=1
                    if(node_set[node_id].control_info[1]==0 or node_set[node_id].control_info[1]==1):
                        node_set[node_id].control_info[1]+=2
                if(x[3]!=0):
                    node_set[node_id].phase[2]=1
                    if(node_set[node_id].control_info[2]==0 or node_set[node_id].control_info[2]==1):
                        node_set[node_id].control_info[2]+=2
        except MyException as e:
            print e
        
        
        
        
        
        
        
        
        
        
        
        
        
    def __readDataFromFile(self,File):
        
        try:
            if(os.path.exists(File)==False):
                raise MyException(File)
        except MyException as e:
            print "File", e, "does not exists!"
        raw_data=xlrd.open_workbook(File).sheets()[0] #only sheet 1 has data
        
        nrows=raw_data.nrows
        ncols=raw_data.ncols
        
        data=[[None]*ncols for i in range(nrows-5)]
        
        for row in range(nrows-5):
            for col in range(ncols):
                if(raw_data.cell_type(row+4,col)!=0):
                    data[row][col]=int(raw_data.cell(row+4,col).value)
                else:
                    data[row][col]=0
        self.rawData=data
        #print self.rawData

        