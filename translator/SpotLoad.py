import xlrd
import os
import numbers

from Exception import *
from Node import *



class Spot_Load:
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
        Update fixed_load and phase info in the node
        """
        try:
            for x in self.rawData:
                if(x[0] not in node_dict.keys()):
                    raise MyException('Bus '+str(x[0])+' does not exist')
                node_id=node_dict[x[0]]
                node_set[node_id].fixed_load=node_set[node_id].fixed_load-np.array([complex(x[2],x[3]),complex(x[4],x[5]),complex(x[6],x[7])])/1000 #convert it to MW.
                if(x[2]!=0 or x[3]!=0):
                    node_set[node_id].phase[0]=1
                if(x[4]!=0 or x[5]!=0):
                    node_set[node_id].phase[1]=1
                if(x[6]!=0 or x[7]!=0):
                    node_set[node_id].phase[2]=1
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
        
        try:
            for row in range(nrows-5):
                for col in range(ncols):
                    if(raw_data.cell_type(row+4,col)==0):
                        raise MyException("row: "+str(row)+" col:"+str(col))
                    if(isinstance(raw_data.cell(row+4,col).value,numbers.Number)):
                        data[row][col]=int(raw_data.cell(row+4,col).value)
                    else:
                        data[row][col]=raw_data.cell(row+4,col).value
        except MyException as e:
            print "Data in", e, "does not exists!"
            
        self.rawData=data
        #print self.rawData

        