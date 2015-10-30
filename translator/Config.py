import numpy as np


class Config:
    """
    Define the constants for each configuration.
    Currently only the impedence per mile for different configuraiton
    """
    
    def __init__(self):
        self.feeder_dict=dict([])
        self.feeder_dict['feeder13']=self.IEEE_feeder13
        self.feeder_dict['feeder2']=self.IEEE_feeder2
        self.feeder_dict['feeder3']=self.IEEE_feeder2
        self.feeder_dict['feeder4']=self.IEEE_feeder2
        
        self.feeder_dict['feeder34']=self.IEEE_feeder34
        self.feeder_dict['feeder123']=self.IEEE_feeder123
        
        self.feeder_sub=dict([])
        self.feeder_sub['feeder13']=650
        self.feeder_sub['feeder2']=650
        self.feeder_sub['feeder3']=650
        self.feeder_sub['feeder4']=650
        self.feeder_sub['feeder34']=800
        self.feeder_sub['feeder123']=150
    
    def IEEE_feeder2(self):
        """
        return a np.array object of the impedence
        """
        z_config=dict([])
        
        z_config[601]=np.matrix([[0.3465+1.0179j,0.1560+0.5017j,0.1580+0.4236j],\
        [0.1560+0.5017j,0.3375+1.0478j,0.1535+0.3849j],[0.1580+0.4236j,0.1535+0.3849j,0.3414+1.0348j]])*100
        
        z_config[603]=np.matrix([[0.3465+1.0179j,0.1560+0.5017j,0.1580+0.4236j],\
        [0.1560+0.5017j,0.3375+1.0478j,0.1535+0.3849j],[0.1580+0.4236j,0.1535+0.3849j,0.3414+1.0348j]])*100

        return z_config
    
    
    def IEEE_feeder13(self):
        """
        return a np.array object of the impedence
        """
        z_config=dict([])
        
        z_config[601]=np.matrix([[0.3465+1.0179j,0.1560+0.5017j,0.1580+0.4236j],\
        [0.1560+0.5017j,0.3375+1.0478j,0.1535+0.3849j],[0.1580+0.4236j,0.1535+0.3849j,0.3414+1.0348j]])
        
        z_config[602]=np.matrix([[0.7526+1.1814j,0.1580+0.4236j,0.1560+0.5017j], \
        [0.1580+0.4236j,0.7475+1.1983j,0.1535+0.3849j],[0.1560+0.5017j,0.1535+0.3849j,0.7436+1.2112j]])
        
        z_config[603]=np.matrix([[0,0,0],[0,1.3294+1.3471j,0.2066+0.4591j],[0,0.2066+0.4591j,1.3238+1.3569j]])
        
        z_config[604]=np.matrix([[1.3238+1.3569j,0,0.2066+0.4591j],[0,0,0],[0.2066+0.4591j,0,1.3294+1.3471j]])
        
        z_config[605]=np.matrix([[0,0,0],[0,0,0],[0,0,1.3292+1.3475j]])
        
        z_config[606]=np.matrix([[0.7982+0.4463j,0.3192+0.0328j,0.2849-0.0143j],\
        [0.3192+0.0328j,0.7891+0.4041j,0.3192+0.0328j],[0.2849-0.0143j,0.3192+0.0328j,0.7982+0.4463j]])
        
        z_config[607]=np.matrix([[1.3425+0.5124j,0,0],[0,0,0],[0,0,0]])
        
        z_config[u'XFM-1']=np.matrix([[0,0,0],[0,0,0],[0,0,0]])
        
        z_config[u'Switch']=np.matrix([[0,0,0],[0,0,0],[0,0,0]])
        
        return z_config
    
    def IEEE_feeder34(self):
        
        z_config=dict([])
        
        z_config[300]=0.1*np.matrix([[1.3368+1.3343j,0.2101+0.5779j,0.2130+0.5015j],
        [0.2101+0.5779j,1.3238+1.3569j,0.2066+0.4591j],[0.2130+0.5015j,0.2066+0.4591j,1.3294+1.3471j]])
                
        z_config[301]=0.1*np.matrix([[1.9300+1.4115j,0.2327+0.6442j,0.2359+0.5691j],\
        [0.2327+0.6442j,1.9157+1.4281j,0.2288+0.5238j],[0.2359+0.5691j,0.2288+0.5238j,1.9219+1.4209j]])
    
        z_config[302]=0.1*np.matrix([[2.7995+1.4855j,0,0],[0,0,0],[0,0,0]])
        
        z_config[303]=0.1*np.matrix([[0,0,0],[0,2.7995+1.4855j,0],[0,0,0]])
        
        z_config[304]=0.1*np.matrix([[0,0,0],[0,1.9217+1.4212j,0],[0,0,0]])


        
        z_config[u'XFM-1']=np.matrix([[0,0,0],[0,0,0],[0,0,0]])
        
        
        return z_config
    
    def IEEE_feeder123(self):
        
        z_config=dict([])
        
        z_config[1]=np.matrix([[0.4576+1.0780j,0.1560+0.5017j,0.1535+0.3849j],\
        [0.1560+0.5017j,0.4666+1.0482j,0.1580+0.4236j],[0.1535+0.3849j,0.1580+0.4236j,0.4615+1.0651j]])
        
        z_config[2]=np.matrix([[0.4666+1.0482j,0.1580+0.4236j,0.1560+0.5017j],\
        [0.1580+0.4236j,0.4615+1.0651j,0.1535+0.3849j],[0.1560+0.5017j,0.1535+0.3849j,0.4576+1.0780j]])
        
        z_config[3]=np.matrix([[0.4615+1.0651j,0.1535+0.3849j,0.1580+0.4236j],\
        [0.1580+0.4236j,0.4576+1.0780j,0.1560+0.5017j],[0.1580+0.4236j,0.1560+0.5017j,0.4666+1.0482j]])
        
        z_config[4]=np.matrix([[0.4615+1.0651j,0.1580+0.4236j,0.1535+0.3849j],\
        [0.1580+0.4236j,0.4666+1.0482j,0.1560+0.5017j],[0.1535+0.3849j,0.1560+0.5017j,0.4576+1.0780j]])
        
        z_config[5]=np.matrix([[0.4666+1.0482j,0.1560+0.5017j,0.1580+0.4236j],\
        [0.1560+0.5017j,0.4576+1.0780j,0.1535+0.3849j],[0.1580+0.4236j,0.1535+0.3849j,0.4615+1.0651j]])
        
        z_config[6]=np.matrix([[0.4576+1.0780j,0.1535+0.3849j,0.1560+0.5017j],\
        [0.1535+0.3849j,0.4615+1.0651j,0.1580+0.4236j],[0.1560+0.5017j,0.1580+0.4236j,0.466+1.0482j]])
        
        z_config[7]=np.matrix([[0.4576+1.0780j,0,0.1535+0.3849j],\
        [0,0,0],[0.1535+0.3849j,0,0.4615+1.0651j]])
        
        z_config[8]=np.matrix([[0.4576+1.0780j,0.1535+0.3849j,0],[0.1535+0.3849j,0.4615+1.0651j,0],[0,0,0]])
        
        z_config[9]=np.matrix([[1.3292+1.3475j,0,0],[0,0,0],[0,0,0]])
        
        z_config[10]=np.matrix([[0,0,0],[0,1.3292+1.3475j,0],[0,0,0]])
        
        z_config[11]=np.matrix([[0,0,0],[0,0,0],[0,0,1.3292+1.3475j]])
        
        z_config[12]=np.matrix([[1.5209+0.7521j,0.5198+0.2775j,0.4924+0.2157j],
        [0.5198+0.2775j,1.5329+0.7162j,0.5198+0.2775j],[0.4924+0.2157j,0.5198+0.2775j,1.5209+0.7521j]])
        
        z_config[u'switch']=np.matrix([[0,0,0],[0,0,0],[0,0,0]])
        
        return z_config
        




if __name__=="__main__":
    x=Config()
    print x.feeder_dict['feeder123']()

        
        