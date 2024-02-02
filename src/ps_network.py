import numpy as np 
import pandas as pd
import io


class PSNetwork:
    '''
    this class loads a data file of microvascular network connectivity from
    Pries et al. 1990 and provides methods to access the data'''

    def __init__(self, seg_data: pd.DataFrame, boundary_nodes: pd.DataFrame, n_vessels):
        self.filename = filename
        self.seg_data = seg_data
        self.boundary_nodes = boundary_nodes
        self.n_vessels = n_vessels

    ### I/O METHODDS ###
    @classmethod
    def from_file(cls, filename: str):
        '''
        read the data in from a file and return a PSNetwork object'''
        # format the data into two dataframes since the initial conformation is... inconvenient

        seg_data = []
        boundary_nodes = []
        reading_seg_data = True

        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if reading_seg_data:
                    if line.startswith('RAT MESENTERY'):
                        # we are reading the title string containing the number of segments
                        # title_str = f.readline()
                        title_str = line
                        title_str.split(' ')
                        n_vessels = title_str[2]
                    elif line.startswith('boundary_nodes'):
                        reading_seg_data = False
                        boundary_nodes.append(line)
                        pass
                    elif line == '':
                        pass
                    else:
                        seg_data.append(line)
                else:
                # not reading segment data, now reading boundary node data
                    boundary_nodes.append(line)

        seg_data = [item.replace('\t', ' ') for item in seg_data]

        # read the data into pandas dataframes
        seg_data = pd.read_csv(io.StringIO('\n'.join(seg_data)), sep='\s+')
        boundary_nodes = pd.read_csv(io.StringIO('\n'.join(boundary_nodes)), sep='\s+')

        return cls(seg_data, boundary_nodes, n_vessels)
    
    def to_file(self, filename: str):
        '''
        write the data to a file, placeholder method at the moment'''
        pass

    ### END OF I/O METHODS ###


    def generate_vessel_map(self):
        '''
        generate a map of the vessels in the network'''
        pass





if __name__ == "__main__":
    filename = 'data/913-segments.txt'
    network = PSNetwork.from_file(filename)