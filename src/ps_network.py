import numpy as np 
import pandas as pd
import os
import io
from svzerodtrees._config_handler import ConfigHandler, Vessel, BoundaryCondition, Junction, SimParams
from svzerodtrees._result_handler import ResultHandler
from svzerodtrees.utils import *


class PSNetwork:
    '''
    this class loads a data file of microvascular network connectivity from
    Pries et al. 1990 and provides methods to access the data'''

    def __init__(self, seg_data: pd.DataFrame, boundary_nodes: list, n_vessels: int):

        self.seg_data = seg_data
        self.boundary_nodes = boundary_nodes
        self.n_vessels = n_vessels

        if 'length' not in seg_data.columns:
            self.add_lengths()

        self.config_handler = self.generate_zerod_model()

        self.result_handler = ResultHandler.from_config_handler(self.config_handler)


        # print(self.config_handler.config)
    ### I/O METHODDS ###
    @classmethod
    def from_file(cls, filename: str):
        '''
        read the data in from a file and return a PSNetwork object'''
        # format the data into two dataframes since the initial conformation is... inconvenient

        seg_data_list = []
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
                        title_str = title_str.split(' ')
                        n_vessels = int(title_str[2])
                    elif line.startswith('boundary_nodes'):
                        reading_seg_data = False
                        pass
                    elif line == '':
                        pass
                    else:
                        seg_data_list.append(line)
                else:
                # not reading segment data, now reading boundary node data
                    if line == '':
                        pass
                    else:
                        boundary_nodes.append(int(line))

        seg_data_list = [item.replace('\t', ' ') for item in seg_data_list]

        # read the data into pandas dataframes
        seg_data = pd.read_csv(io.StringIO('\n'.join(seg_data_list)), sep='\s+')
        # change dtype of columns
        seg_data = seg_data.astype({'segment_name': 'str', 'node_from': 'int', 'node_to': 'int', 'segment_diameter': 'float'})

        return cls(seg_data, boundary_nodes, n_vessels)
    
    def to_json(self, filename: str):
        '''
        write the data to a file, placeholder method at the moment'''
        
        pass

    def add_lengths(self):
        '''
        add the lengths of the segments to the segment data based on the relationship l = 12.4 d^1.1'''

        lengths = []
        for diameter in self.seg_data.segment_diameter:
            lengths.append(12.4 * diameter**1.1)
        
        self.seg_data['length'] = lengths

    ### END OF I/O METHODS ###


    def generate_zerod_model(self):
        '''
        generate a map of the vessels in the network as zerodsolver objects'''
        
        # junction map is not trivial to generate, we will initialize that first
        junction_map = {}
        n_junctions = self.seg_data.node_from.max()
        for i in range(n_junctions):
            # unforunately the index in secomb's data is 1-indexed grrrr
            junction_map[i + 1] = Junction.from_config({'inlet_vessels': [], 
                                                    'outlet_vessels': [],
                                                    'junction_name': 'J' + str(i),
                                                    'junction_type': 'NORMAL_JUNCTION'})

        vessel_map = {}
        bc_map = {}
        for row in self.seg_data.iterrows():
            vessel_map[row[0]] = Vessel.from_config(self.create_vessel_config(row))
            # add the vessel to the junction map
            if row[1]['node_from'] in self.boundary_nodes:
                # add inflow bc
                bc_name = 'inlet' + str(row[1]['node_from'])
                bc_map[bc_name] = BoundaryCondition.from_config({'bc_name': bc_name,
                                                                 'bc_type': 'FLOW',
                                                                 'bc_values': {
                                                                     'Q': [100.0, 100.0], # placeholder
                                                                     't': [0.0, 1.0]
                                                                 }
                                                             })
                vessel_map[row[0]].bc = {"inlet": bc_name}
            else:
                # add to the outelet vessels of a junction
                junction_map[row[1]['node_from']].outlet_branches.append(row[0])

            if row[1]['node_to'] in self.boundary_nodes:
                # add outflow bc
                bc_name = 'outlet' + str(row[1]['node_to'])
                bc_map[bc_name] = BoundaryCondition.from_config({'bc_name': bc_name,
                                                                 'bc_type': 'RESISTANCE',
                                                                 'bc_values': {
                                                                     'R': 10.0, # placeholder
                                                                     'Pd': 0.0
                                                                 }
                                                             })
                vessel_map[row[0]].bc = {"outlet": bc_name}
            else:
                # add to the inlet vessels of a junction
                junction_map[row[1]['node_to']].inlet_branches.append(row[0])
        
        # clean the boundary nodes from the junction map and add bcs
        junction_map = {node_id: junction for node_id, junction in junction_map.items() if junction.outlet_branches != [] and junction.inlet_branches != []}

        # create the sim params
        sim_params = SimParams.from_config({'density': 1.06,
                                            'model_name': str(self.n_vessels) + '-vessels_ps_network',
                                            'number_of_cardiac_cycles': 1,
                                            'number_of_time_pts_per_cardiac_cycle': 10,
                                            'viscosity': 0.04,
                                            })
        
        # create the config handler
        config_handler = ConfigHandler({'junctions': [junction.to_dict() for junction in junction_map.values()],
                                           'vessels': [vessel.to_dict() for vessel in vessel_map.values()],
                                           'boundary_conditions': [bc.to_dict() for bc in bc_map.values()],
                                           'simulation_parameters': sim_params.to_dict()
                                           },
                                           is_pulmonary=False)

        return config_handler


    def create_vessel_config(self, seg_data_row, viscosity=0.04, stenosis_coefficient=0.0):
        '''
        create a 0d vessel config dict from the seg_data_row
        
        :param id: the id of the vessel, from the dataframe row id
        :param name: the segemnt_name parameter
        :param diameter: the diameter of the vessel in mm
        :param length: the length of the vessel in mm
        :param viscosity: the viscosity of the blood in the vessel
        :param stenosis_coefficient: the stenosis coefficient of the vessel'''
        
        # convert from mm to cm
        length = seg_data_row[1]['length'] / 10
        diameter = seg_data_row[1]['segment_diameter'] / 10

        # zerod element values
        R = 8 * viscosity * length / (np.pi * diameter**4)
        C = 0.0 # can implement this later
        L = 0.0 # can implement this later

        # name
        name = 'branch' + str(seg_data_row[0]) + '_seg0'

        return {
                'vessel_id': seg_data_row[0],
                'vessel_length': length,
                'vessel_name': name,
                'zero_d_element_type': "BloodVessel",
                'zero_d_element_values': {
                    'R_poiseuille': R,
                    'C': C,
                    'L': L,
                    'stenosis_coefficient': stenosis_coefficient
                },
            }


    def run_simulation(self):
        '''
        run a 0d simulation of the network
        '''

        print('running simulation...')
        self.config_handler.simulate(self.result_handler, 'preop')

        print('ran simulation!')

    
    def get_network_pressures(self, steady: bool=True):
        '''
        get lists of the inlet, outlet and mean pressure of the vessels in the network
        '''
        
        p_in = []
        p_out = []
        for branch in self.config_handler.vessel_map.keys():
            p_in.append(get_branch_result(self.result_handler.results['preop'], 'pressure_in', branch, steady=steady)) 
            p_out.append(get_branch_result(self.result_handler.results['preop'], 'pressure_out', branch, steady=steady))

        p_mean = [(p + p_out[i]) / 2 for i, p in enumerate(p_in)]

        return p_in, p_out, p_mean
    
    
    def get_network_flows(self, steady: bool=True):
        '''
        get a list of the flow for each vessel in the network
        '''
        
        q = []
        for branch in self.config_handler.vessels_map.keys():
            q.append(get_branch_result(self.result_handler.results['preop'], 'flow_in',branch, steady=steady)) 

        return q


    def get_network_wss(self, steady: bool=True):
        '''
        get a list of the wall shear stress for each vessel in the network
        '''
        
        wss = []
        for branch in self.config_handler.vessel_map.keys():
            wss.append(get_wss(self.result_handler.vessels, self.config_handler.simparams.viscosity, self.result_handler.results['preop'], branch, steady=steady))

        return wss
    

    def plot_vs_p(self, fig_dir: str, steady: bool=True):
        '''
        plot the pressure vs wall shear stress for the network

        :param fig_dir: the directory to save the figure to
        :param steady: whether to use steady or unsteady results
        '''
        
        p_mean = self.get_network_pressures(steady=steady)[2]

        wss = self.get_network_wss(steady=steady)
        diameters = self.seg_data.segment_diameter

        fig, ax = plt.subplots(2, sharex=True)
        
        ax[0].scatter(p_mean, wss)
        ax[1].scatter(p_mean, diameters)
        plt.xlabel('Mean Pressure (mmHg)')
        ax[0].set_ylabel('Wall Shear Stress (Pa)')
        ax[1].set_ylabel('Diameter (mm)')
        plt.suptitle('WSS and diameter plotted against mean pressure')

        if os.path.exists(fig_dir) == False:
            os.mkdir(fig_dir)

        plt.savefig(fig_dir + 'plot_vs_pressure.png')
