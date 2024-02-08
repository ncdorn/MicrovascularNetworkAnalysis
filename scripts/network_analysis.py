import sys
sys.path.append('/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/MicrovascularNetworkAnalysis')
from src.ps_network import PSNetwork
from svzerodtrees.utils import m2d


def analyze_network(filename):

    dirs = filename.split('/')

    fig_dir = dirs[0] + '/' + dirs[1] + '/figs/'

    network = PSNetwork.from_file(filename, q_in=0.1, R_bc=100.0)

    network.config_handler.to_json(filename.replace('.txt', '_config.json'))

    # network.run_simulation()

    network.visualize(fig_dir)

    # network.plot_vs_p(fig_dir)

if __name__ == "__main__":

    filename = 'data/546-segments/546-segments.txt'

    analyze_network(filename)