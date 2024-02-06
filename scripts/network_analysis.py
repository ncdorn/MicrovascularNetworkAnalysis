import sys
sys.path.append('/Users/ndorn/Documents/Stanford/PhD/Marsden_Lab/SimVascular/MicrovascularNetworkAnalysis')
from src.ps_network import PSNetwork


def analyze_network(filename):

    dirs = filename.split('/')

    network = PSNetwork.from_file(filename)

    network.config_handler.to_json(filename.replace('.txt', '_config.json'))

    network.run_simulation()

    network.plot_vs_p(dirs[0] + '/' + dirs[1] + '/figs/')

if __name__ == "__main__":

    filename = 'data/546-segments/546-segments.txt'

    analyze_network(filename)