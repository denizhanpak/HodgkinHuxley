import networkx as nx
import matplotlib
from matplotlib import pyplot as plt
import numpy

def colored_circle(network):
    nx.draw_circular(network, node_color = "red")
