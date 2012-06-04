import sys
sys.path.append('/home/altan/dev/blue-bird/twod')
from util import *
(layers,avg_edge_length,xmin,xmax,tol) = read_layers('/home/altan/dev/blue-bird/twod/layers.inp')
(nodes,edges) = read_geom('/home/altan/dev/blue-bird/twod/geo_pec.inp')

print layers
print nodes
print edges
print xmin, xmax, tol

utfs.init()
utfs.set_tolerance(tol)
utfs.set_x_limits(xmin,xmax)
create_layers(layers, avg_edge_length)

utfs.calculate_green_table()


for x in range(1,1000):
    add_nodes_and_edges(nodes,edges)
    utfs.simulate()
    utfs.clear()
    

