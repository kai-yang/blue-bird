import sys

root = '/home/kyang'
sys.path.append(root + '/blue-bird/twod')
from util import *
(layers,avg_edge_length,xmin,xmax,tol) = read_layers(root + '/blue-bird/twod/layers.inp')
(nodes,edges) = read_geom(root + '/blue-bird/twod/geo_pec.inp')

print layers
print len(nodes)
print len(edges)
print xmin, xmax, tol

utfs.init()
utfs.set_tolerance(tol)
utfs.set_x_limits(xmin,xmax)
create_layers(layers, avg_edge_length)

print 'GF table calculating...'
utfs.calculate_green_table()
print 'GF table done'


for x in range(1,4):
    max_cond_id = add_nodes_and_edges(nodes,edges)
    utfs.simulate()
    for ii in range(1,max_cond_id+1):
        for jj in range(1,max_cond_id+1):
            print ii,jj,utfs.get_cap(ii,jj)
    utfs.clear()
    print x, " done."



