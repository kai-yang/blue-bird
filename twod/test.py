import sys

#root = '/home/kyang'
root = '/Users/altan/dev'
sys.path.append(root + '/blue-bird/twod')
import utfs
import utfs_helper
(layers,avg_edge_length,xmin,xmax,tol) = utfs_helper.read_layers(root + '/blue-bird/twod/layers.inp')
(nodes,edges) = utfs_helper.read_geom(root + '/blue-bird/twod/geo_pec.inp')

print layers
print len(nodes)
print len(edges)
print xmin, xmax, tol

utfs.create_layers(layers, avg_edge_length, xmin, xmax, tol)

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



