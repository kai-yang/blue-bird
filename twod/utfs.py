import utfs_core

def create_layers(layers, avg_edge_length, xmin, xmax, tol):
    utfs_core.init()
    utfs_core.set_tolerance(tol)
    utfs_core.set_x_limits(xmin,xmax)
    utfs_core.num_layers(len(layers))
    index = 1
    for ll in layers:
        utfs_core.set_layer(index, ll[0], ll[1], ll[2])
        index += 1
    utfs_core.init_layers(avg_edge_length * 2.54e-5)

def add_nodes_and_edges(nodes, edges):
    utfs_core.num_node_num_edge(len(nodes), len(edges))
    for nn in nodes:
        utfs_core.add_point(nn[0], nn[1])
    max_cond_id = 0
    for ee in edges:
        utfs_core.add_edge(ee[0], ee[1],ee[2])
        max_cond_id = max(max_cond_id, ee[2])
    utfs_core.num_cond(max_cond_id)
    return max_cond_id

def simulate():
    utfs_core.simulate()

def get_cap(ii,jj):
    return utfs_core.get_cap(ii,jj)

def clear():
    utfs_core.clear();

def calculate_green_table():
    utfs_core.calculate_green_table()
