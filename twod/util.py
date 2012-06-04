import utfs
import re

def convert_to_real(ss):
    return float(str(ss).replace('d','e'))

def convert_to_int(ss):
    if (ss == None): return None
    return int(ss)

def read_layers(fn):
    dd = dict()
    f = open(fn,'r')
    num_layers = int(re.match('(\d+).*', f.readline()).groups()[0])
    #print 'num_layers', num_layers
    layers = list()
    for x in range(1,num_layers+1):
        (eps,height,is_cond) = re.match('([\-\d\.dDeE]+)\s+([\-\d\.dDeE]+)\s+([tTfF]).*', f.readline()).groups()
        (eps,height) = map(convert_to_real, (eps,height))
        if (str(is_cond).upper() == 'T'):
            is_cond = True
        else:
            is_cond = False
        #print 'eps',eps,'height',height
        layers.append((eps,height,is_cond))
    threshold = convert_to_real(re.match('([\-\d\.dDeE]+).*', f.readline()).groups()[0])    
    #print 'threshold',threshold
    avg_edge_length = convert_to_real(re.match('([\-\d\.dDeE]+).*', f.readline()).groups()[0])
    (xmin,xmax) = map(convert_to_real, re.match('\s*([\-\d\.dDeE]+)\s+([\-\d\.dDeE]+).*', f.readline()).groups())
    xmax *= 2.54e-5
    xmin *= 2.54e-5
    f.close()
    return (layers,avg_edge_length,xmin,xmax,threshold)


def read_geom(fn):
    f = open(fn,'r')
    (num_nodes,num_edges) = map(convert_to_int, re.match('\s*([\d]+)\s+([\d]+).*', f.readline()).groups())
    nodes = list()
    edges = list()
    for nn in range(1,num_nodes+1):
        (xx,yy) = map(convert_to_real, re.match('\s*([\-\d\.dDeE]+)\s+([\-\d\.dDeE]+).*', f.readline()).groups())
        nodes.append((xx*2.54e-5,yy*2.54e-5))

    for nn in range(1,num_edges+1):
        (fromm,to, cid) = map(convert_to_int, re.match('\s*([\d]+)\s+([\d]+)\s*([\d]+)?.*', f.readline()).groups())
        edges.append([fromm,to, cid])
    return (nodes, edges)


def assign_edge_to_cond(nodes, edges):
    for edge in edges:
        if (edge[2] != None): continue
        n1 = nodes[edge[0]-1]
        n2 = nodes[edge[1]-1]
        center = (n1[1] + n2[1]) / 2.0
        if (center<2.2e-5 and center>1.27e-5):
            edge[2] = 1
        elif (center>2.5e-5):
            edge[2] = 2
        else:
            edge[2] = 3
    return edges


def create_layers(layers, avg_edge_length):
    utfs.num_layers(len(layers))
    index = 1
    for ll in layers:
        utfs.set_layer(index, ll[0], ll[1], ll[2])
        index += 1
    utfs.init_layers(avg_edge_length * 2.54e-5)

def add_nodes_and_edges(nodes, edges):
    utfs.num_node_num_edge(len(nodes), len(edges))
    for nn in nodes:
        utfs.add_point(nn[0], nn[1])
    max_cond_id = 0
    for ee in edges:
        utfs.add_edge(ee[0], ee[1],ee[2])
        max_cond_id = max(max_cond_id, ee[2])
    utfs.num_cond(max_cond_id)

