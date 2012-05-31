#include <boost/python.hpp>
#include <iostream>

using namespace std;
namespace bp = boost::python;

extern "C" {
  void ky_num_node_num_edge_(int* nnod, int* nedg);
  void ky_add_edge_(int* from, int* to, int* cid);
  void ky_add_point_(double* xx, double* yy);
  void ky_num_layers_(int* num);
  void ky_set_layer_(int* index, double* eps, double* height);
  void ky_simulate_();
};

static void ky_simulate() { ky_simulate_(); }
static void ky_num_node_num_edge(int nnod, int nedg) { ky_num_node_num_edge_(&nnod, &nedg); }
static void ky_add_point(double xx, double yy) { ky_add_point_(&xx, &yy); }
static void ky_add_edge(int xx, int yy, int cid) { ky_add_edge_(&xx, &yy, &cid); }
static void ky_num_layers(int num) { ky_num_layers_(&num); }
static void ky_set_layer(int index, double eps, double height) { ky_set_layer_(&index, &eps, &height); }

BOOST_PYTHON_MODULE(utfs) {
  bp::def("num_node_num_edge", ky_num_node_num_edge);
  bp::def("add_point", ky_add_point);
  bp::def("add_edge", ky_add_edge);
  bp::def("simulate", ky_simulate);
  bp::def("num_layers", ky_num_layers);
  bp::def("set_layer", ky_set_layer);
};

void init_utfs() {
  initutfs();
}

