#include <boost/python.hpp>
#include <iostream>

using namespace std;
namespace bp = boost::python;

extern "C" {
  void ky_num_node_num_edge_(int* nnod, int* nedg);
  void ky_num_cond_(int* ncond_tmp);
  void ky_add_edge_(int* from, int* to, int* cid);
  void ky_add_point_(double* xx, double* yy);
  void ky_num_layers_(int* num);
  void ky_set_layer_(int* index, double* eps, double* height, int* is_cond);
  void ky_simulate_();
  void ky_set_tol_(double* tol);
  void ky_init_();
  void ky_set_x_limits_(double* xx, double* xx2);
  void calculate_green_table_();
  void ky_init_layers_(double* len);
  void ky_clear_local_();
  void ky_get_cap_(int* con1, int* cond2, double* val);
  void ky_compute_one_green_(double* src_x,double* src_y, double* obs_x, double* obs_y, double* outt);
};

static void ky_simulate() { ky_simulate_(); }
static void ky_num_node_num_edge(int nnod, int nedg) { ky_num_node_num_edge_(&nnod, &nedg); }
static void ky_num_cond(int ncond_tmp) { ky_num_cond_(&ncond_tmp); }
static void ky_add_point(double xx, double yy) { ky_add_point_(&xx, &yy); }
static void ky_add_edge(int xx, int yy, int cid) { ky_add_edge_(&xx, &yy, &cid); }
static void ky_num_layers(int num) { ky_num_layers_(&num); }
static void ky_set_layer(int index, double eps, double height, int is_cond) { ky_set_layer_(&index, &eps, &height, &is_cond); }
static void ky_set_tol(double tol) { ky_set_tol_(&tol); }
static void ky_init() { ky_init_(); }
static void ky_set_x_limits(double xx, double mx) { ky_set_x_limits_(&xx, &mx); }
static void calculate_green_table() { calculate_green_table_(); }
static void ky_init_layers(double len) { ky_init_layers_(&len); }
static void ky_clear_local() { ky_clear_local_(); }
static double ky_get_C(int cond1, int cond2) { double val = 0; ky_get_cap_(&cond1,&cond2,&val); return val; }

static double ky_compute_one_green(double src_x,double src_y, double obs_x, double obs_y) {
  double outt = 0;
  ky_compute_one_green_(&src_x,&src_y,&obs_x,&obs_y,&outt);
  return outt;
}

BOOST_PYTHON_MODULE(utfs) {
  bp::def("num_node_num_edge", ky_num_node_num_edge);
  bp::def("num_cond", ky_num_cond);
  bp::def("add_point", ky_add_point);
  bp::def("add_edge", ky_add_edge);
  bp::def("simulate", ky_simulate);
  bp::def("num_layers", ky_num_layers);
  bp::def("set_layer", ky_set_layer);
  bp::def("set_tolerance", ky_set_tol);
  bp::def("set_x_limits", ky_set_x_limits);
  bp::def("init", ky_init);
  bp::def("calculate_green_table", calculate_green_table);
  bp::def("init_layers", ky_init_layers);
  bp::def("clear", ky_clear_local);
  bp::def("get_cap", ky_get_C);
  bp::def("compute_one_green", ky_compute_one_green);
};

void init_utfs() {
  initutfs();
}

