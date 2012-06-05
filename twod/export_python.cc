#include <boost/python.hpp>

extern void init_utfs_core();

namespace bp = boost::python;

int main(int argc, char** argv) {
  Py_Initialize();
  bp::object main = bp::import("__main__");
  bp::object main_namespace = main.attr("__dict__");
  init_utfs_core();
  main_namespace["utfs_core"] = bp::import("utfs_core");
  //bp::exec("from hello import *", main_namespace);
  //bp::exec("a = hello.Goo(); a.doSomething()", main_namespace);
  Py_Main(argc, argv);
  Py_Finalize();
  return 0;
}

