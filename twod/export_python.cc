#include <boost/python.hpp>

extern void init_utfs();

namespace bp = boost::python;

int main(int argc, char** argv) {
  Py_Initialize();
  bp::object main = bp::import("__main__");
  bp::object main_namespace = main.attr("__dict__");
  init_utfs();
  main_namespace["utfs"] = bp::import("utfs");
  //bp::exec("from hello import *", main_namespace);
  //bp::exec("a = hello.Goo(); a.doSomething()", main_namespace);
  Py_Main(argc, argv);
  Py_Finalize();
  return 0;
}

