// Periodic boundary conditions with flexible fibers and particles.

#include "FiberCell2DSimulation.hpp"

int main(int argc, char const *argv[]) {

  FiberCell2DSimulation MD;

  if (argc < 2) {
    MD.result_folder = "./ComprIsoIni";
    MD.check();
  } else if (argc == 2) {
    MD.result_folder = fileTool::GetFilePath(argv[1]) + std::string("/RESULTS");
    MD.loadConf(argv[1]);
    MD.check();
  } else {
    std::cout << "usage: " << argv[0] << " confFile" << std::endl;
  }

  // initial state
  fileTool::create_folder(MD.result_folder);

  MD.saveConf(MD.iconf - 1);
  MD.ResetCloseList(MD.dVerlet);
  MD.integrate();

  return 0;
}
