#pragma once
#include "CageInitializer/CageInitializer.hh"
#include "CageSimplifier/CageSimplifier.hh"

namespace Cage
{
using namespace SimpleUtils;
using namespace Geometry;
using namespace CageInit;
using namespace CageSimp;
namespace SM = SurfaceMesh;
namespace VM = VolumeMesh;

class CageGenerator
{
public:
  // input
  std::unique_ptr<SM::SMeshT> originalMesh;
  ParamCageGenerator param;

  // sub components
  std::unique_ptr<CageInitializer> cageInitializer;
  std::unique_ptr<CageSimplifier> cageSimplifier;

  // middle results
  std::unique_ptr<VM::VMeshT> VMesh;

  // result
  std::unique_ptr<SM::SMeshT> cage;
public:
  void stageInitialize();
  void stageSimplify();

  void generate();
};


}// namespace Cage