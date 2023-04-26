#include "CageGenerator.hh"
#include "omp.h"

namespace Cage
{
void CageGenerator::stageInitialize()
{
  VMesh = std::make_unique<VM::VMeshT>();
  cage = std::make_unique<SM::SMeshT>();

  cageInitializer = std::make_unique<CageInitializer>(
    originalMesh.get(), &param.paramCageInitializer,
    VMesh.get(), cage.get());

  cageInitializer->generate();
}

void CageGenerator::stageSimplify()
{
  cageSimplifier = std::make_unique<CageSimplifier>(
    originalMesh.get(), cage.get(), &param.paramCageSimplifier);

  cageSimplifier->simplify();
}

void CageGenerator::generate()
{
  omp_set_num_threads(12);

  stageInitialize();
  cageInitializer = nullptr;
  VMesh = nullptr;
  stageSimplify();
}

}// namespace Cage