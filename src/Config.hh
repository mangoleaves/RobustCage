#pragma once
#include <vector>
#include <string>
#include "boost/json.hpp"

namespace Cage
{
struct ParamLatticePointsGenerator
{
  /// minimal point's distance to surface mesh, nearer points will be ingored.
  double minDistanceFactor;

  /***** used for "Uniform" *****/

  /// Scale mesh's bounding box to a larger one.
  double bboxScale;
  /// Sample n^3 lattice points in bounding box uniformly.
  size_t pointNumAlongAxis;

  /***** used for "Offset" *****/

  /// good approximation: 4, high efficiency: 16.
  double areaThresholdRate;
};

struct ParamTetrahedralizer
{
  ParamLatticePointsGenerator paramLatticePointsGenerator;
};

struct ParamCageInitializer
{
  std::string fileOutPath;
  std::string fileName;

  ParamTetrahedralizer paramTetrahedralizer;
};

struct ParamFastSimplifier
{
  std::string fileOutPath;
  std::string fileName;

  // target
  size_t targetVerticesNum;

  // parameters
  size_t collapseIter;
  size_t splitIter;
  size_t equalizeValenceIter;
  size_t smoothIter;
  double enlargeTargetLengthRatio;
};

struct ParamCollapseStage
{
  // constraints
  size_t maxValence;

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["maxValence"] = maxValence;
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    maxValence = jo.at("maxValence").as_int64();
  }
};

struct ParamRelocateStage
{
  // simplification parameter
  size_t smoothIter;
};

struct ParamFlipStage
{
  // constraints
  size_t maxValence;

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["maxValence"] = maxValence;
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    maxValence = jo.at("maxValence").as_int64();
  }
};

struct ParamCageSimplifier
{
  // target
  size_t targetVerticesNum;
  // iterations
  size_t maxIter;
  // distance error control
  // error relax iter->Hausdorff Distance: 0->no negtive, 1->initError, 2->initError+errorStep, 3->initError+errorStep*2, ...
  size_t relaxErrorIterStep;
  size_t maxErrorRelaxIter;
  double initError;
  double errorStep;

  // file output
  size_t cageLabel;
  std::string fileOutPath;
  std::string fileName;

  ParamFastSimplifier paramFastSimplifier;
  ParamCollapseStage paramCollapse;
  ParamRelocateStage paramRelocate;
  ParamFlipStage paramFlip;

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["maxIter"] = maxIter;
    jo["relaxErrorIterStep"] = relaxErrorIterStep;
    jo["maxErrorRelaxIter"] = maxErrorRelaxIter;
    jo["initError"] = initError;
    jo["errorStep"] = errorStep;
    jo["paramCollapse"] = paramCollapse.serialize();
    jo["paramFlip"] = paramFlip.serialize();
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    maxIter = jo.at("maxIter").as_int64();
    relaxErrorIterStep = jo.at("relaxErrorIterStep").as_int64();
    maxErrorRelaxIter = jo.at("maxErrorRelaxIter").as_int64();
    initError = jo.at("initError").as_double();
    errorStep = jo.at("errorStep").as_double();
    paramCollapse.deserialize(jo.at("paramCollapse").as_object());
    paramFlip.deserialize(jo.at("paramFlip").as_object());
  }
};

struct ParamCageGenerator
{
  ParamCageInitializer paramCageInitializer;
  ParamCageSimplifier paramCageSimplifier;

  ParamCageGenerator()
  {
    auto& Lpg = paramCageInitializer.paramTetrahedralizer.paramLatticePointsGenerator;
    Lpg.minDistanceFactor = 0.01;
    Lpg.bboxScale = 1.5;
    Lpg.pointNumAlongAxis = 10;
    Lpg.areaThresholdRate = 16;

    auto& simplifier = paramCageSimplifier;
    simplifier.maxIter = 30;
    simplifier.relaxErrorIterStep = 5;
    simplifier.maxErrorRelaxIter = 4;
    simplifier.initError = 0.005;
    simplifier.errorStep = 0.005;
    simplifier.targetVerticesNum = 0;

    auto& fast = paramCageSimplifier.paramFastSimplifier;
    fast.collapseIter = 3;
    fast.splitIter = 1;
    fast.equalizeValenceIter = 1;
    fast.smoothIter = 3;
    fast.enlargeTargetLengthRatio = 1.1;
    fast.targetVerticesNum = 0;

    auto& collapse = paramCageSimplifier.paramCollapse;
    collapse.maxValence = 8;

    auto& relocate = paramCageSimplifier.paramRelocate;
    relocate.smoothIter = 3;

    auto& flip = paramCageSimplifier.paramFlip;
    flip.maxValence = 8;
  }

  void setOutputPath(const std::string& outDir, const std::string& outFile)
  {
    paramCageInitializer.fileOutPath = outDir;
    paramCageInitializer.fileName = outFile;
    paramCageSimplifier.fileOutPath = outDir;
    paramCageSimplifier.fileName = outFile;
    paramCageSimplifier.paramFastSimplifier.fileOutPath = outDir;
    paramCageSimplifier.paramFastSimplifier.fileName = outFile;
  }

  void setCageLabel(const size_t label)
  {
    paramCageSimplifier.cageLabel = label;
  }
  void setTargetNumber(const size_t target)
  {
    paramCageSimplifier.targetVerticesNum = target;
  }
  void setFastTargetNumber(const size_t target)
  {
    paramCageSimplifier.paramFastSimplifier.targetVerticesNum = target;
  }

  boost::json::object serialize()const
  {
    boost::json::object jo;
    jo["paramCageSimplifier"] = paramCageSimplifier.serialize();
    return jo;
  }
  void deserialize(const boost::json::object& jo)
  {
    paramCageSimplifier.deserialize(jo.at("paramCageSimplifier").as_object());
  }
};

inline void tag_invoke(boost::json::value_from_tag, boost::json::value& jv, const ParamCageGenerator& param)
{
  auto& jo = jv.emplace_object();
  jo["paramCageGenerator"] = param.serialize();
}

inline ParamCageGenerator tag_invoke(boost::json::value_to_tag<ParamCageGenerator>, boost::json::value& jv)
{
  ParamCageGenerator param;
  param.deserialize(jv.as_object());
  return param;
}
}