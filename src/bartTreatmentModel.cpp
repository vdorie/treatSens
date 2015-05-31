#include "bartTreatmentModel.hpp"

#include <new>
#include <dbarts/cstdint.hpp>

#include <dbarts/bartFit.hpp>
#include <dbarts/control.hpp>
#include <dbarts/data.hpp>
#include <dbarts/model.hpp>
#include <dbarts/results.hpp>
#include <dbarts/types.hpp>

#include <external/random.h>
#include <external/stats.h>

using std::size_t;
using std::uint32_t;

#define DEFAULT_BART_MAX_NUM_CUTS 100u

namespace {
  using cibart::TreatmentModel;
  
  void* createScratch(TreatmentModel* restrict model, ext_rng* restrict generator, const double* restrict x, size_t numObservations, size_t numPredictors, const double* restrict z);
  void destroyScratch(TreatmentModel* model, void* scratch);
  void updateParameters(TreatmentModel* restrict model, void* restrict scratch, const double* restrict offset);
  void getConditionalProbabilities(TreatmentModel* restrict model, void* restrict scratch, double zetaZ, double* restrict probU0, double* restrict probU1);
  
  void lookupBARTFunctions(cibart::BARTTreatmentModelFunctionTable* functionTable, cibart::voidPtrFunctionLookup functionLookup);
}

namespace cibart {
  struct BARTTreatmentModelFunctionTable {
    void (*initializeFit)(dbarts::BARTFit* fit, dbarts::Control* control, dbarts::Model* model, dbarts::Data* data);
    void (*invalidateFit)(dbarts::BARTFit* fit);
    dbarts::Results* (*runSampler)(dbarts::BARTFit* fit);
    void (*setOffset)(dbarts::BARTFit* fit, const double* offset);
    void (*initializeCGMPrior)(dbarts::CGMPrior* prior, double, double);
    void (*invalidateCGMPrior)(dbarts::CGMPrior* prior);
    void (*initializeNormalPrior)(dbarts::NormalPrior* prior, const dbarts::Control* control, double);
    void (*invalidateNormalPrior)(dbarts::NormalPrior* prior);
    void (*initializeChiSquaredPrior)(dbarts::ChiSquaredPrior* prior, double, double);
    void (*invalidateChiSquaredPrior)(dbarts::ChiSquaredPrior* prior);
  };
  
  BARTTreatmentModel::BARTTreatmentModel(voidPtrFunctionLookup lookup, size_t numTrees, size_t numThin, double nodePriorParameter) :
    numTrees(numTrees), numThin(numThin), nodePriorParameter(nodePriorParameter)
  {
    predictorsIncludeIntercept = false;
    includesLatentVariables = false;
    this->createScratch = &::createScratch;
    this->destroyScratch = &::destroyScratch;
    this->updateParameters = &::updateParameters;
    this->getConditionalProbabilities = &::getConditionalProbabilities;
    
    functionTable = new BARTTreatmentModelFunctionTable;
    lookupBARTFunctions(functionTable, lookup);
  }
  
  BARTTreatmentModel::~BARTTreatmentModel() {
    delete functionTable;
  }
}

namespace {
  using cibart::TreatmentModel;
  using cibart::BARTTreatmentModel;
  typedef cibart::BARTTreatmentModelFunctionTable FunctionTable;
  
  struct Scratch {
    dbarts::BARTFit* bartFit;
    dbarts::Control* bartControl;
    dbarts::Data*    bartData;
    dbarts::Model*   bartModel;
    
    dbarts::Results* bartResults;
  };
  
  dbarts::Control* createBARTControl(BARTTreatmentModel& model, ext_rng* generator);
  dbarts::Data*    createBARTData(const double* x, size_t numObservations, size_t numPredictors, const double* z);
  dbarts::Model*   createBARTModel(BARTTreatmentModel& model, const dbarts::Control* bartControl);
  
  void destroyBARTControl(dbarts::Control* control);
  void destroyBARTData(dbarts::Data* data);
  void destroyBARTModel(BARTTreatmentModel& model, dbarts::Model* bartModel);

  void* createScratch(TreatmentModel* restrict modelPtr, ext_rng* restrict generator, const double* restrict x, size_t numObservations, size_t numPredictors, const double* restrict z)
  {
    BARTTreatmentModel& model(*static_cast<BARTTreatmentModel*>(modelPtr));
    FunctionTable& functionTable(*model.functionTable);
    
    Scratch* scratch = new Scratch;
    
    scratch->bartControl = createBARTControl(model, generator);
        
    scratch->bartData = createBARTData(x, numObservations, numPredictors, z);
    
    scratch->bartModel = createBARTModel(model, scratch->bartControl);
    
    scratch->bartFit = static_cast<dbarts::BARTFit*>(::operator new (sizeof(dbarts::BARTFit)));
    functionTable.initializeFit(scratch->bartFit, scratch->bartControl, scratch->bartModel, scratch->bartData);
    
    scratch->bartResults = NULL;
    
    return scratch;
  }
  
  void destroyScratch(TreatmentModel* modelPtr, void* scratchPtr) {
    BARTTreatmentModel& model(*static_cast<BARTTreatmentModel*>(modelPtr));
    FunctionTable& functionTable(*model.functionTable);
    
    Scratch* scratch = static_cast<Scratch*>(scratchPtr);
    
    if (scratch != NULL) {
      delete scratch->bartResults;

      functionTable.invalidateFit(scratch->bartFit);
      ::operator delete(scratch->bartFit);
      
      destroyBARTModel(model, scratch->bartModel);
      destroyBARTData(scratch->bartData);
      destroyBARTControl(scratch->bartControl);
      
      delete scratch;
    }
  }

  void updateParameters(TreatmentModel* restrict modelPtr, void* restrict scratchPtr, const double* restrict offset)
  {
    FunctionTable& functionTable(*(static_cast<BARTTreatmentModel*>(modelPtr)->functionTable));
    Scratch* scratch = static_cast<Scratch*>(scratchPtr);
    
        
    functionTable.setOffset(scratch->bartFit, offset);
    
    scratch->bartResults = functionTable.runSampler(scratch->bartFit);
    
    size_t numObservations = scratch->bartData->numObservations;
    // subtract out offset from prediction
    for (size_t i = 0; i < numObservations; ++i) {
      scratch->bartResults->trainingSamples[i] -= offset[i];
    }
  }
  
  void getConditionalProbabilities(TreatmentModel* restrict, void* restrict scratchPtr, double zetaZ, double* restrict probZForU0, double* restrict probZForU1)
  {
    Scratch* scratch = static_cast<Scratch*>(scratchPtr);
    
    const double* restrict z = scratch->bartData->y;
    const double* restrict zHat = scratch->bartResults->trainingSamples;
    
    size_t numObservations = scratch->bartData->numObservations;
    for (size_t i = 0; i < numObservations; ++i) {
      double zHatU0 = ext_cumulativeProbabilityOfNormal(zHat[i], 0.0, 1.0);
      double zHatU1 = ext_cumulativeProbabilityOfNormal(zHat[i] + zetaZ, 0.0, 1.0);
      
      probZForU0[i] = (z[i] == 1.0 ? zHatU0 : 1.0 - zHatU0);
      probZForU1[i] = (z[i] == 1.0 ? zHatU1 : 1.0 - zHatU1);
    }
    
    delete scratch->bartResults;
    scratch->bartResults = NULL;
  }
  
  
  void lookupBARTFunctions(FunctionTable* functionTable, cibart::voidPtrFunctionLookup lookupFunction)
  {
    functionTable->initializeFit             = reinterpret_cast<void (*)(dbarts::BARTFit*, dbarts::Control*, dbarts::Model*, dbarts::Data*)>(lookupFunction("dbarts", "initializeFit"));
    functionTable->invalidateFit             = reinterpret_cast<void (*)(dbarts::BARTFit*)>(lookupFunction("dbarts", "invalidateFit"));
    functionTable->runSampler                = reinterpret_cast<dbarts::Results* (*)(dbarts::BARTFit*)>(lookupFunction("dbarts", "runSampler"));
    functionTable->setOffset                 = reinterpret_cast<void (*)(dbarts::BARTFit*, const double*)>(lookupFunction("dbarts", "setOffset"));
    functionTable->initializeCGMPrior        = reinterpret_cast<void (*)(dbarts::CGMPrior*, double, double)>(lookupFunction("dbarts", "initializeCGMPriorFromOptions"));
    functionTable->invalidateCGMPrior        = reinterpret_cast<void (*)(dbarts::CGMPrior*)>(lookupFunction("dbarts", "invalidateCGMPrior"));
    functionTable->initializeNormalPrior     = reinterpret_cast<void (*)(dbarts::NormalPrior*, const dbarts::Control*, double)>(lookupFunction("dbarts", "initializeNormalPriorFromOptions"));
    functionTable->invalidateNormalPrior     = reinterpret_cast<void (*)(dbarts::NormalPrior*)>(lookupFunction("dbarts", "invalidateNormalPrior"));
    functionTable->initializeChiSquaredPrior = reinterpret_cast<void (*)(dbarts::ChiSquaredPrior*, double, double)>(lookupFunction("dbarts", "initializeChiSquaredPriorFromOptions"));
    functionTable->invalidateChiSquaredPrior = reinterpret_cast<void (*)(dbarts::ChiSquaredPrior*)>(lookupFunction("dbarts", "invalidateChiSquaredPrior"));
  }
  
  dbarts::Control* createBARTControl(BARTTreatmentModel& model, ext_rng* generator)
  {
    dbarts::Control* control = new dbarts::Control;
    control->numSamples = 1;
    control->numBurnIn  = 0;
    control->numTrees   = model.numTrees;
    control->treeThinningRate = static_cast<uint32_t>(model.numThin);
    control->responseIsBinary = true;
    control->verbose = false;
    control->numThreads = 1;
    control->rng = generator;
    
    return control;
  }
  
  void destroyBARTControl(dbarts::Control* control)
  {
    delete control;
  }
  
  dbarts::Data* createBARTData(const double* x, size_t numObservations, size_t numPredictors, const double* z)
  {
    
    dbarts::VariableType* variableTypes = new dbarts::VariableType[numPredictors];
    for (size_t i = 0; i < numPredictors; ++i) variableTypes[i] = dbarts::ORDINAL;
    
    uint32_t* maxNumCuts = new uint32_t[numPredictors];
    for (size_t i = 0; i < numPredictors; ++i) maxNumCuts[i] = DEFAULT_BART_MAX_NUM_CUTS;
    
    return new dbarts::Data(z, x, NULL, NULL, NULL, NULL,
                            numObservations, numPredictors, 0, 1.0, variableTypes, maxNumCuts);
  }
  
  void destroyBARTData(dbarts::Data* data)
  {
    const dbarts::VariableType* variableTypes = data->variableTypes;
    const uint32_t* maxNumCuts = data->maxNumCuts;
    
    delete data;
    delete [] maxNumCuts;
    delete [] variableTypes;
  }
  
  dbarts::Model* createBARTModel(BARTTreatmentModel& model, const dbarts::Control* bartControl)
  {
    FunctionTable& functionTable(*model.functionTable);
    
    dbarts::Model* bartModel = new dbarts::Model;
    
    bartModel->treePrior = static_cast<dbarts::CGMPrior*>(::operator new (sizeof(dbarts::CGMPrior)));
    functionTable.initializeCGMPrior(static_cast<dbarts::CGMPrior*>(bartModel->treePrior), DBARTS_DEFAULT_TREE_PRIOR_BASE, DBARTS_DEFAULT_TREE_PRIOR_POWER);
    
    bartModel->muPrior = static_cast<dbarts::NormalPrior*>(::operator new (sizeof(dbarts::NormalPrior)));
    functionTable.initializeNormalPrior(static_cast<dbarts::NormalPrior*>(bartModel->muPrior), bartControl, model.nodePriorParameter);
    
    bartModel->sigmaSqPrior = static_cast<dbarts::ChiSquaredPrior*>(::operator new (sizeof(dbarts::ChiSquaredPrior)));
    functionTable.initializeChiSquaredPrior(static_cast<dbarts::ChiSquaredPrior*>(bartModel->sigmaSqPrior), DBARTS_DEFAULT_CHISQ_PRIOR_DF, DBARTS_DEFAULT_CHISQ_PRIOR_QUANTILE);
    
    return bartModel;
  }
  
  void destroyBARTModel(BARTTreatmentModel& model, dbarts::Model* bartModel)
  {
    FunctionTable& functionTable(*model.functionTable);
    
    functionTable.invalidateChiSquaredPrior(static_cast<dbarts::ChiSquaredPrior*>(bartModel->sigmaSqPrior));
    ::operator delete(bartModel->sigmaSqPrior);
    
    functionTable.invalidateNormalPrior(static_cast<dbarts::NormalPrior*>(bartModel->muPrior));
    ::operator delete(bartModel->muPrior);
    
    functionTable.invalidateCGMPrior(static_cast<dbarts::CGMPrior*>(bartModel->treePrior));
    ::operator delete(bartModel->treePrior);
    
    delete bartModel;
  }
}
