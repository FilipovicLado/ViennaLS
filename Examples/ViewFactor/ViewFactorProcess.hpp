#include <array>
#include <unordered_map>

#include <lsDomain.hpp>
#include <lsSmartPointer.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsWriteVisualizationMesh.hpp>

#include "ViewFactorDistribution.hpp"
#include "ViewFactorProcessData.hpp"
#include "lsBisect.hpp"

template <class T, int D> class ViewFactorProcess {
  using LSPtrType = lsSmartPointer<lsDomain<T, D>>;

  LSPtrType substrate;

  ViewFactorProcessDataType<T> processData;

struct hash {
  private:
    std::size_t hash_combine(std::size_t lhs, std::size_t rhs) const {
      lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
      return lhs;
    }

  public:
    std::size_t operator()(const std::array<T, 3> &v) const {
      using std::hash;
      using std::size_t;
      // using std::string;

      /*
        https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
      */
      std::size_t result = hash<T>()(v[0]);
      result = hash_combine(result, hash<T>()(v[1]));
//      if (D == 3) {
        result = hash_combine(result, hash<T>()(v[2]));
//      }
      return result;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:

      // size_t result = hash<T>()(v[0]);
      // result ^= hash<T>()(v[1]) << 1;
      // if (D == 3) {
      //   result = (result >> 1) ^ (hash<T>()(v[2]) << 1);
      // }
      // return result;
    }
  };

public:

  ViewFactorProcess() {}

  ViewFactorProcess(LSPtrType passedSubstrate)
      : substrate(passedSubstrate) {
    processData.gridDelta = substrate->getGrid().getGridDelta();
  }

  void setSubstrate(LSPtrType levelSet) {
    substrate = levelSet;
    processData.gridDelta = substrate->getGrid().getGridDelta();
  }
  
  void setTrenchDiameter(T trenchDiameter) {
    processData.trenchDiameter = trenchDiameter;
  }
  
  void setTrenchDepth(T trenchDepth) {
    processData.trenchDepth = trenchDepth;
  }

  void setTaperAngle(hrleVectorType<double, 2> taperAngle) {
    processData.taperAngle = taperAngle;
  }

  void setTopRate(T topRate) {
    processData.topRate = topRate;
  }
  
  void setProcessTime(T processTime) {
    processData.processTime = processTime;
  }

  void setTimeStep(T timeStep) {
    processData.timeStep = timeStep;
  }
  
  void apply() {
    using mapType = std::unordered_map<std::array<hrleCoordType,3>, double, hash>;
    mapType map;
    auto dist = lsSmartPointer<ViewFactorDistribution<T, D, mapType>>::New(processData, map);
    lsGeometricAdvect<T, D>(substrate, dist).apply();
    
  }

};