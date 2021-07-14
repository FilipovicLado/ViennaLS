#include <iostream>
#include <chrono>

#include <lsBooleanOperation.hpp>
#include <lsExpand.hpp>
#include <lsGeometricAdvect.hpp>
#include <lsMakeGeometry.hpp>
#include <lsToDiskMesh.hpp>
#include <lsToMesh.hpp>
#include <lsToSurfaceMesh.hpp>
#include <lsVTKWriter.hpp>
#include <lsWriteVisualizationMesh.hpp>

#include "ViewFactorProcess.hpp"

void makeTaperedTrench(lsSmartPointer<lsMesh<>> mesh,
                       hrleVectorType<double, 2> center,
                       hrleVectorType<double, 2> taperAngle,
                       double diameter,
                       double depth)
{
  auto cloud = lsSmartPointer<lsPointCloud<double, 2>>::New();
  {
    //top left
    hrleVectorType<double, 2> point1(-diameter/2., 0.);
    cloud->insertNextPoint(point1);
    //top right
    hrleVectorType<double, 2> point2(diameter/2., 0.);
    cloud->insertNextPoint(point2);
    //bottom right
    hrleVectorType<double, 2> point3(diameter/2.-(depth * taperAngle[1]/taperAngle[0]), -depth);
    cloud->insertNextPoint(point3);
    //bottom left
    hrleVectorType<double, 2> point4(-diameter/2.+(depth * taperAngle[1]/taperAngle[0]), -depth);
    cloud->insertNextPoint(point4);
  }
  lsConvexHull<double, 2>(mesh, cloud).apply();
}


int main() {

  omp_set_num_threads(1);

  constexpr int D = 2;
  typedef double NumericType;
  double gridDelta = 1;
  // Process parameters

  double extent = 20;
  double bounds[2 * D] = {-extent, extent, -3*extent, 3*extent};
  if constexpr (D == 3) {
    bounds[4] = -extent;
    bounds[5] = extent;
  }

  typename lsDomain<NumericType, D>::BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D - 1; ++i) {
    boundaryCons[i] =
        lsDomain<NumericType, D>::BoundaryType::REFLECTIVE_BOUNDARY;
  }
  boundaryCons[D - 1] =
      lsDomain<NumericType, D>::BoundaryType::INFINITE_BOUNDARY;

  hrleVectorType<double, 2> center(0., 0.);
  hrleVectorType<double, 2> normSide(1, 0);
  double diameter = 20.;
  double depth = 40.;

  auto substrate = lsSmartPointer<lsDomain<NumericType, D>>::New(
      bounds, boundaryCons, gridDelta);
  {
    NumericType origin[D] = {0., 0.};
    NumericType planeNormal[D] = {0., 1.};
    auto plane =
        lsSmartPointer<lsPlane<NumericType, D>>::New(origin, planeNormal);
    lsMakeGeometry<NumericType, D>(substrate, plane).apply();
  }
  
  //make LS from trench mesh and remove from substrate
  auto trench = lsSmartPointer<lsDomain<NumericType, D>>::New(
        bounds, boundaryCons, gridDelta);
  {
    //create trench
    auto trenchMesh = lsSmartPointer<lsMesh<>>::New();
    makeTaperedTrench(trenchMesh, center, normSide, diameter, depth);
    lsFromSurfaceMesh<double, D>(trench, trenchMesh, false).apply();
    lsBooleanOperation<double, D>(substrate, trench,
                                  lsBooleanOperationEnum::RELATIVE_COMPLEMENT)
        .apply();
  }
 
  std::cout << "Output initial" << std::endl;
  auto mesh = lsSmartPointer<lsMesh<>>::New();
  lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
  lsVTKWriter(mesh, "surface_i.vtk").apply();
  auto volumeMeshingi = lsSmartPointer<lsWriteVisualizationMesh<NumericType, D>>::New();
  volumeMeshingi->insertNextLevelSet(substrate);
  volumeMeshingi->setFileName("volume_i");
  volumeMeshingi->apply();


  ViewFactorProcess<NumericType, D> processKernel(substrate);
  processKernel.setTrenchDiameter(diameter);
  processKernel.setTrenchDepth(depth);
  processKernel.setTaperAngle(normSide);
  processKernel.setTopRate(1.);
  processKernel.setProcessTime(10.);
  processKernel.setTimeStep(1e-4);
  
  // Run advection
  auto start = std::chrono::high_resolution_clock::now();
  processKernel.apply();
  auto stop = std::chrono::high_resolution_clock::now();
  std::cout << "Geometric advect took: " << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count() << " ms" << std::endl;
  std::cout << "Final structure has " << substrate->getNumberOfPoints() << " LS points" << std::endl;

  lsToSurfaceMesh<NumericType, D>(substrate, mesh).apply();
  lsVTKWriter(mesh, "surface_f.vtk").apply();

  std::cout << "Making volume output..." << std::endl;

  auto volumeMeshingf = lsSmartPointer<lsWriteVisualizationMesh<NumericType, D>>::New();
  volumeMeshingf->insertNextLevelSet(substrate);
  volumeMeshingf->setFileName("volume_f");
  volumeMeshingf->apply();

  return 0;
}
