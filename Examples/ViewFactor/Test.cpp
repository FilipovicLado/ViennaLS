#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>

#include <hrleCoordType.hpp>

#include <lsMesh.hpp>
#include <lsSmartPointer.hpp>
#include <lsVTKWriter.hpp>

#include "../BoschDistribution.hpp"

// template<class T, unsigned D>
std::ostream &operator<<(std::ostream &os,
                         const std::array<hrleCoordType, 3> &arr) {
  os << "[";
  os << std::fixed << std::setprecision(2);
  for (unsigned i = 0; i < arr.size() - 1; ++i) {
    os << arr[i] << ", ";
  }
  os << *(arr.end()) << "]";
  return os;
}

template <class T, class U>
double calcDist(const U &dist, const T &origin, const T &coord) {
  double distance = dist.getSignedDistance(origin, coord);
  // std::cout << "Testing: " << coord << " = " << distance << std::endl;
  return distance;
}

double sum(double a, double b, double n) {
  return a * (1 - std::pow(b, n)) / (1 - b);
}

int main() {
  // // series tests to check calculation
  // {
  //   std::cout << std::fixed << std::setprecision(4);
  //   double z_0 = 0;
  //   double z_n = 0;
  //   double z_old = 0;
  //   const double d = 3.0;
  //   const double D = 50.;
  //   const double b = 1 - d/D; //1 / (1 + d / D);
  //   const double a = d; //d * b;
  //   const double denom = std::log(1.0 - d/D);
  //   for(unsigned i = 0; i <= 100; ++i) {
  //     // double z = sum(a, b, i);
  //     double error_z = 0.99 * z_n;
  //     double n = std::log(1 - error_z/D) / denom;

  //     std::cout << "n" << i << ": " << std::endl;
  //     std::cout << "calculated n: " << n << std::endl;
  //     double z = sum(a, b, std::round(n));

  //     std::cout << "z_n: " << z_n << std::endl;
  //     std::cout << "z: " << z << std::endl;

  //     std::cout << "Error:      " << std::abs(z_n - z) << std::endl;
  //     // std::scientific << "dn: " << (double(i)/n - 1);
  //     std::cout << std::fixed << std::endl;
  //     z_old = z_n;
  //     z_n = a + b * z_old;
  //   }

  //   return 0;
  // }

  // check distribution
  {
    constexpr double gridDelta = 0.5;
    constexpr double radius = 3.0;
    const double trenchBottom = -120;
    const double taperStart = -60;
    BoschDistribution<double, 2> dist(-radius, taperStart, trenchBottom,
                                      gridDelta, -0.6, 10.0);

    std::cout << std::fixed << std::setprecision(4);
    for (double z = -10; z > -120; z -= 0.5) {
      // std::cout << "z: " << z << "--> " <<
      double radius = dist.getRadius(z); // << std::endl;
      // std::cout << "radius: " << radius << std::endl;
    }

    return 0;
  }

  const std::array<hrleCoordType, 3> origin{};
  constexpr double gridDelta = 0.5;
  constexpr double radius = 3.0;
  const double trenchBottom = -238;
  const double taperStart = trenchBottom * 1.5;
  BoschDistribution<double, 2> dist(-radius, taperStart, trenchBottom,
                                    gridDelta, -0.6, 10.0);
  // lsSphereDistribution<double, 2> dist(3.0, gridDelta);

  std::array<hrleCoordType, 3> coord{};

  auto mesh = lsSmartPointer<lsMesh>::New();
  std::vector<double> distanceData;

  for (double y = -2.5 * radius; y <= 2.5 * radius; y += 0.5) {
    std::cout << "y = " << y << std::endl;
    coord[1] = y;
    for (double x = -1.5 * radius; x < 1.5 * radius; x += 0.5) {
      coord[0] = x;
      double distance = calcDist(dist, origin, coord);

      if (std::abs(y) < 0.1) {
        std::cout << coord << ": " << distance << std::endl;
      }

      mesh->insertNextNode(coord);
      std::array<unsigned, 1> vertex{
          static_cast<unsigned int>(mesh->vertices.size())};
      mesh->insertNextVertex(vertex);
      if (std::abs(distance) > 3 * gridDelta) {
        distance = 3 * gridDelta * (std::signbit(distance) ? -1 : 1);
      }
      distanceData.push_back(distance);
    }
  }

  mesh->insertNextScalarData(distanceData, "LSValues");

  lsVTKWriter(mesh, "testMesh.vtk").apply();

  return 0;
}