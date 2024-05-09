#include <iostream>
#include <iomanip>
#include <cassert>
#include <fstream>

#include <omp.h>

#include <Eigen/Dense>
// #include <Eigen/Eigenvalues>

#include "typedefs.hpp"
#include "constants.hpp"
#include "header.hpp"

// ======================================



// int main(){
int main(int argc, char **argv){
#ifdef _OPENMP
  omp_set_num_threads( nparallel );
#endif

  if (argc>1){
    nu = atoi(argv[1]);
    // printf("%s\n", argv[i]);
  }

  const std::string description = "Lx"+std::to_string(Lx)+"Ly"+std::to_string(Ly)+"nu"+std::to_string(nu);

  set_all();

  std::cout << "ell0 = " << ell0[0] << ", " << ell0[1] << std::endl
            << "ell1 = " << ell1[0] << ", " << ell1[1] << std::endl
            << "ell2 = " << ell2[0] << ", " << ell2[1] << std::endl;

  std::cout << "ell = " << ell[0] << ", " << ell[1] << ", " << ell[2] << std::endl;

  std::cout << "kappa = " << kappa[0] << ", " << kappa[1] << ", " << kappa[2] << std::endl;

  std::cout << "ell0* = " << ell_star0[0] << ", " << ell_star0[1] << std::endl
            << "ell1* = " << ell_star1[0] << ", " << ell_star1[1] << std::endl
            << "ell2* = " << ell_star2[0] << ", " << ell_star2[1] << std::endl;

  std::cout << "e0 = " << e0[0] << ", " << e0[1] << std::endl
            << "e1 = " << e1[0] << ", " << e1[1] << std::endl
            << "e2 = " << e2[0] << ", " << e2[1] << std::endl;


  {
    Vect init = Eigen::VectorXcd::Zero(2*Lx*Ly);
    const int xx = 0, yy = 0;

    Vect e0 = Eigen::VectorXcd::Zero(2*Lx*Ly);
    e0( 2*idx(xx, yy) ) = 1.0;
    e0 = multDdagger_eigen(e0);
    // e0 = multD_eigen(e0);

    std::cout << "e0 = " << e0 << std::endl;

    // Eigen::MatrixXcd D = get_Dirac_matrix();
    // // std::cout << "D = " << std::endl
    // //           << D << std::endl;


    // e0 = Eigen::VectorXcd::Zero(2*Lx*Ly);
    // e0( 2*idx(xx, yy) ) = 1.0;
    // e0 = D*e0;
    // std::cout << "e0 = " << e0.transpose() << std::endl;


    // e0 = Eigen::VectorXcd::Zero(2*Lx*Ly);
    // e0( 2*idx(xx, yy) ) = 1.0;
    // // e0 = D*e0;
    // e0 = multDdagger_eigen(e0);
    // std::cout << "e0 = " << e0.transpose() << std::endl;


    // e0 = Eigen::VectorXcd::Zero(2*Lx*Ly);
    // e0( 2*idx(xx, yy) ) = 1.0;
    // e0 = D.adjoint()*e0;
    // std::cout << "e0 = " << e0.transpose() << std::endl;


    // return 1;

    const Vect Dinv0 = CG(init, e0);

    std::cout << "Dinv0 = " << Dinv0 << std::endl;
    // std::cout << "D Dinv0 = " << multD_eigen(Dinv0) << std::endl;

    {
      std::ofstream of( dir_data+description+"Dinv_0_0_0.dat",
                        std::ios::out | std::ios::binary | std::ios::trunc);
      if(!of) assert(false);

      double tmp = 0.0;
      for(Idx i=0; i<2*Lx*Ly; i++){
        tmp = Dinv0(i).real();
        of.write((char*) &tmp, sizeof(double) );

        tmp = Dinv0(i).imag();
        of.write((char*) &tmp, sizeof(double) );
      }
    }

    Vect e1 = Eigen::VectorXcd::Zero(2*Lx*Ly);
    e1( 2*idx(xx, yy)+1 ) = 1.0;
    e1 = multDdagger_eigen(e1);
    const Vect Dinv1 = CG(init, e1);

    {
      std::ofstream of( dir_data+description+"Dinv_0_0_1.dat",
                        std::ios::out | std::ios::binary | std::ios::trunc);
      if(!of) assert(false);

      double tmp = 0.0;
      for(Idx i=0; i<2*Lx*Ly; i++){
        tmp = Dinv1(i).real();
        of.write((char*) &tmp, sizeof(double) );

        tmp = Dinv1(i).imag();
        of.write((char*) &tmp, sizeof(double) );
      }

    }
  }



  {
    Vect init = Eigen::VectorXcd::Zero(2*Lx*Ly);
    const int xx = -1, yy = 0;

    Vect e0 = Eigen::VectorXcd::Zero(2*Lx*Ly);
    e0( 2*idx(xx, yy) ) = 1.0;
    e0 = multDdagger_eigen(e0);

    // std::cout << "e0 = " << e0.transpose() << std::endl;

    const Vect Dinv0 = CG(init, e0);

    // std::cout << "Dinv0 = " << Dinv0 << std::endl;
    // std::cout << "D Dinv0 = " << multD_eigen(Dinv0) << std::endl;

    {
      std::ofstream of( dir_data+description+"Dinv_m1_0_0.dat",
                        std::ios::out | std::ios::binary | std::ios::trunc);
      if(!of) assert(false);

      double tmp = 0.0;
      for(Idx i=0; i<2*Lx*Ly; i++){
        tmp = Dinv0(i).real();
        of.write((char*) &tmp, sizeof(double) );

        tmp = Dinv0(i).imag();
        of.write((char*) &tmp, sizeof(double) );
      }
    }

    Vect e1 = Eigen::VectorXcd::Zero(2*Lx*Ly);
    e1( 2*idx(xx, yy)+1 ) = 1.0;
    e1 = multDdagger_eigen(e1);
    const Vect Dinv1 = CG(init, e1);

    {
      std::ofstream of( dir_data+description+"Dinv_m1_0_1.dat",
                        std::ios::out | std::ios::binary | std::ios::trunc);
      if(!of) assert(false);

      double tmp = 0.0;
      for(Idx i=0; i<2*Lx*Ly; i++){
        tmp = Dinv1(i).real();
        of.write((char*) &tmp, sizeof(double) );

        tmp = Dinv1(i).imag();
        of.write((char*) &tmp, sizeof(double) );
      }

    }
  }


  {
    Vect init = Eigen::VectorXcd::Zero(2*Lx*Ly);
    const int xx = 1, yy = -1;

    Vect e0 = Eigen::VectorXcd::Zero(2*Lx*Ly);
    e0( 2*idx(xx, yy) ) = 1.0;
    e0 = multDdagger_eigen(e0);

    // std::cout << "e0 = " << e0.transpose() << std::endl;

    const Vect Dinv0 = CG(init, e0);

    // std::cout << "Dinv0 = " << Dinv0 << std::endl;
    // std::cout << "D Dinv0 = " << multD_eigen(Dinv0) << std::endl;

    {
      std::ofstream of( dir_data+description+"Dinv_1_m1_0.dat",
                        std::ios::out | std::ios::binary | std::ios::trunc);
      if(!of) assert(false);

      double tmp = 0.0;
      for(Idx i=0; i<2*Lx*Ly; i++){
        tmp = Dinv0(i).real();
        of.write((char*) &tmp, sizeof(double) );

        tmp = Dinv0(i).imag();
        of.write((char*) &tmp, sizeof(double) );
      }
    }

    Vect e1 = Eigen::VectorXcd::Zero(2*Lx*Ly);
    e1( 2*idx(xx, yy)+1 ) = 1.0;
    e1 = multDdagger_eigen(e1);
    const Vect Dinv1 = CG(init, e1);

    {
      std::ofstream of( dir_data+description+"Dinv_1_m1_1.dat",
                        std::ios::out | std::ios::binary | std::ios::trunc);
      if(!of) assert(false);

      double tmp = 0.0;
      for(Idx i=0; i<2*Lx*Ly; i++){
        tmp = Dinv1(i).real();
        of.write((char*) &tmp, sizeof(double) );

        tmp = Dinv1(i).imag();
        of.write((char*) &tmp, sizeof(double) );
      }

    }
  }


  {
    Vect init = Eigen::VectorXcd::Zero(2*Lx*Ly);
    const int xx = 0, yy = 1;

    Vect e0 = Eigen::VectorXcd::Zero(2*Lx*Ly);
    e0( 2*idx(xx, yy) ) = 1.0;
    e0 = multDdagger_eigen(e0);

    // std::cout << "e0 = " << e0.transpose() << std::endl;

    const Vect Dinv0 = CG(init, e0);

    // std::cout << "Dinv0 = " << Dinv0 << std::endl;
    // std::cout << "D Dinv0 = " << multD_eigen(Dinv0) << std::endl;

    {
      std::ofstream of( dir_data+description+"Dinv_0_1_0.dat",
                        std::ios::out | std::ios::binary | std::ios::trunc);
      if(!of) assert(false);

      double tmp = 0.0;
      for(Idx i=0; i<2*Lx*Ly; i++){
        tmp = Dinv0(i).real();
        of.write((char*) &tmp, sizeof(double) );

        tmp = Dinv0(i).imag();
        of.write((char*) &tmp, sizeof(double) );
      }
    }

    Vect e1 = Eigen::VectorXcd::Zero(2*Lx*Ly);
    e1( 2*idx(xx, yy)+1 ) = 1.0;
    e1 = multDdagger_eigen(e1);
    const Vect Dinv1 = CG(init, e1);

    {
      std::ofstream of( dir_data+description+"Dinv_0_1_1.dat",
                        std::ios::out | std::ios::binary | std::ios::trunc);
      if(!of) assert(false);

      double tmp = 0.0;
      for(Idx i=0; i<2*Lx*Ly; i++){
        tmp = Dinv1(i).real();
        of.write((char*) &tmp, sizeof(double) );

        tmp = Dinv1(i).imag();
        of.write((char*) &tmp, sizeof(double) );
      }

    }
  }



  return 0;
}

