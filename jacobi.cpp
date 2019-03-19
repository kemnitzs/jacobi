/*
******************************************************************
* Subroutine HelmholtzJ
* Solves poisson equation on rectangular grid assuming :
* (1) Uniform discretization in each direction, and
* (2) Dirichlect boundary conditions
*
* Jacobi method is used in this routine
*
* Input : n,m   Number of grid points in the X/Y directions
*         dx,dy Grid spacing in the X/Y directions
*         omega Relaxation factor
*         f(n,m) Right hand side function
*         u(n,m) Dependent variable/Solution
*         tolerance Tolerance for iterative solver
*         maxit  Maximum number of iterations
*
* Output : u(n,m) - Solution
*****************************************************************
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <omp.h>
#include "jacobi.h"

void JacobiData::run() {
  std::vector<double> uold_data(n_cols * n_rows);

  auto uold = [&](size_t j, size_t i) {
    return uold_data[((j)-first_row) * n_cols + (i)];
  };

  const double ax = 1.0 / (dx * dx);	 /* X-direction coef */
  const double ay = 1.0 / (dy * dy);	 /* Y_direction coef */
  const double b = -2.0 * (ax + ay); /* Central coeff */
  double residual = 10.0 * tolerance;

  while (effective_iter_count < max_iterations && residual > tolerance) {
    residual = 0.0;
    std::swap(U, uold_data);

#pragma omp parallel reduction(+ : residual)
    {
      // compute stencil, residual and update 
#pragma omp for
      for (int j = first_row + 1; j <= last_row - 1; j++) {
	for (int i = 1; i <= n_cols - 2; i++) {
	  double fLRes = (ax * (uold(j, i - 1) + uold(j, i + 1)) +
			  ay * (uold(j - 1, i) + uold(j + 1, i)) +
			  b * uold(j, i) - get_f(j, i)) /
			 b;

	  // update solution 
	  get_u(j, i) = uold(j, i) - relax * fLRes;

	  // accumulate residual error 
	  residual += pow(fLRes,2);
	}
      }
    } // omp parallel 
    
    // output of solution U every out_iter steps 
    if (effective_iter_count%out_iter == 0 || effective_iter_count == 1){
 	std::string filename = iter2filename("out_",effective_iter_count,".dat"); 
      	JacobiData::out(U,filename);
    }
    // error check 
    effective_iter_count++;
    residual = sqrt(residual) / (n_cols * n_rows);
  } // while 

  std::string last_filename = iter2filename("out_",effective_iter_count,".dat"); 
  JacobiData::out(U,last_filename);

  long latice_site = (last_row - 1 - first_row + 1) * (n_cols - 2 - 1);
  total = latice_site * (effective_iter_count - 1);
  std::cout << "total " << total << std::endl;
  this->residual = residual;
}

JacobiData::JacobiData(){
// default medium 
        n_cols      = 500;
        n_rows      = 500;
        relax     = 1.0;
        tolerance = 1e-10;
        max_iterations   = 1000000;
	out_iter = (int) max_iterations/20;
#ifdef READ_INPUT
        printf("Input n - matrix size in x direction:                 ");
        scanf("%d", &n_cols);
        printf("\nInput m - matrix size in y direction:               ");
        scanf("%d", &n_rows);
        printf("\nInput relax - Successive over-relaxation parameter: ");
        scanf("%lf", &relax);
        printf("\nInput tol - error tolerance for iterrative solver:  ");
        scanf("%lf", &tolerance);
        printf("\nInput mits - Maximum iterations for solver:         ");
        scanf("%d", &max_iterations);
	out_iter= (int) max_iterations/20;
#elif defined DATA_LARGE
        n_cols      = 7000;
        n_rows      = 7000;
        relax     = 1.0;
        tolerance = 1e-12;
        max_iterations   = 2;
	out_iter = 1;
#elif defined DATA_SMALL
        n_cols      = 200;
        n_rows      = 200;
        relax     = 1.0;
        tolerance = 1e-7;
        max_iterations   = 1000;
	out_iter = 50;
#endif
        printf("\n-> matrix size: %dx%d"
               "\n-> relax: %f"
               "\n-> tolerance: %e"
               "\n-> #of iterations: %d \n\n",
               n_cols, n_rows, relax,
               tolerance, max_iterations);


    /* MPI values, set to defaults to avoid data inconsistency */
    my_rank   = 0;
    n_processes = 1;
    first_row = 0;
    last_row  = n_rows - 1;

    /* memory allocation for serial & omp */
    U.resize(n_rows*n_cols);
    F.resize(n_rows*n_cols);

    /* calculate dx and dy */
    dx = 2.0 / (n_cols - 1);
    dy = 2.0 / (n_rows - 1);

    effective_iter_count = 0;
    init_matrix();
}

void JacobiData::init_matrix(){
    /* Initialize initial condition and RHS */

#pragma omp parallel for  
    for (int j = first_row; j <= last_row; j++)
    {
        for (int i = 0; i < n_cols; i++)
        {
            double xx = -1.0 + dx * i;
            double yy = -1.0 + dy * j;

            double xx2 = xx * xx;
            double yy2 = yy * yy;

            get_u(j,i) = 0.0;
	//    get_f(j,i)= - (1.0 - xx2) * (1.0 - yy2) + 2.0 * (-2.0 + xx2 + yy2);
	    get_f(j,i) = - 2.0 * (1.0 - xx2) - 2.0 * (1.0 - yy2);
	//    get_f(j,i) = -M_PI*M_PI*(xx2+yy2)*sin(M_PI*xx*yy);
	}
    }
    max_threads=omp_get_num_threads();	  	
  
    JacobiData::out(F,"charge.dat");
}

void JacobiData::out(std::vector<double> out_array,std::string filename){
      
    std::ofstream out_file;
    out_file.open (filename);
 
    for (int i = 0; i < n_rows; ++i){
      for (int j = 0; j < n_cols; ++j){
	out_file << out_array[j * n_cols + i] << "\t";
      }
      out_file << std::endl;
    } 
    
}

std::string JacobiData::iter2filename(std::string pre_iter, int iter, std::string post_iter){
	  
	std::stringstream ss;	
      	ss << std::setw(6) << std::setfill('0') << iter;
      	std::string step = ss.str();
	std::string filename = pre_iter + step + post_iter;
	
	return filename;	
  	
}
