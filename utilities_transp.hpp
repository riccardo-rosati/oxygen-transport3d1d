/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   utilities_transp.hpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016 - May 2018.
  @brief  Declaration of some miscelleanous auxiliary routines.
 */


#ifndef M3D1D_UTILITIES_TRANSP_HPP_
#define M3D1D_UTILITIES_TRANSP_HPP_

#include <gmm/gmm.h>
#include <defines.hpp>
#include "../utilities/muparser/include/muParser.h"

namespace gmm {

//Calculate the product between two vector, component by component.
//That means: C=A**B   --> C[i]= A[i]* B[i]

//As a matter of fact, it returns: B= A**B
template
<typename VEC1,typename VEC2>
void vscale (const VEC1 & A, VEC2 & B){

//calculate the length of the two vectors

int lengthA = gmm::mat_nrows(gmm::col_vector(A));
int lengthB = gmm::mat_nrows(gmm::col_vector(B));
GMM_ASSERT1(lengthA==lengthB, "impossible to scale the vectors: different lengths");

for(int i=0; i<lengthA; i++)
  B[i]= A[i]*B[i];
  
} /* end of vscale */



//Calculate the reciprocal of a vector, component by component.
//That means: B=A^^-1   --> B[i]= 1/A[i]

//As a matter of fact, it returns: A= A^^-1
template
<typename VEC>
void reciprocal (VEC & A){

//calculate the length of the vector

int lengthA = gmm::mat_nrows(gmm::col_vector(A));

for(int i=0; i<lengthA; i++)
  A[i]= 1/A[i];
  
} /* end of reciprocal */


//Calculate D=A*B*C.

template
<typename MAT1,typename MAT2,typename MAT3,typename MAT4>
void mult3 (const MAT1 & A,const MAT2 & B,const MAT3 & C, MAT4 & D){

//calculate the length of the vector

GMM_ASSERT1(mat_ncols(A)==mat_nrows(B), "first and second matrix dimensions don't match");
GMM_ASSERT1(mat_ncols(B)==mat_nrows(C), "second and third matrix dimensions don't match");
GMM_ASSERT1(mat_ncols(D)==mat_ncols(C), "final matrix has not the same number of columns of third matrix");
GMM_ASSERT1(mat_nrows(D)==mat_nrows(A), "final matrix has not the same number of rows of first matrix");

MAT4 AB(mat_nrows(A), mat_ncols(B)); clear(AB);
mult(A,B,AB); //AB = A*B
mult(AB,C,D); //D = AB*C = A*B*C
  
} /* end of mult3 */


} /* end of gmm namespace */


namespace getfem {

//Aux function for computing Peclet number
//Calculate the max value of a scalar or a vector field contained in a vector V.
// dim is te dimesion of te field, e.g. dim=1 is a scalar field and dim=3 is a vectorial 3-dimensional field
template
<typename VEC>
scalar_type max_vec (const VEC & V, size_type dim){


GMM_ASSERT1(dim>0, "wrong dimension of field: must be greater or equal to 1");

scalar_type max=0;
size_type size = gmm::vect_size(V);
if(dim>1){  // compute the norm of the vector field, then find the maximum
GMM_ASSERT1((size%dim)==0, "ERROR: the dimension of the field and the size of the vector do not match!");
VEC Vtemp(size/dim, 0);
for(int i=0; i<size/dim; i++){
   for( int j=0; j<dim; j++){
   Vtemp[i]+=(V[i+j])*(V[i+j]);
    }
   Vtemp[i]= sqrt(Vtemp[i]);
   if(Vtemp[i]>max) max= Vtemp[i]; //update maximum
  }
}
else if (dim==1){ //just find the maximum
for(int i=0; i<size; i++){
if(std::abs(V[i])>max) max= std::abs(V[i]);
}

} 
return max;
}  /* end of max_vec*/


//! Compute the Peclèt Number, defined as: 
//! @f$ \mathcal{P}e =  \frac{U~h}{A} @f$
//! that is the ratio between the advection flux and the diffusion coefficient
template
<typename VEC>
scalar_type peclet(const mesh & mesh,const VEC & U, const scalar_type & A, size_type  dim){


	#ifdef M3D1D_VERBOSE_
	cout <<"computing Peclet number...  "<<endl;
	#endif

// find Umax
scalar_type Umax= max_vec (U, dim);

// find h
scalar_type h=0;
scalar_type temp=0;
for(dal::bv_visitor i(mesh.convex_index()); !i.finished(); ++i){
	if(dim==1)
		temp= estimate_h(mesh, i);
	else if(dim==3)
		temp= 2*mesh.convex_radius_estimate(i);
	if(temp>h) h=temp;
	}

// compute peclet
scalar_type peclet= Umax*h/A;

	#ifdef M3D1D_VERBOSE_
	cout <<"U:   "<<Umax<<endl;
	cout <<"h:   "<<h<<endl;
	cout <<"A:   "<<A<<endl;
	cout <<"Peclet:   "<<peclet<<endl;
	#endif

return peclet;

} /* end of peclet_vessel */




} /* end of getfem namespace */






#endif
