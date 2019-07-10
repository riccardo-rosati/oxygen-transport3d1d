/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                         A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   assembling3d1d_transp.hpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016 -May 2018.
  @brief  Miscelleanous assembly routines for the 3D/1D coupling for transport problem.
 */


#ifndef M3D1D_ASSEMBLING_3D1D_TRANSP_HPP_
#define M3D1D_ASSEMBLING_3D1D_TRANSP_HPP_

#include <defines.hpp>
#include <utilities.hpp>
#include <utilities_transp.hpp>

namespace getfem {

/*!
	Build the averaging matrix @f$\bar{\Pi}_{tv}@f$ and the interpolation
	matrix @f${\Pi}_{tv}@f$. 
	Used to build the exchange matrices @f$B_{ii}@f$.
	@ingroup asm
 */
template<typename MAT, typename VEC>
void 
asm_exchange_aux_mat_transp
	(MAT & Mbar, MAT & Mlin,
	 const getfem::mesh_im & mim,
	 const getfem::mesh_fem & mf_t,
	 const getfem::mesh_fem & mf_v,
	 const getfem::mesh_fem & mf_coefv,
	 const VEC & RAD,
	 const size_type NInt,
	 const size_type nb_branches
	 ) 
{
	gmm::clear(Mbar); gmm::clear(Mlin);
	// Aux params
	const scalar_type Pi = 2*acos(0.0);
	size_type nb_dof_t     = mf_t.nb_dof();
	size_type nb_dof_v     = mf_v.nb_dof();

	// Linear interpolation map mf_t --> mf_v
	getfem::interpolation(mf_t, mf_v, Mlin);  
	
	// Interpolate RAD from mf_coefv to mf_v
	vector_type RADIUS(nb_dof_v);
	getfem::interpolation(mf_coefv, mf_v, RAD,RADIUS);  

	// Aux vectors for local interpolation
	std::vector<scalar_type> Pbari(NInt); 
	std::vector<scalar_type> Pt(nb_dof_t); 
	size_type counter = 0;
	size_type first_el_branch = 0; //counter for the first element of the branch
	for (size_type b = 0; b < nb_branches; b++){
		//mesh_region &rg_branch = mf_v.linked_mesh().region(b); // branch region
		dal::bit_vector dofs_b= mf_v.basic_dof_on_region(b);  // list of the indexes of all the dofs in the branch b
		first_el_branch = 0;
	for (dal::bv_visitor i(dofs_b); !i.finished(); ++i){
		counter++;
		first_el_branch++;
		if (counter*100 >= nb_dof_v) {
			counter = 0; 
			#ifdef M3D1D_VERBOSE_
			cout << "*"; cout.flush();
			#endif
		}
		// We need the following interpolation tool <getfem_interpolation.h>      
		getfem::mesh_trans_inv mti(mf_t.linked_mesh());
		// Build the list of point on the i-th circle:
		// ... first consider an orthonormal system v0, v1, v2:
		base_node v0;
		if (first_el_branch==1) {
			v0 = mf_v.point_of_basic_dof(i+1) - mf_v.point_of_basic_dof(i);			
		} else {
			v0 = mf_v.point_of_basic_dof(i) - mf_v.point_of_basic_dof(i-1);
		}
		base_node v1(0.0, -v0[2], v0[1]);
		base_node v2(v0[1]*v0[1] +v0[2]*v0[2], -v0[0]*v0[1], -v0[0]*v0[2]);
		if (gmm::vect_norm2(v2) < 1.0e-8 * gmm::vect_norm2(v0)) {
			v1[0] = -v0[1]; v1[1] = v0[0]; v1[2] = 0.0;
			v2[0] = -v0[0]*v0[2]; v2[1] = -v0[1]*v0[2]; v2[2] = v0[0]*v0[0] +v0[1]*v0[1];
		}
		v1 = v1 / gmm::vect_norm2(v1);
		v2 = v2 / gmm::vect_norm2(v2);
		// ... then parametrize the circle:
		for (size_type j = 0; j < NInt; j++){ 	
			mti.add_point( mf_v.point_of_basic_dof(i) + 
						   RADIUS[i]*(cos(2*Pi*j/NInt)*v1 + sin(2*Pi*j/NInt)*v2) );
		}
		// Get the local interpolation matrix Mbari
		MAT Mbari(NInt, nb_dof_t); gmm::clear(Mbari);
		interpolation(mf_t, mti, Pt, Pbari, Mbari, 1);
		scalar_type sum_row = 0.0;
		for (size_type j=0; j < NInt; ++j) {
			typename gmm::linalg_traits<MAT>::const_sub_row_type 
				row = mat_const_row(Mbari,j);
			gmm::linalg_traits< gmm::rsvector<scalar_type> >::const_iterator 
				it_nz = vect_const_begin(row);
			gmm::linalg_traits< gmm::rsvector<scalar_type> >::const_iterator 
				ite_nz = vect_const_end(row);
			for (; it_nz != ite_nz ; ++it_nz) {
				Mbar(i, it_nz.index()) += (*it_nz);
				sum_row += (*it_nz);
			}
		}
		typename gmm::linalg_traits<MAT>::sub_row_type 
			row = mat_row(Mbar,i);
		gmm::linalg_traits< gmm::rsvector<scalar_type> >::iterator
			it_nz = vect_begin(row);
		gmm::linalg_traits< gmm::rsvector<scalar_type> >::iterator
			ite_nz = vect_end(row);
		for (; it_nz != ite_nz ; ++it_nz) {
			(*it_nz)/= sum_row; //2*pi*RADIUS[i]:NO!!! //sum_row:notaro?  //Nint:pare quasi uguale a notaro
		}

	} /* end of convexes loop */

	} /* end of branch loop */




} /* end of build_aux_matrices */





//! Compute mass matrix with two parameters, that is:
//! @f$ M = \int_{\Omega} (~p_1~+~p_2 )~u~v~dx @f$ and
/*!
	@param M		Computed mass matrix
	@param mim		The integration metod used
	@param mf_c		The finite element method for @f$ u @f$ and @f$ v @f$
	@param mf_data1		The finite element method for parameter @f$ p_1 @f$
	@param mf_data2		The finite element method for parameter @f$ p_2 @f$
	@param PARAM1		First parameter @f$ p_1 @f$
	@param PARAM2		Second parameter @f$ p_2 @f$
	@param rg		The region where to integrate

*/
template<typename MAT, typename VEC1, typename VEC2>
void 
asm_mass_matrix_two_param
	(MAT & M, 
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_data1,
	 const mesh_fem & mf_data2,
	 const VEC1 & PARAM1,
	 const VEC2 & PARAM2,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{

	generic_assembly 
	assem("param1=data$1(#1); param2=data$2(#2); "
		"M$1(#3,#3)+=comp(Base(#3).Base(#3).Base(#1))(:,:,i).param1(i);" 			
		"M$1(#3,#3)+=comp(Base(#3).Base(#3).Base(#2))(:,:,i).param2(i);");
	assem.push_mi(mim);
	assem.push_mf(mf_data1);
	assem.push_mf(mf_data2);
	assem.push_mf(mf_c);
	assem.push_data(PARAM1);
	assem.push_data(PARAM2);
	assem.push_mat(M);
	assem.assembly(rg);
}



/*!
	Build the exchange matrices
	@f$B_{tt} = \left( - A + B \right)~ \Pi^T_{tv} M_{vv} \bar{\Pi}_{tv}@f$,    
	@f$B_{tv} = \left( - A - B \right)~ \Pi^T_{tv} M_{vv}@f$,    
	@f$B_{vt} = \left( + A - B \right)~ M_{vv} \bar{\Pi}_{tv}@f$,    
	@f$B_{vv} = \left( + A + B \right)~ M_{vv}@f$,    
	where @f$ A @f$ is the oncotic term and @f$ B @f$ is the vessel permeability term.   
	If ALT_FORM==true we substitute @f${\Pi}_{tv}@f$ with @f$\bar{\Pi}_{tv}@f$.
	@ingroup asm
 */
template<typename MAT, typename VEC>
void 
asm_exchange_mat_transp
	(MAT & Btt, MAT & Btv, MAT & Bvt, MAT & Bvv, 
	 const getfem::mesh_im & mim,
	 const getfem::mesh_fem & mf_c, 
	 const getfem::mesh_fem & mf_coefv,
	 const getfem::mesh_fem & mf_pv,
	 const MAT & Mbar, const MAT & Mlin,
	 const VEC & ONCOTIC,
	 const VEC & PERM,
	 const bool ALT_FORM
	 ) 
{


	MAT Bvv_temp(mf_c.nb_dof(),mf_c.nb_dof()); gmm::clear(Bvv_temp);

	#ifdef M3D1D_VERBOSE_
	cout << "    Assembling Bvv ..." << endl;  
	#endif
	getfem::asm_mass_matrix_two_param(Bvv, mim, 
					    mf_c, mf_pv, mf_coefv,
					gmm::scaled(ONCOTIC, +1),
					gmm::scaled(PERM, +1)); 
	#ifdef M3D1D_VERBOSE_
	cout << "    Assembling Bvt ..." << endl;
	#endif
	getfem::asm_mass_matrix_param(Bvv_temp, mim, 
					    mf_c, mf_pv, mf_coefv,
					gmm::scaled(PERM, -1)); 
	gmm::mult(Bvv_temp, Mbar, Bvt);
	gmm::clear(Bvv_temp);

	if (ALT_FORM){
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btv (alternative form) ..." << endl;
		#endif
		
		getfem::asm_mass_matrix_two_param(Bvv_temp, mim, 
					    mf_c, mf_pv, mf_coefv,
					gmm::scaled(ONCOTIC, -1),
					gmm::scaled(PERM, -1) ); 
		gmm::mult(gmm::transposed(Mbar), Bvv_temp, Btv); 
		gmm::clear(Bvv_temp);
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btt (alternative form) ..." << endl;
		#endif
		getfem::asm_mass_matrix_param(Bvv_temp, mim, 
					    mf_c, mf_pv, mf_coefv,
					gmm::scaled(PERM, +1) ); 

		gmm::mult3(gmm::transposed(Mbar),Bvv_temp,Mbar, Btt);
		gmm::clear(Bvv_temp);
	}
	else{
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btv (alternative form) ..." << endl;
		#endif
		
		getfem::asm_mass_matrix_two_param(Bvv_temp, mim, 
					    mf_c, mf_pv, mf_coefv,
					gmm::scaled(ONCOTIC, -1),
					gmm::scaled(PERM, -1) ); 
		gmm::mult(gmm::transposed(Mlin), Bvv_temp, Btv); 
		gmm::clear(Bvv_temp);
		#ifdef M3D1D_VERBOSE_
		cout << "    Assembling Btt (alternative form) ..." << endl;
		#endif
		getfem::asm_mass_matrix_param(Bvv_temp, mim, 
					    mf_c, mf_pv, mf_coefv,
					gmm::scaled(PERM, +1) ); 

		gmm::mult3(gmm::transposed(Mlin),Bvv_temp,Mbar, Btt);
		gmm::clear(Bvv_temp);
	}
	
} /* end of build_exchange_matrices */



  //! Build a single tissue Dirichlet condition on the network (modify @f$ B_{vt} @f$)
  /* @f{equation}{
	B = 
	@f{bmatrix} {B_{tt}& B_{tv}\\
		     B_{vt}& B_{vv}  @f}
	}
  */
  //! @f$ BU = F @f$

  /*!

	@param B	The monolitic matrix of the system
	@param F	The right hand side vector of the system  
	@param mf1 	The finite element method on tissue
	@param mf2 	The finite element method on network
	@param boundary	The index of the boundary region to be buildt
	@param DIR	The vector containing the Dirichlet condition (should be of the same dimension of F)

     @ingroup asm
  */
 template<typename MATRM, typename VECT1, typename VECT2>
  void assembling_Dirichlet_condition_coupled_tissue
  (MATRM &B, VECT1 &F, const mesh_fem &mf1, const mesh_fem &mf2, size_type boundary,
   const VECT2 &DIR) {
   

    size_type Q1=mf1.get_qdim();
    size_type Q2=mf2.get_qdim();

    size_type nb_dof1=mf1.nb_dof();
    size_type nb_dof2=mf2.nb_dof();    
    
    GMM_ASSERT1(!(mf1.is_reduced()), "This function is not adapted to "
		"reduced finite element methods"); 
    GMM_ASSERT1(!(mf2.is_reduced()), "This function is not adapted to "
		"reduced finite element methods"); 
						
    dal::bit_vector nndof = mf1.basic_dof_on_region(boundary);
    pfem pf1;
    
    for (dal::bv_visitor cv(mf1.convex_index()); !cv.finished(); ++cv) {	 	//per tutti i convessi cv della mesh 1
    
      pf1 = mf1.fem_of_element(cv);
      pdof_description ldof = lagrange_dof(pf1->dim());
      size_type nbd = pf1->nb_dof(cv);	      
      for (size_type i = 0; i < nbd; i++) {					//per tutti i dof i del convesso cv
	size_type dof1 = mf1.ind_basic_dof_of_element(cv)[i*Q1];				//trova l'indice delle colonne riferite all
	if (nndof.is_in(dof1) && pf1->dof_types()[i] == ldof) {			//se il dof i del convesso cv è in "boundary"
  
	  for (size_type j = nb_dof1; j < nb_dof1+ nb_dof2; j++) {				//allora per tutti i dof j della mesh 2
		for (size_type l = 0; l < Q1; ++l) {
			F[j] -= B(j, dof1+l) * DIR[dof1+l];
	    		B(j, dof1+l) =  0;
	    		B(dof1+l, j) =  0;
	    	}
	    }
	  } 
	}
     }
   } /* end of assembling_Dirichlet_condition_coupled_tissue*/

  //! Build a single vessel Dirichlet condition on the tissue (modify @f$ B_{tv} @f$)
  /* @f{equation}{
	B = 
	@f{bmatrix} {B_{tt} & B_{tv}\\ B_{vt} & B_{vv}  @f}
	}
  */
  //! @f$ BU = F @f$

  /*!

	@param B	The monolitic matrix of the system
	@param F	The right hand side vector of the system  
	@param mf1 	The finite element method on tissue
	@param mf2 	The finite element method on network
	@param boundary	The index of the boundary region to be buildt
	@param DIR	The vector containing the Dirichlet condition (should be of the same dimension of F)

     @ingroup asm
  */
 template<typename MATRM, typename VECT1, typename VECT2>
  void assembling_Dirichlet_condition_coupled_vessel
  (MATRM &B, VECT1 &F, const mesh_fem &mf1, const mesh_fem &mf2, size_type boundary,
   const VECT2 &DIR) {
   

    size_type Q1=mf1.get_qdim();
    size_type Q2=mf2.get_qdim();

    size_type nb_dof1=mf1.nb_dof();
    size_type nb_dof2=mf2.nb_dof();    
    
    GMM_ASSERT1(!(mf1.is_reduced()), "This function is not adapted to "
		"reduced finite element methods"); 
    GMM_ASSERT1(!(mf2.is_reduced()), "This function is not adapted to "
		"reduced finite element methods"); 
					
    dal::bit_vector nndof = mf2.basic_dof_on_region(boundary);
    pfem pf2;
    
    for (dal::bv_visitor cv(mf2.convex_index()); !cv.finished(); ++cv) {	 	//per tutti i convessi cv della mesh 1
    
      pf2 = mf2.fem_of_element(cv);
      pdof_description ldof = lagrange_dof(pf2->dim());
      size_type nbd = pf2->nb_dof(cv);	      
      for (size_type i = 0; i < nbd; i++) {					//per tutti i dof i del convesso cv
	size_type dof2 = mf2.ind_basic_dof_of_element(cv)[i*Q2];				//trova l'indice delle colonne riferite all
	if (nndof.is_in(dof2) && pf2->dof_types()[i] == ldof) {			//se il dof i del convesso cv è in "boundary"
	  
	  for (size_type j = 0; j < nb_dof1; j++) {				//allora per tutti i dof j della mesh 2
		for (size_type l = 0; l < Q2; ++l) {
			F[j] -= B(j, nb_dof1 + dof2+l) * DIR[nb_dof1 + dof2+l];
	    		B(j, nb_dof1 + dof2+l) =  0;
	    		B(nb_dof1 + dof2+l, j) =  0;
	    	}
	    }
	  } 
	}
     }
   } /* end of assembling_Dirichlet_condition_coupled_vessel*/
   

  //! Build all the Dirichlet conditions on the coupling matrixes (modify both @f$ B_{tv} @f$ and @f$ B_{vt} @f$))
  /* @f{equation}{
	B = 
	@f{bmatrix} {B_{tt} & B_{tv}\\ B_{vt} & B_{vv}  @f}
	}
  */
  //! @f$ BU = F @f$

  /*!
	@param M		The monolitic matrix of the system
	@param F		The right hand side vector of the system  
	@param mf_ct 		The finite element method on tissue
	@param mf_cv 		The finite element method on network
	@param BC_tissue	List of nodes of the boundary regions in tissue
	@param BC_vessel	List of nodes of the boundary regions in network

     @ingroup asm
  */
template<typename MAT, typename VEC>
void
asm_coupled_bc_transp
	(MAT & M,
 	 VEC & F,
	 const mesh_fem & mf_ct,
	 const mesh_fem & mf_cv,
	 const std::vector<getfem::node> & BC_tissue,
	 const std::vector<getfem::node> & BC_vessel
	 )
{

	
	
	GMM_ASSERT1(mf_ct.get_qdim()==1,  "invalid data mesh fem (Qdim=1 required)");
	GMM_ASSERT1(mf_cv.get_qdim()==1,  "invalid data mesh fem (Qdim=1 required)");


	//cycle over the tissue boundary nodes
	for (size_type bc=0; bc < BC_tissue.size(); ++bc) {
		GMM_ASSERT1(mf_ct.linked_mesh().has_region(bc), "missed mesh region" << bc);
		if (BC_tissue[bc].label=="DIR") { // Dirichlet BC
			VEC BC_temp(mf_ct.nb_dof(), BC_tissue[bc].value);
			getfem::assembling_Dirichlet_condition_coupled_tissue(M, F, mf_ct, mf_cv, BC_tissue[bc].rg, BC_temp);
			gmm::clear(BC_temp);				
		} 
	}
	
	//cycle over the vessels boundary nodes
	for (size_type bc=0; bc < BC_vessel.size(); ++bc) {
		GMM_ASSERT1(mf_cv.linked_mesh().has_region(bc), "missed mesh region" << bc);
		if (BC_vessel[bc].label=="DIR") { // Dirichlet BC
			VEC BC_temp(mf_cv.nb_dof(), BC_vessel[bc].value);
			getfem::assembling_Dirichlet_condition_coupled_vessel(M, F, mf_ct, mf_cv, BC_vessel[bc].rg, BC_temp);
			gmm::clear(BC_temp);				
		} 
	}


} /* end of asm_coupled_bc_transp */


/*!
	Build the averaging matrix @f$\bar{\bar{\Pi}}_{tv}@f$.
	@ingroup asm
 */
template<typename MAT, typename VEC>
void 
asm_exchange_aux_mat_bar_bar
	(MAT & Mbarbar,
	 const getfem::mesh_im & mim,
	 const getfem::mesh_fem & mf_t,
	 const getfem::mesh_fem & mf_v,
	 const getfem::mesh_fem & mf_coefv,
	 const VEC & RAD,
	 const size_type NInt,
	 const size_type NIntA,
	 const size_type nb_branches
	 ) 
{
	gmm::clear(Mbarbar);
	// Aux params
	const scalar_type Pi = 2*acos(0.0);
	size_type nb_dof_t     = mf_t.nb_dof();
	size_type nb_dof_v     = mf_v.nb_dof();

	size_type Ntot= NInt*NIntA+1;
	// Interpolate RAD from mf_coefv to mf_v
	vector_type RADIUS(nb_dof_v);
	getfem::interpolation(mf_coefv, mf_v, RAD,RADIUS);  

	// Aux vectors for local interpolation
	std::vector<scalar_type> Pbari(Ntot); 
	std::vector<scalar_type> Pt(nb_dof_t); 
	size_type counter = 0;
	size_type first_el_branch = 0; //counter for the first element of the branch
	for (size_type b = 0; b < nb_branches; b++){
		//mesh_region &rg_branch = mf_v.linked_mesh().region(b); // branch region
		dal::bit_vector dofs_b= mf_v.basic_dof_on_region(b);  // list of the indexes of all the dofs in the branch b
		first_el_branch = 0;
	for (dal::bv_visitor i(dofs_b); !i.finished(); ++i){
		counter++;
		first_el_branch++;
		if (counter*100 >= nb_dof_v) {
			counter = 0; 
			#ifdef M3D1D_VERBOSE_
			cout << "*"; cout.flush();
			#endif
		}
		// We need the following interpolation tool <getfem_interpolation.h>      
		getfem::mesh_trans_inv mti(mf_t.linked_mesh());
		// Build the list of point on the i-th circle:
		// ... first consider an orthonormal system v0, v1, v2:
		base_node v0;
		if (first_el_branch==1) {
			v0 = mf_v.point_of_basic_dof(i+1) - mf_v.point_of_basic_dof(i);			
		} else {
			v0 = mf_v.point_of_basic_dof(i) - mf_v.point_of_basic_dof(i-1);
		}
		base_node v1(0.0, -v0[2], v0[1]);
		base_node v2(v0[1]*v0[1] +v0[2]*v0[2], -v0[0]*v0[1], -v0[0]*v0[2]);
		if (gmm::vect_norm2(v2) < 1.0e-8 * gmm::vect_norm2(v0)) {
			v1[0] = -v0[1]; v1[1] = v0[0]; v1[2] = 0.0;
			v2[0] = -v0[0]*v0[2]; v2[1] = -v0[1]*v0[2]; v2[2] = v0[0]*v0[0] +v0[1]*v0[1];
		}
		v1 = v1 / gmm::vect_norm2(v1);
		v2 = v2 / gmm::vect_norm2(v2);
		// ... then parametrize the circle:
		mti.add_point( mf_v.point_of_basic_dof(i));
		for (size_type j = 0; j < NInt; j++){ 
			for (size_type k = 1; k <= NIntA; k++){ 
				scalar_type partial_Radius = RADIUS[i]*k/NIntA;
				scalar_type partial_Angle = 2*Pi*j/NInt;				
	//if (j==3 && k==4) cout <<"j==3, k==4 : partial_radius = "<<partial_Radius <<"; partial_Angle = "<<partial_Angle<<endl;
				mti.add_point( mf_v.point_of_basic_dof(i) + 
						   partial_Radius*(cos(partial_Angle)*v1 + sin(partial_Angle)*v2) );
			}
		}
		
		// Get the local interpolation matrix Mbari
		MAT Mbari(Ntot, nb_dof_t); gmm::clear(Mbari);
		interpolation(mf_t, mti, Pt, Pbari, Mbari, 1);
		scalar_type sum_row = 0.0;
		for (size_type j=0; j < Ntot; ++j) {
			typename gmm::linalg_traits<MAT>::const_sub_row_type 
				row = mat_const_row(Mbari,j);
			gmm::linalg_traits< gmm::rsvector<scalar_type> >::const_iterator 
				it_nz = vect_const_begin(row);
			gmm::linalg_traits< gmm::rsvector<scalar_type> >::const_iterator 
				ite_nz = vect_const_end(row);
			for (; it_nz != ite_nz ; ++it_nz) {
				Mbarbar(i, it_nz.index()) += (*it_nz);
				//cout <<"riga i="<<i<<",colonna j="<<j<<"; it_nz = "<<(*it_nz)<<endl;
				sum_row += (*it_nz);
			}
		}
		//cout <<"riga i = "<< i <<": sum_row = "<<sum_row<<endl;
		typename gmm::linalg_traits<MAT>::sub_row_type 
			row = mat_row(Mbarbar,i);
		gmm::linalg_traits< gmm::rsvector<scalar_type> >::iterator
			it_nz = vect_begin(row);
		gmm::linalg_traits< gmm::rsvector<scalar_type> >::iterator
			ite_nz = vect_end(row);
		for (; it_nz != ite_nz ; ++it_nz) {
			(*it_nz)/= sum_row; 
		}

	} /* end of convexes loop */

	} /* end of branch loop */




} /* end of build_aux_matrices */



} /* end of namespace */

#endif
