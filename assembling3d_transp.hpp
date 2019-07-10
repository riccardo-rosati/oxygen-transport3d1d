 /* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   assembling3d_transp.hpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016 - May 2018.
  @brief  Miscelleanous assembly routines for the 3D transport problem.
 */


#ifndef M3D1D_ASSEMBLING_3D_TRANSP_HPP_
#define M3D1D_ASSEMBLING_3D_TRANSP_HPP_
#include <defines.hpp>
#include <utilities.hpp>

namespace getfem {

//! Build the mass, reaction and diffusion matrices for the 3D transport problem,
//! @f$ M = \int_{\Omega}   c~v~dx @f$ and
//! @f$ R = \int_{\Omega} r~c~v~dx @f$ and
//! @f$ D = \int_{\Omega}  d~\nabla c \cdot \nabla  v~dx @f$
/*! 
	@param M            Computed mass (time derivative) matrix
	@param D            Computed diffusion matrix
	@param R            Computed reaction matrix
	@param mim          The integration method to be used
	@param mf_c         The finite element method for the concentration @f$ c @f$
	@param mf_coef	    The finite element method for the coefficients
	@param reac_data    Coefficients for reaction term
	@param diff_data    Coefficients for diffusion term
	@param rg           The region where to integrate

	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void 
asm_tissue_transp
	(MAT & M, MAT & D,MAT & L, MAT & R,
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_coef,
	 const VEC & diff_data,
	 const VEC & linf_data,
	 const VEC & reac_data,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_c.get_qdim() == 1, 
		"invalid data mesh fem for pressure (Qdim=1 required)");
	// Build the mass matrix Mt (consumption)
	getfem::asm_mass_matrix_param(R, mim, mf_c, mf_coef, reac_data, rg);
	// Build the mass matrix Tt for time derivative 
	getfem::asm_mass_matrix(M, mim, mf_c, rg);
	// Build the divergence matrix Dtt
	getfem::asm_stiffness_matrix_for_laplacian(D,mim,mf_c, mf_coef, diff_data, rg); 
	// Build the mass matrix Lt (linfatic drainage)
	getfem::asm_mass_matrix_param(R, mim, mf_c, mf_coef, linf_data, rg);
	
} /* end of asm_tissue_transp*/


//! Build the advection matrice for the 3D transport problem,
//! @f$ B = \int_{\Omega} \mathbf{u} \cdot \nabla c~v~dx  ~+~  \int_{\Omega} \nabla \cdot \mathbf{u}   c~v~dx @f$
/*! 
	@param B            Advection matrix
	@param mim          The integration method to be used
	@param mf           The finite element method for the concentration @f$ c @f$
	@param mf vel       The finite element method for the advection field @f$ u @f$
	@param vel          The advection field
	@param rg           The region where to integrate

	@ingroup asm
 */  
  template<typename MAT, typename VECT>
  void asm_advection_tissue(MAT &B, const getfem::mesh_im &mim,
			    const getfem::mesh_fem &mf,
                            const getfem::mesh_fem &mfvel,
                            const VECT &vel,
                            const mesh_region & rg = mesh_region::all_convexes()                           
                            	 ) {
    getfem::generic_assembly
      assem1("vel=data(#2);"
            "M$1(#1,#1) += comp(Base(#1).Grad(#1).vBase(#2)) (:, :,i, k,i).vel(k);");
    assem1.push_mi(mim);
    assem1.push_mf(mf);
    assem1.push_mf(mfvel);
    assem1.push_data(vel);
    assem1.push_mat(B);
    assem1.assembly(rg);
    
    
    getfem::generic_assembly
      assem2("vel=data(#2);"
            "M$1(#1,#1) += comp( Base(#1).Base(#1).vGrad(#2) )(:, :,k, p,p).vel(k);");
    assem2.push_mi(mim);
    assem2.push_mf(mf);
    assem2.push_mf(mfvel);
    assem2.push_data(vel);
    assem2.push_mat(B);
    assem2.assembly(rg);
  }  /* end of asm_advection_tissue*/


/*! Build the Mixed boundary conditions (weak form) and Dirichlet (strong form) for vessels
    @f$ M=\int_{\Gamma_{MIX}} \beta~c~v~d\sigma@f$ and
    @f$ F=\int_{\Gamma_{MIX}} \beta~c_0~v~d\sigma@f$
 */
/*!
	@param F        BC contribution to rhs
	@param M        BC contribution to mass matrix
	@param mim      The integration method to be used
	@param mf_c     The finite element method for the concentration @f$ c @f$
	@param mf_data  The finite element method for the coefficients
	@param BC       Array of values of network boundary points
	@param beta     The beta value for mix condition
	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void
asm_tissue_bc_transp
	(VEC & F,
	 MAT & M,
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_data,
	 const std::vector<getfem::node> & BC,
	 const scalar_type beta
	 )
{

	
	
	GMM_ASSERT1(mf_c.get_qdim()==1,  "invalid data mesh fem (Qdim=1 required)");
	GMM_ASSERT1(mf_data.get_qdim()==1, "invalid data mesh fem (Qdim=1 required)");


	for (size_type bc=0; bc < BC.size(); ++bc) {
		GMM_ASSERT1(mf_c.linked_mesh().has_region(bc), "missed mesh region" << bc);
		if (BC[bc].label=="DIR") { // Dirichlet BC
			VEC BC_temp(mf_c.nb_dof(), BC[bc].value);
			getfem::assembling_Dirichlet_condition(M, F, mf_c, BC[bc].rg, BC_temp);
			gmm::clear(BC_temp);				
		} 
		else if (BC[bc].label=="MIX") { // Robin BC
			VEC BETA(mf_data.nb_dof(), beta);
			getfem::asm_mass_matrix_param(M, mim, mf_c, mf_data, BETA,mf_c.linked_mesh().region(BC[bc].rg) );
			
			VEC BETA_C0(mf_data.nb_dof(), beta*BC[bc].value);
			asm_source_term(F,mim, mf_c, mf_data,BETA_C0);
			
		}
		else if (BC[bc].label=="INT") { // Internal Node
			GMM_WARNING1("internal node passed as boundary.");
		}
		else if (BC[bc].label=="JUN") { // Junction Node
			GMM_WARNING1("junction node passed as boundary.");
		}
		else {
			GMM_ASSERT1(0, "Unknown Boundary Condition " << BC[bc].label << endl);
		}
	}

} /* end of asm_tissue_bc */



} /* end of namespace */

#endif
