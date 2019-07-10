 /* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   assembling1d_transp.hpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016 - May 2018
  @brief  Miscelleanous assembly routines for the 1D transport problem.
 */
/** @defgroup asm Assembly routines */

#ifndef M3D1D_ASSEMBLING_1D_TRANSP_HPP_
#define M3D1D_ASSEMBLING_1D_TRANSP_HPP_

#include <defines.hpp>
#include <node.hpp>
#include <utilities.hpp>
#include <utilities_transp.hpp>

namespace getfem {

//! Build the mass and diffusion matrices for the 1D Poiseuille's problem,
//! @f$ M = \int_{\Lambda} \pi~R^2~c~v~ds @f$ and
//! @f$ D = \int_{\Lambda} \pi~R^2~A~\frac{\partial~c}{\partial~s} \, \frac{\partial~v}{\partial~s} ~ds @f$
/*!
	@param M         Computed mass (time derivative) matrix
	@param D         Computed diffusion matrix
	@param mim       The integration method to be used
	@param mf_c      The finite element method for the concentration @f$ c @f$
	@param mf_data   The finite element method for the diffusion coefficient
	@param diff      The diffusion coefficient for D
	@param R	 Network radii
	@param rg        The region where to integrate

	@ingroup asm
 */ 
template<typename MAT, typename VEC, typename VEC2>
void 
asm_network_transp
	(MAT & M, MAT & D, 
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_data,
	 const VEC & diff,
	 const VEC2 & R,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf_c.get_qdim() == 1 ,
		"invalid data mesh fem (Qdim=1 required)");
	//build mass matrix Mv for time derivative
	VEC param(mf_data.nb_dof()); gmm::clear(param);
	gmm::add(R, param);
	gmm::vscale(R, param);
	gmm::scale(param, pi); //param = pi*R^2
	getfem::asm_mass_matrix_param(M, mim, mf_c, mf_data, param, rg);
	// Build the diffusion matrix Dv
	gmm::vscale(diff, param); //param= pi*R^2*Av
	getfem::asm_stiffness_matrix_for_laplacian(D,mim,mf_c,mf_data, param, rg);

} //end of asm_network_transp


//! Build the advection matrice for the 1D Poiseuille's problem,
//! @f$ B = \int_{\Lambda} ~u  ~\nabla c \cdot \mathbf{\lambda} \,v ~ds ~ + ~  \int_{\Lambda} ~c  ~\nabla\cdot \left( u ~\mathbf{\lambda}\right) \,v~ds @f$
/*!
	@param B         Computed advection matrix
	@param mim       The integration method to be used
	@param mf_c      The finite element method for the concentration @f$ c @f$
	@param mf_data   The finite element method for the the tangent versor on @f$ \Lambda @f$
	@param mf_u      The finite element method for the advection field  @f$ u @f$
	@param mf_R	 The finite element method for the network radii
	@param U         The advection field
	@param lambdax   First cartesian component of the tangent versor  @f$ \mathbf{\lambda} @f$
	@param lambday   Second cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param lambdaz   Third cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param R	 Network radii	
	@param rg        The region where to integrate

	@ingroup asm
 */ 
 
template<typename MAT, typename VEC>
void 
asm_advection_network
	(MAT & B,
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_data,
	 const mesh_fem & mf_u,
	 const mesh_fem & mf_R,
	 const VEC & U,
	 const VEC & lambdax, const VEC & lambday, const VEC & lambdaz,
	 const VEC & R,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 	
	 
	 {

	generic_assembly 
	assem1("l1=data$1(#2); l2=data$2(#2); l3=data$3(#2);  u=data$4(#3); R=data$5(#4);"
		  "t=comp(Base(#1).Grad(#1).Base(#2).Base(#3).Base(#4).Base(#4));"
		  "M$1(#1,#1)+=t(:,:,1,i,p,m,n).l1(i).u(p).R(m).R(n)+t(:,:,2,i,p,m,n).l2(i).u(p).R(m).R(n)+t(:,:,3,i,p,m,n).l3(i).u(p).R(m).R(n);");
	assem1.push_mi(mim);
	assem1.push_mf(mf_c);
	assem1.push_mf(mf_data);
	assem1.push_mf(mf_u);
	assem1.push_mf(mf_R);
	assem1.push_data(lambdax);
	assem1.push_data(lambday);
	assem1.push_data(lambdaz);
	assem1.push_data(U);
	assem1.push_data(R);
	assem1.push_mat(B);
	assem1.assembly(rg);
	
		
	
	generic_assembly 
	assem2("l1=data$1(#2); l2=data$2(#2); l3=data$3(#2);  u=data$4(#3);  R=data$5(#4);"
		  "t=comp(Base(#1).Base(#1).Base(#2).Grad(#3).Base(#4).Base(#4));"
		  "M$1(#1,#1)+=t(:,:,i,p,1,m,n).l1(i).u(p).R(m).R(n)+t(:,:,i,p,2,m,n).l2(i).u(p).R(m).R(n)+t(:,:,i,p,3,m,n).l3(i).u(p).R(m).R(n);"); 
	assem2.push_mi(mim);
	assem2.push_mf(mf_c);
	assem2.push_mf(mf_data);
	assem2.push_mf(mf_u);
	assem2.push_mf(mf_R);
	assem2.push_data(lambdax);
	assem2.push_data(lambday);
	assem2.push_data(lambdaz);
	assem2.push_data(U);
	assem2.push_data(R);
	assem2.push_mat(B);
	assem2.assembly(rg);
} //end of asm_advection_network

//RR: sto assemblando i due vettori ottenuti da pi*R*R*d(uv*psi)/dS

template<typename MAT, typename VEC>
void 
asm_hemoadvection_rhs_network
	(VEC & Ov
			const mesh_im & mim,
			const mesh_fem & mf_c,
			const mesh_fem & mf_data,
			const mesh_fem & mf_u,
			const mesh_fem & mf_R,
			const mesh_fem & mf_H,
			const VEC & U,
			const VEC & lambdax, const VEC & lambday, const VEC & lambdaz,
			const VEC & R,
			const VEC & psi,
			const mesh_region & rg = mesh_region::all_convexes())	
	 
	{
generic_assembly
assem1("l1=data$1(#2); l2= data$2(#2); l3=data$3(#2); u=data$4(#3); R=data$5(#4); psi=data$6(#5);"
	"t=comp(Base(#1).Base(#2).Base(#3).Grad(#5).Base(#4).Base(#4));"
	"V$1(#1)+=t(:,i,p,d,1,m,n).l1(i).u(p).psi(d).R(m).R(n)+t(:,i,p,d,2,m,n).l2(i).u(p).psi(d).R(m).R(n)+t(:,i,p,d,3,m,n).l3(i).u(p).psi(d).R(m).R(n);");
assem1.push_mi(mim);
assem1.push_mf(mf_c);
assem1.push_mf(mf_data);
assem1.push_mf(mf_u);
assem1.push_mf(mf_R);
assem1.push_mf(mf_H);
assem1.push_data(lambdax);
assem1.push_data(lambday);
assem1.push_data(lambdaz);
assem1.push_data(U);
assem1.push_data(R);
assem1.push_data(psi);
assem1.push_vec(Ov);
assem1.assembly(rg);

generic_assembly
assem2("l1=data$1(#2); l2= data$2(#2); l3=data$3(#2); u=data$4(#3); R=data$5(#4); psi=data$6(#5);"
	"t=comp(Base(#1).Base(#2).Grad(#3).Base(#5).Base(#4).Base(#4));"
	"V$1(#1)+=t(:,i,p,1,d,m,n).l1(i).u(p).psi(d).R(m).R(n)+t(:,i,p,2,d,m,n).l2(i).u(p).psi(d).R(m).R(n)+t(:,i,p,3,d,m,n).l3(i).u(p).psi(d).R(m).R(n);");
assem2.push_mi(mim);
assem2.push_mf(mf_c);
assem2.push_mf(mf_data);
assem2.push_mf(mf_u);
assem2.push_mf(mf_R);
assem2.push_mf(mf_H);
assem2.push_data(lambdax);
assem2.push_data(lambday);
assem2.push_data(lambdaz);
assem2.push_data(U);
assem2.push_data(R);
assem2.push_data(psi);
assem2.push_vec(Ov);
assem2.assembly(rg);

} //end of asm_hemoadvection_rhs_network


/*! Build the mixed boundary conditions (weak form) and dirichlet (strong form) for vessels
    @f$ M=\int_{\mathcal{E}_{MIX}} \beta~c~v~d\sigma@f$ and
    @f$ F=\int_{\mathcal{E}_{MIX}} \beta~c_0~v~d\sigma@f$
 */
/*!
	@param F        BC contribution to rhs
	@param M        BC contribution to mass matrix
	@param mim      The integration method to be used
	@param mf_c     The finite element method for the concentration @f$ c @f$
	@param mf_data  The finite element method for the coefficients
	@param BC       Array of values of network boundary points
	@param beta     The beta value for mix condition @f$ p_0 @f$
	@param R	Network radii
	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void
asm_network_bc_transp
	(VEC & F, MAT & M, 
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_fem & mf_data,
	 const std::vector<getfem::node> & BC,
	 const scalar_type beta,
	 const VEC & R) 
{ 
	GMM_ASSERT1(mf_c.get_qdim()==1,  "invalid data mesh fem (Qdim=1 required)");
	GMM_ASSERT1(mf_data.get_qdim()==1, "invalid data mesh fem (Qdim=1 required)");


	for (size_type bc=0; bc < BC.size(); bc++) { 
		GMM_ASSERT1(mf_c.linked_mesh().has_region(bc), "missed mesh region" << bc);
		if (BC[bc].label=="DIR") { // Dirichlet BC
			VEC BC_temp(mf_c.nb_dof(), BC[bc].value);
			getfem::assembling_Dirichlet_condition(M, F, mf_c, BC[bc].rg, BC_temp);
			gmm::clear(BC_temp);			
		} 
		else if (BC[bc].label=="MIX") { // Robin BC
			VEC BETA(mf_data.nb_dof(), beta*pi);
			gmm::vscale(R, BETA); gmm::vscale(R, BETA);
			getfem::asm_mass_matrix_param(M, mim, mf_c, mf_data, BETA,mf_c.linked_mesh().region(BC[bc].rg) ); //int(beta*cv*bv)
			
			VEC BETA_C0(mf_data.nb_dof(), pi*beta*BC[bc].value);
			gmm::vscale(R, BETA_C0); gmm::vscale(R, BETA_C0);
			asm_source_term(F,mim, mf_c, mf_data,BETA_C0); //int(beta*c0*bv)
		}
		else if (BC[bc].label=="INT") { // Internal Node
			GMM_WARNING1("internal node passed as boundary.");
		}
		else if (BC[bc].label=="JUN") { // Junction Node
			GMM_WARNING1("junction node passed as boundary.");
		}
		else {
			GMM_ASSERT1(0, "Unknown Boundary Condition"<< BC[bc].label << endl);
		}
	}

}


} /* end of namespace */

#endif
