/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   transport3d1d.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016.
  @brief  Definition of the aux class for the number of degrees of freedom.
 */
#ifndef M3D1D_DOF3D1D_TRANSP_HPP_
#define M3D1D_DOF3D1D_TRANSP_HPP_

namespace getfem {

//! Class to store the number of degrees of freedom of used FEMs
struct dof3d1d_transp {


	//! Number of dof of the interstitial concentration FEM mf_Pt
	size_type Ct_;

	//! Number of dof of the vessel concentration FEM mf_Pv
	size_type Cv_;

	//! Number of dof of the interstitial concentration FEM mf_Pt
	size_type Ct_Omega_;
	//! Number of dof of the interstitial concentration FEM mf_Pt
	size_type Ct_Sigma_;

	//! Total number of dof of the tissue problem
	size_type tissue_;
	//! Total number of dof of the vessel problem
	size_type vessel_;
	//! Total number of dof of the coupled problem
	size_type tot_;
	
	//! Compute the number of dof of given FEM
	void set (
			const getfem::mesh_fem & mf_Ct,	const getfem::mesh_fem & mf_Cv,
			const getfem::mesh_fem & mf_Ct_Omega, const getfem::mesh_fem & mf_Ct_Sigma
			)
	{
		Ct_ = mf_Ct.nb_dof();
		Cv_ = mf_Cv.nb_dof();

		Ct_Omega_ = mf_Ct_Omega.nb_dof();
		Ct_Sigma_ = mf_Ct_Sigma.nb_dof();
;
		
		tissue_ = Ct_; 
		vessel_ = Cv_ ;
		tot_    = tissue_+ vessel_;
	}

	void set (
			const getfem::mesh_fem & mf_Ct,	const getfem::mesh_fem & mf_Cv
			)
	{
		Ct_ = mf_Ct.nb_dof();
		Cv_ = mf_Cv.nb_dof();
		
		tissue_ = Ct_; 
		vessel_ = Cv_ ;
		tot_    = tissue_+ vessel_;
	}


	//! Accessor to the number of dof of mf_Ct
	inline size_type Ct (void) { return Ct_; } const

	//! Accessor to the number of dof of mf_Cv
	inline size_type Cv (void) { return Cv_; } const

	//! Accessor to the number of dof of mf_Ct_Omega
	inline size_type Ct_Omega (void) { return Ct_Omega_; } const
	//! Accessor to the number of dof of mf_Ct_Sigma
	inline size_type Ct_Sigma (void) { return Ct_Sigma_; } const

	//! Accessor to the number of dof of tissue problem
	inline size_type tissue (void) { return tissue_; } const
	//! Accessor to the number of dof of vessel problem
	inline size_type vessel (void) { return vessel_; } const
	//! Accessor to the number of dof of coupled problem
	inline size_type tot (void) { return tot_; } const

	//! Overloading of the output operator
	friend std::ostream & operator << (
			std::ostream & out, const dof3d1d_transp & dof
			)
	{ 
		out << "--- DEGREES OF FREEDOM --- " << endl;
		out << "  nb_dof_Ct     : " 		 << dof.Ct_  << endl;
		out << "  nb_dof_Cv     : " 		 << dof.Cv_  << endl;
		out << "  nb_dof_tot    : " 		 << dof.tot_ << endl;
		out << "-------------------------- " << endl;
		return out;            
	}

};

} /* end of namespace */

#endif
