/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   param3d1d_transp.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016.
  @brief  Definition of the aux class for physical parameters.
 */
 
#ifndef M3D1D_PARAM3D1D_TRANSP_HPP_
#define M3D1D_PARAM3D1D_TRANSP_HPP_

#include <mesh1d.hpp>    // import_network_radius
#include <utilities.hpp> // compute_radius

namespace getfem {

//! Class to handle the physical parameter of the coupled 3D/1D model
/*!
	\ingroup input
 */
struct param3d1d_transp {

	// Dimensional physical parameters (microcirc applications)
	//Diffusivity in the tissue [m^2/s]
	scalar_type Dt_;
	//Diffusivity in the vessels [m^2/s]
	scalar_type Dv_;
	//rate of metabolization [1/s]
	scalar_type m_;
	//Permeability of the vessel wall [m/s]
	scalar_type Perm_;
	//hydraulic conductivity of the lymphatic wall [s * m^2/kg]
	scalar_type Lp_LF_;
	// surface area of lymphatic vessels per unit volume of tissue [1/m]
	scalar_type SV_;

//MODIFICHE	
	//Max rate of oxygen metabolization [1/s]
	scalar_type m0_;
	//Tissue Concentration guess [kg/m^3]
	scalar_type Ct0_;
	//Vessel Concentration guess [kg/m^3]
	scalar_type Cv0_;
	//Partial pressure at half max rate of metabolization [mmHg]
	scalar_type Pm_50_;
	//Solubility of oxygen in the tissue [kg/(m^3*mmHg]
	scalar_type alpha_t_;
	//Hill constant [-]
	scalar_type delta_;
	//Mean Corpuscolar Hemoglobin Concentration [-]
	scalar_type MCHC_;
	//Hufner factor [-];
	scalar_type N_;
	//Oxygen tissue solubility [kg/(m^3*mmHg);
	scalar_type alpha_pl_;
	//Partial Pressure of oxygen at half saturation [mmHg]
	scalar_type Ps_50_; 

	
	
	// Dimensionless physical parameters (test-cases)
	//Inverse of Peclet number for tissue
	vector_type At_;
	//Inverse of Peclet number for vessel
	vector_type Av_;
	//Damkohler number (metabolism VS diffusion)
	vector_type Dalpha_;
	//Magnitude of leakage from the capillary bed
	vector_type Y_;
	//lymphatic drainage
	vector_type Q_pl_;
	
	//time parameters 
	// simulation time length
	scalar_type T_;
	// time step
	scalar_type dt_;
	// initial concentration in tissue
	scalar_type C0t_;
	// initial concentration in network
	scalar_type C0v_;
	
	

	// Utils
	//! File .param
	ftool::md_param FILE_;
	//! Finite Element Method for tissue data
	getfem::mesh_fem mf_datat_;
	//! Finite Element Method for vessel data
	getfem::mesh_fem mf_datav_;
	// Methods
	//! Build the arrays of dimensionless parameters
	void build(ftool::md_param & fname,
			const getfem::mesh_fem & mf_datat,
			 const getfem::mesh_fem & mf_datav
			) 
	{
		FILE_ = fname;
		mf_datat_ = mf_datat;
		mf_datav_ = mf_datav;
		size_type dof_datat = mf_datat_.nb_dof();
		size_type dof_datav = mf_datav_.nb_dof();
		 
//		bool IMPORT_RADIUS = FILE_.int_value("IMPORT_RADIUS");
		bool NONDIM_PARAM  = FILE_.int_value("TEST_PARAM");
		bool EXPORT_PARAM  = FILE_.int_value("EXPORT_PARAM");
		
		
		#ifdef M3D1D_VERBOSE_
		cout << "  Assembling dimensionless parameters Dt, Dv, Dalpha, Q_pl ... "   << endl;
		#endif
		if (NONDIM_PARAM) {
			// Import dimensionless params from FILE_
			scalar_type Atval = FILE_.real_value("At"); 
			scalar_type Avval = FILE_.real_value("Av"); 
			scalar_type Dalphaval = FILE_.real_value("D_alpha"); 
			scalar_type Yval = FILE_.real_value("Y"); 
			scalar_type Q_plval  = FILE_.real_value("Q_pl"); 
			// Fill the data arrays
			 At_.assign(dof_datat,  Atval);
			 Av_.assign(dof_datav,  Avval);
			 Dalpha_.assign(dof_datat,  Dalphaval);
			 Y_.assign(dof_datav,  Yval);
			 Q_pl_.assign(dof_datat,  Q_plval);

			T_   = FILE_.real_value("T","Simulation time length [s]");
			dt_   = FILE_.real_value("dt","Time step [s]");			
			C0t_   = FILE_.real_value("C0t","Initial concentration in tissue []");	
			C0v_   = FILE_.real_value("C0v","Initial concentration in network []");			
		} 
		else { 
			// Import dimensional params from FILE_
			scalar_type P_  = FILE_.real_value("P", "average interstitial pressure [Pa]"); 
			scalar_type U_  = FILE_.real_value("U", "characteristic flow speed in the capillary bed [m/s]"); 
			scalar_type d_  = FILE_.real_value("d", "characteristic length of the problem [m]"); 
			scalar_type k_  = FILE_.real_value("k", "permeability of the interstitium [m^2]"); 
			scalar_type Lp_ = FILE_.real_value("Lp", "Hydraulic conductivity of the capillary walls [m^2 s/kg]"); 
			
			Dt_   = FILE_.real_value("Dt","Diffusivity in the tissue [m^2/s]");
			Dv_   = FILE_.real_value("Dv","Diffusivity in the vessels [m^2/s]");

		//RR: MICHEALIS-MENTEN: modifica termine di reazione: da m_ a m0_ (consumo massimo):
			m0_    = FILE_.real_value("m0","Max rate of oxygen metabolization [1/s]");
			Ct0_	= FILE_.real_value("Ct0","Tissue Concentration guess [kg/m^3]");
			Cv0_	= FILE_.real_value("Cv0","Vessel Concentration guess [kg/m^3]");
			Pm_50_	= FILE_.real_value("Pm_50","Partial pressure at half max rate of metabolization [mmHg]");
			alpha_t_ = FILE_.real_value("alpha_t","Solubility of oxygen in the tissue [kg/(m^3*mmHg]");

	
		//RR: Parametri OSSIEMOGLOBINA
			MCHC_	= FILE_.real_value("MCHC","Mean Corpuscolar Hematocrit Concentration [-]");
			N_ = FILE_.real_value("N","Hufner factor [-]");
			delta_ = FILE_.real_value("delta","Hill constant [-]");
			Ps_50_ = FILE_.real_value("Ps_50","Partial pressure at half saturation [mmHg]");
			alpha_pl_ = FILE_.real_value("alpha_pl","oxygen solubility in the plasma [kg/(m^3*mmHg)]");


			Perm_ = FILE_.real_value("Perm","Permeability of the capillary walls [m/s]");
			Lp_LF_ = FILE_.real_value("Lp_LF","hydraulic conductivity of the lymphatic wall [s * m^2/kg]");
			SV_ = FILE_.real_value("SV","surface area of lymphatic vessels per unit volume of tissue [1/m]");

			T_   = FILE_.real_value("T","Simulation time length [s]");
			dt_   = FILE_.real_value("dt","Time step [s]");	
			C0t_   = FILE_.real_value("C0t","Initial concentration in tissue []");	
			C0v_   = FILE_.real_value("C0v","Initial concentration in network []");				
			// Compute the dimentionless params
			At_.assign(dof_datat, Dt_/d_/U_);
			Av_.assign(dof_datav, Dv_/d_/U_);

		//modifica coefficiente di reazione
			//Dalpha_.assign(dof_datat, (m0_/((C0_+(Pm_50_*alpha_T_)))*U_*d_);
			Y_.assign(dof_datav, Perm_/U_);
			Q_pl_.assign(dof_datat,Lp_LF_*SV_*P_*d_/U_);
						
		}
		// Check values
		GMM_ASSERT1(At_[0] != 0, "wrong tissue diffusivity (At>0 required)"); 
		GMM_ASSERT1(Av_[0] != 0, "wrong vessel bed diffusivity (Av>0 required)");
		//if (Q_[0] == 0) cout << "Warning: uncoupled problem (Q=0)" << endl;
		
		if (EXPORT_PARAM){
/*			std::string ODIR = FILE_.string_value("OutputDir","OutputDirectory");
			getfem::vtk_export exp(ODIR+"radius.vtk");
			exp.exporting(mf_datav_);
			exp.write_mesh();
			exp.write_point_data(mf_datav_, R_, "R");
			getfem::vtk_export expQ(ODIR+"conductivity.vtk");
			expQ.exporting(mf_datav_);
			expQ.write_mesh();
			expQ.write_point_data(mf_datav_, Q_, "Q"); */
		}

	}
	//! Get the radius at a given dof
	//inline scalar_type R  (size_type i) { return R_[i];  } const
	//! Get the tissue diffusivity at a given dof
	inline scalar_type At (size_type i) { return At_[i]; } const
	//! Get the vessel diffusivity at a given dof
	inline scalar_type Av (size_type i) { return Av_[i]; } const
	//! Get the linphatic drainage at a given dof
	inline scalar_type Q_pl  (size_type i) { return Q_pl_[i];  } const
	//! Get the Dahmkholer number at a given dof
	inline scalar_type Dalpha  (size_type i) { return Dalpha_[i];  } const
	//! Get the leakage of the capillary bed at a given dof
	inline scalar_type Y  (size_type i) { return Y_[i];  } const
	//! Get the simulation time length
	inline scalar_type T  () { return T_;  } const
	//! Get the time step
	inline scalar_type dt  () { return dt_;  } const
	//! Get the sinitial concentration in tissue
	inline scalar_type C0t  () { return C0t_;  } const
	//! Get the sinitial concentration in network
	inline scalar_type C0v  () { return C0v_;  } const
	//! Get the radius at a given mesh_region
	//scalar_type R  (const getfem::mesh_im & mim, const size_type rg) { 
	//	return compute_radius(mim, mf_datav_, R_, rg);  
	//}
	//! Get the radius
	//vector_type & R (void) { return R_; }
	//! Get the vessel wall permeabilities
	vector_type & Q_pl (void) { return Q_pl_; }
	//! Get the Dahmkholer number 
	vector_type & Dalpha (void) { return Dalpha_; }
	//! Get the tissue diffusivity 
	vector_type & At (void) { return At_; }
	//! Get the vessel diffusivity 
	vector_type & Av (void) { return Av_; }
	//! Get the leakage of the capillary bed
	vector_type & Y  (void) { return Y_;  } const
	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const param3d1d_transp & param
		)
	{ 
		out << "--- PHYSICAL PARAMS ------" << endl;
		out << "  At     : "                << param.At_[0] << endl; 
		out << "  Av : "                << param.Av_[0] << endl; 
		out << "  Y      : "                << param.Y_[0] << endl; 
		out << "  Q_pl : "                << param.Q_pl_[0] << endl; 
		out << "  D_alpha : "                << param.Dalpha_[0] << endl; 
		out << "  T : "                << param.T_ << endl; 
		out << "  dt : "                << param.dt_ << endl; 
		out << "--------------------------" << endl;

		return out;            
	}

}; /* end of class */

} /* end of namespace */

#endif
