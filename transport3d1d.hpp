/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   transport3d1d.hpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016 - September 2018.
  @brief  Declaration of the main class for the 3D/1D coupled transport problem.
 */
 
#ifndef M3D1D_TRANSPORT3D1D_HPP_
#define M3D1D_TRANSPORT3D1D_HPP_


// GetFem++ libraries

#include <gmm/gmm_matrix.h>
#include <getfem/dal_bit_vector.h>

// Project headers
#include <problemHT.hpp>

#include <assembling1d_transp.hpp>          
#include <assembling3d_transp.hpp> 
#include <assembling3d1d_transp.hpp>        
#include <dof3d1d_transp.hpp>
#include <descr3d1d_transp.hpp>
#include <param3d1d_transp.hpp>
#include <utilities_transp.hpp>
#include <node_transp.hpp>
#include "../utilities/muparser/include/muParser.h"


 namespace getfem {

//!	Main class defining the coupled 3D/1D transport problem.
class transport3d1d: public problemHT { 

public:
	transport3d1d(void) : 
		mf_Ct(mesht), mf_Cv(meshv),mf_Ct_Omega(mesht),mf_Ct_Sigma(mesht){} 
	
	// Main methods of class: implement standard and complete transport problem
	//! Initialize the transport problem
	void init_transp (int argc, char *argv[]);	
	//! Assemble the transport problem
	void assembly_transp (void);
	//! Solve the transport problem
	bool solve_transp (void);
	//! Export the transport solution
	const void export_vtk_transp (const string & time_suff = "",const string & suff = "");
	//! Compute residuals for mass balance at each junction
	void mass_balance (void);	

	//! Getter for solution
	inline vector_type get_UM(void) {return UM_transp;};

	//Aux methods for interface with problem3d1d class
	//! Initialize the fluid problem
	void init_fluid (int argc, char *argv[]);
	//! Assemble the fluid problem
	void assembly_fluid (void);
	//! Solve the fluid problem
	bool solve_fluid (void);
	//! Export the fluid solution
	const void export_vtk_fluid (const string & suff = "");


	// Methods for reduced problem with exact solution (compute convergence error)
	//! Assemble the reduced problem (with exact solution)
	void assembly_reduced_transp (void); 
	//! Solve the reduced problem (with exact solution)
	bool solve_reduced_transp (void);
 	//! Export the reduced problem (with exact solution)
	void export_vtk_reduced_transp (const string & suff = ""); 
 	//! Compute the norm of error (with exact solution)
	void compute_error_reduced_transp (void); 


	// Methods for reduced problem  (compute model error with dual solution)
	//Method for reduced problem
	void model_error(int argc, char *argv[]); 
	//Method for reference problem, assumption A0 (That is the initial 3D-3D reference problem)
	void model_error_A0(int argc, char *argv[]); 
	//Method for reduced problem, assumption A1
	void model_error_A1(int argc, char *argv[]); 
	//Method for reduced problem, assumption A2
	void model_error_A2(int argc, char *argv[]); 
	//Method for reduced problem, assumption A3 (That is the final 3D-1D reduced problem)
	void model_error_A3(int argc, char *argv[]); 


	//! Initialize reduced problem
	void init_model(int argc, char *argv[]);
	//! Assemble reduced primal problem
	void assembly_model(const size_type ASSUMPTION, const size_type VERSION);
	//! Solve reduced problem 
	bool solve_model(const size_type ASSUMPTION, const size_type VERSION);
	//! Export reduced problem
	void export_model(const size_type ASSUMPTION, const size_type VERSION, const string & suff = "");

	//! Compute model error of A1	
	void compute_model_error_1(const size_type DUALMODEL);
	//! Compute model error of A2	
	void compute_model_error_2(const size_type DUALMODEL);
	//! Compute model error of A3	
	void compute_model_error_3(const size_type DUALMODEL);


	
protected:
	 
	//! Finite Element Method for the tissue concentration @f$c_t@f$
	mesh_fem mf_Ct; 
	//! Finite Element Method for the vessel concentration @f$c_v@f$
	mesh_fem mf_Cv; 

	//! Finite Element Method for the tissue concentration @f$c_t@f$
	mesh_fem mf_Ct_Omega; 
	//! Finite Element Method for the tissue concentration @f$c_t@f$
	mesh_fem mf_Ct_Sigma;

	//! Algorithm description strings (mesh files, FEM types, solver info, ...) 
	descr3d1d_transp descr_transp;
	//! Physical parameters
	param3d1d_transp param_transp;
	//! Number of degrees of freedom
	dof3d1d_transp dof_transp;
		
	
	//! List of BC nodes of the network
	vector< node > BCv_transp;	
	//! List of BC nodes of the tissue
	vector< node > BCt_transp;
	//! List of junction nodes of the network
	vector< node_transp > Jv_transp;

		
	//! Monolithic matrix for the coupled problem
	sparse_matrix_type AM_transp;
	//! Monolithic array of unknowns for the coupled problem
	vector_type        UM_transp;
	//! Monolithic right hand side for the coupled problem
	vector_type        FM_transp;
	
	
	//! Monolithic temporary matrix for update
	sparse_matrix_type AM_temp;
	//! Monolithic temporary right hand side for update
	vector_type        FM_temp; 
	
	// Primal Solution, reference problem
	vector_type U0;	
	// Primal Solution, assumption A1
	vector_type U1;
	// Primal Solution, assumption A2
	vector_type U2;
	// Primal Solution, assumption A3
	vector_type U3;

	
	// Dual Solution, reference problem
	vector_type Z0;
	// Dual Solution, assumption A1
	vector_type Z1;
	// Dual Solution, assumption A2
	vector_type Z2;
	// Dual Solution, assumption A3
	vector_type Z3;

	// Global residuals
	vector_type ETA_RESIDUAL;
	// Global estimator
	vector_type ETA;

	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type MLIN;
	// Aux tissue-to-vessel average matrix (circumsference)
	sparse_matrix_type MBAR;
	// Aux tissue-to-vessel average matrix (section) 
	sparse_matrix_type MBARBAR;


	// Aux methods for init
	//! Import algorithm specifications
	void import_data_transp(void);
	//! Import mesh for tissue (3D) and vessel (1D)  
	void build_mesh_transp(void); 
	//! Set finite elements methods and integration methods 
	void set_im_and_fem_transp(void);
	//! Build problem parameters
	void build_param_transp(void);
	//! Build the list of tissue boundary data 
	/*!	Face numbering:
		  0 : {x = 0 }  "back"
		  1 : {x = Lx}  "front"
		  2 : {y = 0 }  "left"
		  3 : {y = Ly}  "right"
		  4 : {z = 0 }  "bottom"
		  5 : {z = Lz}  "top"
	 */
	void build_tissue_boundary_transp(void);
	//! Build the list of vessel boundary (and junctions) data 
	void build_vessel_boundary_transp(void);
	
	//Aux method for assembly
	//! Build the monolithic matrix AM_transp by blocks
	void assembly_mat_transp(void);
	//! Build the monolithic rhs FM_transp by blocks
	void assembly_rhs_transp(void);
	
	//Aux method for solve 
	//! Aux function for update of rhs at each time step
	void update_transp(void);
	//! Aux function for solver: contains the different solving methods and actually solve the system
	bool solver (	const size_type dof1=0, 
			const size_type dof2=0, 
			const size_type dof3=0,
			const size_type dof4=0);	

	//! Set finite elements and integration methods: mf_t should be defined only on Omega_plus	
	void set_im_and_fem_model(void);



	
}; //end of class trasport3d1d


}  //end of namespace

#endif
