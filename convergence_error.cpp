/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   convergence_error.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016 - May 2018
  @brief  Methods of transport3d1d for computing convergence error
 */
 
 #include <transport3d1d.hpp>
 #include <AMG_Interface.hpp>
 #include <cmath>
 #include "gmm/gmm_inoutput.h"
 #include "getfem/getfem_import.h"


 namespace getfem {


   /* Reduced model: analysis of convergence error and model error
      Use an exact solution! */


// If you look for declaration of c_exact_function
// (That is the exact solution used for convergence)
// you'll find it in utilities_transp.hpp
///////////////////////////////////////////
// Function declaration for exact solution, source terms, etc
///////////////////////////////////////////
// declaration of some useful parameters for exact solution, source terms, etc
double C=1.0, k=1.0, R=1.0;

/////////////////////////////////////////////////////


// Used in convergence_error:
//! Exact concentration
double c_exact_func(const bgeot::base_node & x){
	if( sqrt((x[1]-0)*(x[1]-0)+(x[2]-0)*(x[2]-0)) <= R){
	return C*k/(1+k);}
	else{
	return C*k/(1+k)*(1-R*log(1/R*sqrt((x[1]-0)*(x[1]-0)+(x[2]-0)*(x[2]-0))));}
}



	//Assemble the reduced problem (with exact solution)
	void transport3d1d::assembly_reduced_transp (void){

	#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM_transp, UM_transp, FM_transp ..." << endl;
	#endif
	gmm::resize(AM_transp, dof_transp.Ct(), dof_transp.Ct());	gmm::clear(AM_transp);
	gmm::resize(UM_transp, dof_transp.Ct()); 			gmm::clear(UM_transp);
	gmm::resize(FM_transp, dof_transp.Ct()); 			gmm::clear(FM_transp);
	
	
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM_transp ..." << endl;
	#endif
	// Diffusion matrix for the interstitial problem
	sparse_matrix_type Dt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Dt);

	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Btt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Btt);
	// Vessel-to-tissue exchange matrix
	sparse_matrix_type Btv(dof_transp.Ct(), dof_transp.Cv());gmm::clear(Btv);
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Bvt);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Bvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Bvv);
	// Aux tissue-to-vessel averaging matrix
	sparse_matrix_type Mbar(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mbar);
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin(dof_transp.Cv(), dof_transp.Ct());gmm::clear(Mlin);
	


	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the stiffness matrix Dt ..." << endl;
	#endif	
	// Assemble Dt, stiffness matrix for laplacian
	getfem::asm_stiffness_matrix_for_homogeneous_laplacian(Dt,mimt,mf_Ct); 
	
	gmm::add(Dt,  
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 



	// Assemble Btt, Btv, Bvt and Bvv, coupling terms
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling aux exchange matrices Mbar and Mlin ..." << endl;
	#endif
	if(PARAM.int_value("couple", "flag for coupling function (notaro 0, brambilla 1)"))
	asm_exchange_aux_mat_transp(Mbar, Mlin, 
			mimv, mf_Ct, mf_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);

	if(!PARAM.int_value("couple", "flag for coupling function (notaro 0, brambilla 1)"))
	asm_exchange_aux_mat(Mbar, Mlin, 
			mimv, mf_Ct, mf_Cv, param.R(), descr.NInt);

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
	#endif
	k = PARAM.real_value("kk", "value of permeability for reduced problem");
	vector_type K(dof.coefv());
	gmm::copy(param.R(), K);
	gmm::scale(K, 2*pi*k);
	C = PARAM.real_value("Cv", "value of concentration in the vessel");
	R = param.R(0);
	bool NEWFORM = PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");
	
	asm_exchange_mat(Btt, Btv, Bvt, Bvv,
			mimv, mf_Cv, mf_coefv, Mbar, Mlin, K, NEWFORM);



	vector_type Cv(dof_transp.Cv(),C);
	gmm::mult(Btv, Cv, FM_transp);   // F = Btv*U

	// Copying Btt
	gmm::add(Btt,			 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); //A = Dt + Btt
	

	// De-allocate memory
	gmm::clear(Dt);    
	gmm::clear(Mbar);  gmm::clear(Mlin);
	gmm::clear(Btt);   gmm::clear(Btv);
	gmm::clear(Bvt);   gmm::clear(Bvv);

	// Assemble boundary conditions on tissue
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling boundary conditions ..." << endl;
	#endif
	scalar_type beta_t  = PARAM.real_value("BETAtissue_transp", "Coefficient for mixed BC for transport problem in tissue");
	vector_type c_ex(dof_transp.Ct());
	interpolation_function(mf_Ct, c_ex, c_exact_func ); // utilities_transp.hpp


		for (size_type bc=0; bc < BCt_transp.size(); ++bc) {
			if (BCt_transp[bc].label=="MIX") { // Robin BC
				vector_type BETA(mf_coeft.nb_dof(), beta_t);
				getfem::asm_mass_matrix_param(AM_transp, mimt, mf_Ct, mf_coeft, BETA,mf_Ct.linked_mesh().region(BCt_transp[bc].rg) );
			
			vector_type BETA_C0(mf_coeft.nb_dof(), beta_t*BCt_transp[bc].value);
			asm_source_term(FM_transp,mimt, mf_Ct, mf_coeft,BETA_C0);
			}
		}
		for (size_type bc=0; bc < BCt_transp.size(); ++bc) {
			if (BCt_transp[bc].label=="DIR") { // Dirichlet BC
				getfem::assembling_Dirichlet_condition(AM_transp, FM_transp, mf_Ct, BCt_transp[bc].rg, c_ex);				
			}
		} 
	}; //end of assembly_reduced_transp	


	 bool transport3d1d::solve_reduced_transp (void)
 	{

  	#ifdef M3D1D_VERBOSE_
	cout << "Solving the monolithic system ... " << endl;
	#endif
	double time = gmm::uclock_sec();


	gmm::csc_matrix<scalar_type> A_transp;
	gmm::clean(AM_transp, 1E-12);
	gmm::copy(AM_transp, A_transp);
	
	vector_type F_transp(gmm::vect_size(FM_transp));
	gmm::clean(FM_transp, 1E-12);
	gmm::copy(FM_transp, F_transp);
	
		
	if ( descr_transp.SOLVE_METHOD == "SuperLU" ) { // direct solver //
		#ifdef M3D1D_VERBOSE_
		cout << "  Applying the SuperLU method ... " << endl;
		#endif
		scalar_type cond;
		gmm::SuperLU_solve(A_transp, UM_transp, F_transp, cond);
		cout << "  Condition number (transport problem): " << cond << endl;
	}
else if (descr_transp.SOLVE_METHOD == "SAMG"){
	#ifdef WITH_SAMG	
	#ifdef M3D1D_VERBOSE_
		cout << "Solving the monolithic system ... " << endl;
	#endif



	//////////////////////////////////////AMG INTERFACE
	#ifdef M3D1D_VERBOSE_
	std::cout<<"converting A"<<std::endl;
	#endif
	gmm::csr_matrix<scalar_type> A_csr;
	gmm::clean(AM_transp, 1E-12);


	int dim_matrix=dof_transp.Ct();
	gmm::copy(gmm::sub_matrix(AM_transp,
			gmm::sub_interval(0 , dim_matrix),
			gmm::sub_interval(0 , dim_matrix)), A_csr);
	#ifdef M3D1D_VERBOSE_
	std::cout<<"converting X"<<std::endl;
	#endif
	std::vector<scalar_type> X,  B;

	gmm::resize(X,dim_matrix); gmm::clean(X, 1E-12);
	gmm::copy(gmm::sub_vector(UM_transp,gmm::sub_interval(0,dim_matrix)),X);

	#ifdef M3D1D_VERBOSE_
	std::cout<<"converting B"<<std::endl;
	#endif
	gmm::resize(B,dim_matrix);gmm::clean(B, 1E-12);
	gmm::copy(gmm::sub_vector(FM_transp,gmm::sub_interval(0,dim_matrix)),B);



	AMG amg("3d1d");
	amg.set_dof(dof_transp.Ct(),0,0,0);

	
	amg.convert_matrix(A_csr);
	amg.solve(A_csr, X , B , 1);
	gmm::copy(amg.getsol(),gmm::sub_vector(UM_transp,gmm::sub_interval(0,dim_matrix)));



	#ifdef SPARSE_INTERFACE
				for(int i = 0 ; i < nrows ; i++ ){U_1[i]=u[i];
					UM_transp[i]=u[i];	
				}
				gmm::copy(U_1, UM_transp);
	#endif
	#ifdef CSC_INTERFACE
				for(int i = 0 ; i < nnu ; i++ ){
					U_2[i]=u_samg[i];UM_transp[i]=u_samg[i];}
				gmm::copy(U_2,UM_transp);
	#endif
				
	#ifdef CSR_INTERFACE
				// for(int i = 0 ; i < nnu ; i++ ){
				//	U_2[i]=u_samg[i];UM[i]=u_samg[i];}
				 // gmm::copy(U_2,UM);
	#endif
	
	#else // with_samg=0
	std::cout<< "ERROR: you are trying to solve with samg, but WITH_SAMG=0"<<std::endl;
	std::cout<< "--> Call 'source configure.sh'and install the library again"<<std::endl;
	
	#endif
}
	else { // Iterative solver //

		// Iterations
		gmm::iteration iter(descr_transp.RES);  // iteration object with the max residu
		iter.set_noisy(1);               // output of iterations (2: sub-iteration)
		iter.set_maxiter(descr_transp.MAXITER); // maximum number of iterations

		// Preconditioners
		//! \todo Add preconditioner choice to param file
		// See \link http://download.gna.org/getfem/html/homepage/gmm/iter.html
		gmm::identity_matrix PM; // no precond
		//gmm::diagonal_precond<sparse_matrix_type> PM(AM); // diagonal preocond
		//gmm::ilu_precond<sparse_matrix_type> PM(AM);
		// ...
		//gmm::clear(AM);
		// See <http://download.gna.org/getfem/doc/gmmuser.pdf>, pag 15
	
		if ( descr_transp.SOLVE_METHOD == "CG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Conjugate Gradient method ... " << endl;
			#endif
			gmm::identity_matrix PS;  // optional scalar product
			gmm::cg(AM_transp, UM_transp, F_transp, PS, PM, iter);
		}
		else if ( descr_transp.SOLVE_METHOD == "BiCGstab" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the BiConjugate Gradient Stabilized method ... " << endl;
			#endif
			gmm::bicgstab(AM, UM, FM, PM, iter);
		}
		else if ( descr_transp.SOLVE_METHOD == "GMRES" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Generalized Minimum Residual method ... " << endl;
			#endif
			size_type restart = 50;
			gmm::gmres(A_transp, UM, FM, PM, restart, iter);
		}
		else if ( descr_transp.SOLVE_METHOD == "QMR" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the Quasi-Minimal Residual method ... " << endl;
			#endif
			gmm::qmr(AM, UM, FM, PM, iter);
		}
		else if ( descr_transp.SOLVE_METHOD == "LSCG" ) {
			#ifdef M3D1D_VERBOSE_
			cout << "  Applying the unpreconditionned Least Square CG method ... " << endl;
			#endif
			gmm::least_squares_cg(AM, UM, FM, iter);
		}
		// Check
		if (iter.converged())
			cout << "  ... converged in " << iter.get_iteration() << " iterations." << endl;
		else if (iter.get_iteration() == descr_transp.MAXITER)
			cerr << "  ... reached the maximum number of iterations!" << endl;

	}
	
	//export solution
	#ifdef M3D1D_VERBOSE_
	std::cout<<"solved! going to export..."<<std::endl;
	#endif	
	
	cout << endl<<"... time to solve " << gmm::uclock_sec() - time << " seconds\n";				
	
	return true;
 	}; // end of solve_transp









 void transport3d1d::export_vtk_reduced_transp (const string & suff)
 {
  if (PARAM.int_value("VTK_EXPORT"))
  {
	#ifdef M3D1D_VERBOSE_
	cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Saving the results from the monolithic unknown vector ... " << endl;
	#endif
	
	// Array of unknown dof of the interstitial velocity
	vector_type Ct(dof_transp.Ct()); 

	// Array of unknown dof of the network velocity
	scalar_type C = PARAM.real_value("Cv", "value of concentration in the vessel");
	vector_type Cv(dof_transp.Cv(),C);

	//Copy solution
	gmm::copy(gmm::sub_vector(UM_transp, 
		gmm::sub_interval(0, dof_transp.Ct())), Ct);


	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Ct ..." << endl;
	#endif
	vtk_export exp_Ct(descr_transp.OUTPUT+"Ct"+suff+"_t"+".vtk");
	exp_Ct.exporting(mf_Ct);
	exp_Ct.write_mesh();
	exp_Ct.write_point_data(mf_Ct, Ct, "Ct");



	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Cv ..." << endl;
	#endif
	vtk_export exp_Cv(descr_transp.OUTPUT+"Cv"+suff+"_t"+".vtk");
	exp_Cv.exporting(mf_Cv);
	exp_Cv.write_mesh();
	exp_Cv.write_point_data(mf_Cv, Cv, "Cv");

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
	#endif

	
  }
 }; // end of export_transp

 	// Compute the norm of error (with exact solution)
	void transport3d1d::compute_error_reduced_transp (void){

	double time = gmm::uclock_sec();
	
	std::string FEM_TYPET_CEX = PARAM.string_value("FEM_TYPET_CEX","FEM 3D tissue - exact concentration");
	pfem pf_Ct_ex = fem_descriptor(FEM_TYPET_CEX);
	mesh_fem mf_Ct_ex(mesht);
	mf_Ct_ex.set_finite_element(mesht.convex_index(), pf_Ct_ex);

	vector_type c_ex(mf_Ct_ex.nb_dof());
	interpolation_function(mf_Ct_ex, c_ex, c_exact_func ); // utilities_transp.hpp

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Cv ..." << endl;
	#endif
	vtk_export exp_Ct_ex(descr_transp.OUTPUT+"Ct_ex"+".vtk");
	exp_Ct_ex.exporting(mf_Ct_ex);
	exp_Ct_ex.write_mesh();
	exp_Ct_ex.write_point_data(mf_Ct_ex, c_ex, "Ct_exact");

	mesh_im mimt_ex(mesht);
	pintegration_method pim_t_ex = int_method_descriptor(PARAM.string_value("IM_TYPET_CEX","Integration method for L2/H1 norm of exact solution"));
	mimt_ex.set_integration_method(mesht.convex_index(), pim_t_ex);


	cout << endl<<"... time to interpolate the exact solution:  " << gmm::uclock_sec() - time << " second\n";

	time = gmm::uclock_sec();

if (PARAM.int_value("NORM","flag if you want to compute the norm only on a face or in the volume")==0)
{ // compute the norm on a single face, with ch in P1, cex in P3
std::vector<scalar_type> L2v(1);
std::vector<scalar_type> H1v(1);

	getfem::generic_assembly
	assemL2("Ch=data$1(#1);" "Cex=data$2(#2);"
		"t1=comp(Base(#1).Base(#1))(i,j).Ch(i).Ch(j);"
		"t2=comp(Base(#2).Base(#2))(i,j).Cex(i).Cex(j);"
		"t3=comp(Base(#1).Base(#2))(i,j).Ch(i).Cex(j);"
	"V$1() += t1 + t2 - t3-t3;");
	assemL2.push_mi(mimt_ex);
	assemL2.push_mf(mf_Ct);
	assemL2.push_mf(mf_Ct_ex);
	assemL2.push_data(UM_transp);
	assemL2.push_data(c_ex);
	assemL2.push_vec(L2v);
	assemL2.assembly(1);   //region 0 is the face x=0			
	

cout << endl<<"... time to compute the L2 norm:  " << gmm::uclock_sec() - time << " second\n";	
time = gmm::uclock_sec();

	getfem::generic_assembly
	assemH1("Ch=data$1(#1);" "Cex=data$2(#2);"
		"t1=comp(Grad(#1).Grad(#1))(i,p,j,p).Ch(i).Ch(j);"
		"t2=comp(Grad(#2).Grad(#2))(i,p,j,p).Cex(i).Cex(j);"
		"t3=comp(Grad(#1).Grad(#2))(i,p,j,p).Ch(i).Cex(j);"
	"V$1() += t1 + t2 - t3-t3;");
	assemH1.push_mi(mimt_ex);
	assemH1.push_mf(mf_Ct);
	assemH1.push_mf(mf_Ct_ex);
	assemH1.push_data(UM_transp);
	assemH1.push_data(c_ex);
	assemH1.push_vec(H1v);
	assemH1.assembly(1);   //region 0 is the face x=0	

cout << endl<<"... time to compute the H1 norm:  " << gmm::uclock_sec() - time << " second\n\n";	
	scalar_type L2_norm=L2v[0];
	scalar_type H1_norm=H1v[0];
	H1_norm+=L2_norm;

	L2_norm=sqrt(L2_norm);
	H1_norm=sqrt(H1_norm);
	
	cout <<"norm computed on a face: "<<endl;
	cout << "L2 norm error = " << L2_norm << endl;
	cout << "H1 norm error = " << H1_norm << endl; 

}
else if (PARAM.int_value("NORM","flag if you want to compute the norm only on a face or in the volume")==1)
{
// compute the norm on all volume, with ch in P1, cex in P3
std::vector<scalar_type> L2v(1);
std::vector<scalar_type> H1v(1);

	getfem::generic_assembly
	assemL2("Ch=data$1(#1);" "Cex=data$2(#2);"
		"t1=comp(Base(#1).Base(#1))(i,j).Ch(i).Ch(j);"
		"t2=comp(Base(#2).Base(#2))(i,j).Cex(i).Cex(j);"
		"t3=comp(Base(#1).Base(#2))(i,j).Ch(i).Cex(j);"
	"V$1() += t1 + t2 - t3-t3;");
	assemL2.push_mi(mimt_ex);
	assemL2.push_mf(mf_Ct);
	assemL2.push_mf(mf_Ct_ex);
	assemL2.push_data(UM_transp);
	assemL2.push_data(c_ex);
	assemL2.push_vec(L2v);
	assemL2.assembly();   //region 0 is the face x=0			
	

cout << endl<<"... time to compute the L2 norm:  " << gmm::uclock_sec() - time << " second\n";	
time = gmm::uclock_sec();

	getfem::generic_assembly
	assemH1("Ch=data$1(#1);" "Cex=data$2(#2);"
		"t1=comp(Grad(#1).Grad(#1))(i,p,j,p).Ch(i).Ch(j);"
		"t2=comp(Grad(#2).Grad(#2))(i,p,j,p).Cex(i).Cex(j);"
		"t3=comp(Grad(#1).Grad(#2))(i,p,j,p).Ch(i).Cex(j);"
	"V$1() += t1 + t2 - t3-t3;");
	assemH1.push_mi(mimt_ex);
	assemH1.push_mf(mf_Ct);
	assemH1.push_mf(mf_Ct_ex);
	assemH1.push_data(UM_transp);
	assemH1.push_data(c_ex);
	assemH1.push_vec(H1v);
	assemH1.assembly();   //region 0 is the face x=0	

cout << endl<<"... time to compute the H1 norm:  " << gmm::uclock_sec() - time << " second\n\n";	
	scalar_type L2_norm=L2v[0];
	scalar_type H1_norm=H1v[0];
	H1_norm+=L2_norm;

	L2_norm=sqrt(L2_norm);
	H1_norm=sqrt(H1_norm);
	
	cout <<"norm computed on the volume: "<<endl;
	cout << "L2 norm error = " << L2_norm << endl;
	cout << "H1 norm error = " << H1_norm << endl; 
}
else if (PARAM.int_value("NORM","flag if you want to compute the norm only on a face or in the volume")==2)
{
// compute the norm on a single face, with ch in P3, cex in P3
	vector_type c_error(mf_Ct_ex.nb_dof());
	getfem::interpolation(mf_Ct, mf_Ct_ex, UM_transp,c_error);
	gmm::add(c_error, gmm::scaled(c_ex, -1.0), c_error); // V1 - 1.0 * V2 --> V2	
	
std::vector<scalar_type> L2v(1);
std::vector<scalar_type> H1v(1);

	getfem::generic_assembly
	assemL2("C_error=data$1(#1);"
	"V$1() += comp(Base(#1).Base(#1))(i,j).C_error(i).C_error(j);");
	assemL2.push_mi(mimt_ex);
	assemL2.push_mf(mf_Ct_ex);
	assemL2.push_data(c_error);
	assemL2.push_vec(L2v);
	assemL2.assembly(1);   //region 0 is the face x=0			
	

cout << endl<<"... time to compute the L2 norm:  " << gmm::uclock_sec() - time << " second\n";	
time = gmm::uclock_sec();

	getfem::generic_assembly
	assemH1("C_error=data$1(#1);"
	"V$1() += comp(Grad(#1).Grad(#1))(i,p,j,p).C_error(i).C_error(j);");
	assemH1.push_mi(mimt_ex);
	assemH1.push_mf(mf_Ct_ex);
	assemH1.push_data(c_error);
	assemH1.push_vec(H1v);
	assemH1.assembly(1);   //region 0 is the face x=0	

cout << endl<<"... time to compute the H1 norm:  " << gmm::uclock_sec() - time << " second\n\n";	
	scalar_type L2_norm=L2v[0];
	scalar_type H1_norm=H1v[0];
	H1_norm+=L2_norm;

	L2_norm=sqrt(L2_norm);
	H1_norm=sqrt(H1_norm);
	
	cout <<"norm computed on a face: "<<endl;
	cout << "L2 norm error = " << L2_norm << endl;
	cout << "H1 norm error = " << H1_norm << endl; 
}
else if (PARAM.int_value("NORM","flag if you want to compute the norm only on a face or in the volume")==3)
{
// compute the norm on all the volume, with ch in P3, cex in P3
	vector_type c_error(mf_Ct_ex.nb_dof());
	getfem::interpolation(mf_Ct, mf_Ct_ex, UM_transp,c_error);
	gmm::add(c_error, gmm::scaled(c_ex, -1.0), c_error); // V1 - 1.0 * V2 --> V2	
	
std::vector<scalar_type> L2v(1);
std::vector<scalar_type> H1v(1);

	getfem::generic_assembly
	assemL2("C_error=data$1(#1);"
	"V$1() += comp(Base(#1).Base(#1))(i,j).C_error(i).C_error(j);");
	assemL2.push_mi(mimt_ex);
	assemL2.push_mf(mf_Ct_ex);
	assemL2.push_data(c_error);
	assemL2.push_vec(L2v);
	assemL2.assembly();   //region 0 is the face x=0			
	

cout << endl<<"... time to compute the L2 norm:  " << gmm::uclock_sec() - time << " second\n";	
time = gmm::uclock_sec();

	getfem::generic_assembly
	assemH1("C_error=data$1(#1);"
	"V$1() += comp(Grad(#1).Grad(#1))(i,p,j,p).C_error(i).C_error(j);");
	assemH1.push_mi(mimt_ex);
	assemH1.push_mf(mf_Ct_ex);
	assemH1.push_data(c_error);
	assemH1.push_vec(H1v);
	assemH1.assembly();   //region 0 is the face x=0	

cout << endl<<"... time to compute the H1 norm:  " << gmm::uclock_sec() - time << " second\n\n";	
	scalar_type L2_norm=L2v[0];
	scalar_type H1_norm=H1v[0];
	H1_norm+=L2_norm;

	L2_norm=sqrt(L2_norm);
	H1_norm=sqrt(H1_norm);
	
	cout <<"norm computed on the volume: "<<endl;
	cout << "L2 norm error = " << L2_norm << endl;
	cout << "H1 norm error = " << H1_norm << endl; 
}

//interpolate exact function
	vector_type c_ex_P1(mf_Ct.nb_dof());
	interpolation_function(mf_Ct, c_ex_P1, c_exact_func );

	gmm::add(UM_transp, gmm::scaled(c_ex_P1, -1.0), c_ex_P1); // V1 - 1.0 * V2 --> V2
//export ch-cex
	vtk_export exp_err(descr_transp.OUTPUT+"Ct_error"+".vtk");
	exp_err.exporting(mf_Ct);
	exp_err.write_mesh();
	exp_err.write_point_data(mf_Ct, c_ex_P1, "Ct_error");


}; 





 
 } // end of namespace
