/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2015-2016
                  
                Copyright (C) 2018 Stefano Brambilla
======================================================================*/
/*! 
  @file   model_error.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   July 2018
  @brief  Methods of the main class for computing model error.
 */
 
 #include <transport3d1d.hpp>
 #include <AMG_Interface.hpp>
 #include <cmath>
 #include "gmm/gmm_inoutput.h"
 #include "getfem/getfem_import.h"

 namespace getfem {

 
//////////////////////////////////////////////////
// Some global variables
/* PRIMAL A0 = reference solution (after A0 = small radius)
 * PRIMAL A1 = reduced solution of A1 (vessel 1D)
 * PRIMAL A2 = reduced solution of A2 (3D domain = Omega, not Omega_plus)
 * PRIMAL A3 = reduced solution of A3 (transmission condition. U and u solutions)

 * DUAL Z0 = reference dual solution z1ref (vessel 3D)
 * DUAL Z1 = reference dual solution z2ref (3d domain = Omega_plus)
 * DUAL Z2 = reference dual solution z3ref (3d domain = Omega; transmission not neglected)
 * DUAL Z3 = reduced dual suolutions z and Z (after all assumptions)

 * NB: Ai-th solution means that we have already considered all the reductions from 0 to i. 
 */
const int A0=0;
const int A1=1;
const int A2=2;
const int A3=3;

const int PRIMAL=1;
const int DUAL=2;

const int REFERENCE =1;
const int REDUCED = 2;


///////////////////////////////////////////
// declaration of some useful parameters for exact solution, source terms, etc
double CC=1.0, kk=1.0, RR=1.0, constant_value=1;
std::string expr="", constant_expr="CONSTANT_VALUE";

// Function declaration for source terms, etc

/* In this functions, we use the useful muParser library.
 * For further informations, check the following link: 
 * http://beltoforion.de/article.php?a=muparser
 * Anyway, you can pass whatever expression do you need, 
 * with a wide range of possiblities, such as:
 * sin,cos,tan,log,ln,log2,log10,exp,sqrt,sign,abs,
 * rint (round to integer),min,max,sum,avg,
 * =,==,&&,||,<=,>=,+,-,*,/,^.
 */

//! Function for source term in vessels
double g_function(const bgeot::base_node & x){

	//Define p
	mu::Parser p;
	//Define variables
		/* We define both coordinates x,y,z,
		 * which are needed for the function,
		 * and other useful parameters, namely C,R,k,
		 * that you can pass in the .param file.
		 * Note that you can add in this file any other 
		 * variable you want. */ 
	double xx=x[0], yy=x[1], zz=x[2];
	p.DefineVar("x", &xx);
	p.DefineVar("y", &yy);
	p.DefineVar("z", &zz);
	p.DefineVar("C", &CC);	
	p.DefineVar("k", &kk);
	p.DefineVar("R", &RR);
	//p.DefineVar("new variable", double variable);	
	
	//Define the expression
		/* The string should be already modified
		 * in the code */
	p.SetExpr(expr);
	return p.Eval();
}

// Used in model_error:
//! Function for source term in tissue 
double f_function(const bgeot::base_node & x){

	//Define p
	mu::Parser p;
	//Define variables
		/* We define both coordinates x,y,z,
		 * which are needed for the function,
		 * and other useful parameters, namely C,R,k,
		 * that you can pass in the .param file.
		 * Note that you can add in this file any other 
		 * variable you want. */ 
	double xx=x[0], yy=x[1], zz=x[2];
	p.DefineVar("x", &xx);
	p.DefineVar("y", &yy);
	p.DefineVar("z", &zz);
	p.DefineVar("C", &CC);	
	p.DefineVar("k", &kk);
	p.DefineVar("R", &RR);
	//p.DefineVar("new variable", double variable);	
	
	//Define the expression
		/* The string should be already modified
		 * in the code */
	p.SetExpr(expr);
	return p.Eval();
}

//! Function for source term in tissue 
double constant_function(const bgeot::base_node & x){

	//Define p
	mu::Parser p;
	//Define variables
		/* We define both coordinates x,y,z,
		 * which are needed for the function,
		 * and other useful parameters, namely C,R,k,
		 * that you can pass in the .param file.
		 * Note that you can add in this file any other 
		 * variable you want. */ 
	double xx=x[0], yy=x[1], zz=x[2];
	p.DefineVar("x", &xx);
	p.DefineVar("y", &yy);
	p.DefineVar("z", &zz);
	p.DefineVar("CONSTANT_VALUE", &constant_value);
	//p.DefineVar("new variable", double variable);	
	
	//Define the expression
		/* The string should be already modified
		 * in the code */
	p.SetExpr(constant_expr);
	return p.Eval();
}

	/////////////////////////////////////////////////////
	// Methods for computing reduced (primal and dual) solution

	//! Initialize problem U1
	void transport3d1d::init_model(int argc, char *argv[]){

	PARAM.read_command_line(argc, argv);
	//1. Import data (algorithm specifications, boundary conditions, ...)	
	import_data_transp();
	//2. Import mesh for tissue (3D) and vessel network (1D)
	build_mesh_transp();
	//3. Set finite elements and integration methods
	set_im_and_fem_model();
 	//4. Build problem parameters
	build_param_transp();
	//5. Build the list of tissue boundary data
	build_tissue_boundary_transp();
	//6. Build the list of tissue boundary (and junction) data
 	build_vessel_boundary_transp();

	gmm::resize(ETA_RESIDUAL, dof_transp.Ct());
	gmm::resize(ETA		, dof_transp.Ct());
};
	//3. Set finite elements and integration methods
	void transport3d1d::set_im_and_fem_model(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEMs for tissue and vessel problems ..." << endl;
	#endif
	
	
	pfem pf_Ct = fem_descriptor(descr_transp.FEM_TYPET_C);
	pfem pf_Cv = fem_descriptor(descr_transp.FEM_TYPEV_C); 

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for tissue ..." << endl;
	#endif
		
/*! \todo Check the elements of Gamma.
Using gmesh, one can easily build the parallelepiped Omega, on which are defined two regions,
Sigma and Omega_plus, that are the physical vessel and the tissue surrounding it. (use Physical_volume(rg))
Nevertheless, GetFEM goes banana if you use a mesh composed of both thetrahedra AND triangles: 
so you cannot use Physical_surface(rg) to build the surface Gamma, that is the vessel wall, made of the triangles
(faces) in common between Sigma and Omega_plus.
The problem is that in GetFEM the faces have not a general index, as a matter of fact only the elements (tetrahedra)
has an index for storage. For example, a mesh of 1000 tetrahedra has an index going from 0 to 999. Every tetrahedron
has his 4 faces numbered from 0 to three.
If tetrahedron 15 has a face in common with tetrahedron 347, those are, maybe, faces 2 of element 15 and
face 1 of element 347. GetFEM doesn't know they are the same face; furthermore, they do not share the same basises!!
Integral on those two faces (for GetFEM are different faces!) they result in different values.
For out interest, the "rightest" thing to do seemed to create a region Gamma disposing of both the version of each face (coming from the inner and the outer tetrahedron). This might be incorrect, though. The following code builds different integrals using different definitions of the surface Gamma; those results convinced us that it was better considering both inner and outer faces. (@s.brambilla93@gmail.com)
EDIT: the phenomenon described before arises when integrating only on Sigma, only on Omega_plus, or on the whole Omega.
If we define Gamma to be the faces of the elements of Sigma, the integrals on Sigma and on Omega coincides, but 
all the integrals on Omega_plus are nulls.
Viceversa, if we define Gamma to be the faces of the elements of Omega_plus, the integrals on Omega_plus and on Omega coincides, but all the integrals on Sigma are nulls.
Finally, defining Gamma to be both the faces of elements in Sigma and Omega_t, give correct results on integrals on Sigma and Omega_t, but the integrals on the whole Omega are exactly* doubled. Therefore, as the following test validate, it is correct the definition of Gamma with both the surfaces, multiplying for 0.5 every matrix coming from
the integral on Gamma from the whole Omega domain. 

*Altough "exactly" seems a strong word, looking at F2b_ and F2c_ in check=2, we can be sure enough of this assumption.

*/

	//Define mf_Ct_Omega and mf_Ct_Sigma
	mf_Ct.set_finite_element(mesht.convex_index(), pf_Ct);
	mf_Ct_Omega.set_finite_element(mesht.region(descr_transp.OMEGA).index(), pf_Ct);
	mf_Ct_Sigma.set_finite_element(mesht.region(descr_transp.SIGMA).index(), pf_Ct); 

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for vessel branches ..." << endl;
	#endif

	mf_Cv.set_finite_element(meshv.convex_index(), pf_Cv);

	
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEM dimensions for tissue and vessel problems ..." << endl;
	#endif

	dof_transp.set(mf_Ct, mf_Cv, mf_Ct_Omega, mf_Ct_Sigma);
	#ifdef M3D1D_VERBOSE_
	cout << std::scientific << dof_transp;
	#endif

	mimv.clear();
	mf_Uvi.clear();
	mf_Pv.clear();
	mf_coefv.clear();
	mf_coefvi.clear();

	mimt.clear();
	mf_Ut.clear();
	mf_Pt.clear();
	mf_coeft.clear();
	problem3d1d::set_im_and_fem();


};

template<typename VEC>
void 
asm_flux
	(VEC & V,
	 const mesh_im & mim,
	 const mesh_fem & mf_c,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 	
	 
	 {
	generic_assembly 
	assem1("V$1(#1)+=comp(Grad(#1).Normal())(:,k,k);");
	assem1.push_mi(mim);
	assem1.push_mf(mf_c);
	assem1.push_vec(V);
	assem1.assembly(rg);
};



//! Assemble elliptic problem for compute modeling error, for each assumption and version
	void transport3d1d::assembly_model(const size_type ASSUMPTION,const size_type VERSION=PRIMAL){

	// define dofs:
	/* A0: 3d defined on mf_Ct_Omega;	1d defined on mf_Ct_Sigma
	 * A1: 3d defined on mf_Ct_Omega;	1d defined on mf_Cv
	 * A2: 3d defined on mf_Ct;		1d defined on mf_Cv 
	 * A3: 3d defined on mf_Ct;		1d defined on mf_Cv
	 */

	size_type dof_t;
	size_type dof_v;

	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		dof_t  = dof_transp.Ct_Omega();	}
	else{
		dof_t  = dof_transp.Ct();	}

 	if(ASSUMPTION==A0){
		dof_v  = dof_transp.Ct_Sigma();	}
	else{
		dof_v  = dof_transp.Cv();	}

	size_type dof_tot= dof_t+dof_v;

	#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM_transp, UM_transp, FM_transp ..." << endl;
	#endif
	gmm::resize(AM_transp, dof_tot, dof_tot);	gmm::clear(AM_transp);
	gmm::resize(UM_transp, dof_tot);		gmm::clear(UM_transp);
	gmm::resize(FM_transp, dof_tot);	 	gmm::clear(FM_transp);
	
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM_transp ..." << endl;
	#endif
	// Mass(time derivative)  matrix for the interstitial problem
	sparse_matrix_type At(dof_t, dof_t);gmm::clear(At);
	// Mass (time derivative)  matrix for the network problem
	sparse_matrix_type Av(dof_v, dof_v);gmm::clear(Av);	

	// Tissue-to-tissue exchange matrix on Gamma
	sparse_matrix_type Mtt(dof_t, dof_t);gmm::clear(Mtt);
	// Tissue-to-vessel exchange matrix on Gamma
	sparse_matrix_type Mtv(dof_t, dof_v);gmm::clear(Mtv);
	// Vessel-to-tissue exchange matrix on Gamma
	sparse_matrix_type Mvt(dof_v, dof_t);gmm::clear(Mtv);
	// Vessel-to-vessel exchange matrix on Gamma
	sparse_matrix_type Mvv(dof_v, dof_v);gmm::clear(Mvv);

	// Tissue-to-tissue exchange matrix
	sparse_matrix_type Btt(dof_t, dof_t);gmm::clear(Btt);
	// Vessel-to-tissue exchange matrix
	sparse_matrix_type Btv(dof_t, dof_v);gmm::clear(Btv);
	// Tissue-to-vessel exchange matrix
	sparse_matrix_type Bvt(dof_v, dof_t);gmm::clear(Bvt);
	// Vessel-to-vessel exchange matrix
	sparse_matrix_type Bvv(dof_v, dof_v);gmm::clear(Bvv);

	// Aux tissue-to-vessel averaging matrix (Circumference)
	sparse_matrix_type Mbar(dof_v, dof_t);gmm::clear(Mbar);
	// Aux tissue-to-vessel interpolation matrix
	sparse_matrix_type Mlin(dof_v, dof_t);gmm::clear(Mlin);
	// Aux tissue-to-vessel averaging matrix (Section)
	sparse_matrix_type Mbarbar(dof_v, dof_transp.Ct());gmm::clear(Mbarbar); //Whatever the case, this matrix should be defined on the whole 3d mesh
	// Aux tissue source vector
	vector_type Ft(dof_t); gmm::clear(Ft);
	// Aux vessels source vector
	vector_type Fv(dof_v); gmm::clear(Fv);
	// Aux tissue source vector
	vector_type Jt(dof_t); gmm::clear(Jt);
	// Aux vessels source vector
	vector_type Jv(dof_v); gmm::clear(Jv);

	//Assemble Perimeter, Area, and build Kt and Kv
	#ifdef M3D1D_VERBOSE_
	cout << "Assemble Perimeter, Area, and build Kt and Kv ..." << endl;
	#endif
	//Permeability
	scalar_type k = PARAM.real_value("kk", "value of permeability for reduced problem");
	vector_type Kv(dof.coefv());
	vector_type Kt(dof.coeft());
	constant_value=k;
	interpolation_function(mf_coefv, Kv, constant_function ); 
	interpolation_function(mf_coeft, Kt, constant_function ); 

	//Perimeter = 2*pi*R*k 
	vector_type Perimeter(dof.coefv());
	gmm::copy(param.R(), Perimeter);
	gmm::scale(Perimeter, 2*pi*k);

	//Area = pi*R^2
	vector_type Area(dof.coefv());
	gmm::copy(param.R(), Area);
	gmm::vscale(param.R(), Area);
	gmm::scale(Area, pi); 


	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the stiffness matrix At ..." << endl;
	#endif	
	// Assemble At, stiffness matrix for laplacian in tissue
if(VERSION==PRIMAL){
	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At,mimt,mf_Ct_Omega); 	}
	else{
		getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At,mimt,mf_Ct      );    }
}else if(VERSION==DUAL)	{
	if(ASSUMPTION==A0){
		getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At,mimt,mf_Ct_Omega);
		/*do nothing! A0 only vessel */  /* add something non-zero*/;			}
	else if(ASSUMPTION==A1){
		getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At,mimt,mf_Ct_Omega); 	}
	else{
		getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At,mimt,mf_Ct      );    }
 

}
	gmm::add(At,  
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t), 
					gmm::sub_interval(0, dof_t))); 


	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the stiffness matrix Av ..." << endl;
	#endif	
	// Assemble Av, stiffness matrix for laplacian in vessel

if(VERSION==PRIMAL){
	if(ASSUMPTION==A0){
		getfem::asm_stiffness_matrix_for_homogeneous_laplacian(Av,mimt,mf_Ct_Sigma); 	}
	else{
		getfem::asm_stiffness_matrix_for_laplacian(Av,mimv,mf_Cv,mf_coefv,Area);  	}

}else if(VERSION==DUAL)	{
	if(ASSUMPTION==A0){
		getfem::asm_stiffness_matrix_for_homogeneous_laplacian(Av,mimt,mf_Ct_Sigma);	}
	else if(ASSUMPTION==A3){
		getfem::asm_stiffness_matrix_for_laplacian(Av,mimv,mf_Cv,mf_coefv,Area); 	}
	else{
		getfem::asm_stiffness_matrix_for_laplacian(Av,mimv,mf_Cv,mf_coefv,Area);
		/*do nothing! A1 and A2 only tissue */ /* add something non-zero*/; }
}
	
	gmm::add(Av,  
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_t, dof_v), 
					gmm::sub_interval(dof_t, dof_v))); 


	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
	#endif
if (VERSION==PRIMAL){

	if(ASSUMPTION==A0){
	cout <<"---------------------------------------------------------"<<endl;
	cout<<"WARNING!! Primal reference problem is in developing phase!"<<std::endl;
	cout <<" ... U0t and U0v uncoupled! Do not trust this solution!"<<endl;
	cout <<"---------------------------------------------------------"<<endl;

	/* For exchange terms on Gamma, I must use interpolation matrix between Omega and Sigma.
	 * In fact, in the mesht, are stored only the elements, not their faces.
	 * i.e. elements from 1 to 10.000, and faces from 0 to 3 of each element. (there is not a unic ID for a face)
	 * Therefore, when i build Gamma to be the outer faces of Sigma, it takes the faces of the elements stored in Sigma. 
	 * If i buildt Gamma to be the outer faces of Omega (beware, you should take care of the exterior faces of the cube),
	 * this region would contain exactly the same faces, but with different ID.
	 * i.e. face 2 of element 120 in Omega == face 1 of element 6300 in Sigma, but GetFEM DOESN'T KNOW THAT!

	 * SOLUTION: I build region Gamma+1 to be the outer faces of Sigma AND outer faces of Omega.
	 * I build the mesh_fem of Gamma+1, and build the mass matrix on Gamma (M_Gamma).
	 * I assemble the interpolation matrixes from Gamma to Sigma and from Gamma to Omega (INT_Gamma_Sigma and INT_Gamma_Omega).
	 * (u_omega, v_sigma)_Gamma = INT_Gamma_Omega* (INT_Gamma_Sigma * M_Gamma)^T = Btv = (Bvt)^T 
	 */
		//The mass terms on Gamma are easy to compute
		getfem::asm_mass_matrix_param (Btt, mimt, mf_Ct_Omega,mf_Ct_Omega, mf_coeft, Kt, descr_transp.GAMMA);
		getfem::asm_mass_matrix_param (Bvv, mimt, mf_Ct_Sigma,mf_Ct_Sigma, mf_coeft, Kt, descr_transp.GAMMA+1);
	}
	else if(ASSUMPTION==A1){
		if( (gmm::mat_nrows(MBAR)==dof_v) && (gmm::mat_ncols(MBAR)==dof_t) 
		&&  (gmm::mat_nrows(MLIN)==dof_v) && (gmm::mat_ncols(MLIN)==dof_t)
		&&  (gmm::mat_maxnorm(MBAR)>=1e-16) && (gmm::mat_maxnorm(MLIN)>=1e-16) ) //Mbar and Mlin are already defined
			{gmm::copy(MBAR, Mbar); gmm::copy(MLIN, Mlin); }
		else{ //Mbar and Mlin are not already defined
			asm_exchange_aux_mat_transp(Mbar, Mlin, 
				mimv, mf_Ct_Omega,       mf_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);
			gmm::resize(MLIN, dof_transp.Cv(), dof_transp.Ct());
			gmm::resize(MBAR, dof_transp.Cv(), dof_transp.Ct());
			gmm::copy(Mbar,MBAR); gmm::copy(Mlin,MLIN);
		}
	}
	else{
		if( (gmm::mat_nrows(MBAR)==dof_v) && (gmm::mat_ncols(MBAR)==dof_t) 
		&&  (gmm::mat_nrows(MLIN)==dof_v) && (gmm::mat_ncols(MLIN)==dof_t)
		&&  (gmm::mat_maxnorm(MBAR)>=1e-16) && (gmm::mat_maxnorm(MLIN)>=1e-16) ) //Mbar and Mlin are already defined
			{gmm::copy(MBAR, Mbar); gmm::copy(MLIN, Mlin); }
		else{ //Mbar and Mlin are not already defined
			asm_exchange_aux_mat_transp(Mbar, Mlin, 
				mimv, mf_Ct,       mf_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);
			gmm::resize(MLIN, dof_transp.Cv(), dof_transp.Ct());
			gmm::resize(MBAR, dof_transp.Cv(), dof_transp.Ct());
			gmm::copy(Mbar,MBAR); gmm::copy(Mlin,MLIN);
		}
	}

	bool NEWFORM = 1;//PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");	
	if(ASSUMPTION != A0){									   
		asm_exchange_mat(Btt, Btv, Bvt, Bvv,
			mimv, mf_Cv, mf_coefv, Mbar, Mlin, Perimeter, NEWFORM);}
	// Copying Bvt
	gmm::add(scaled(Bvt,-1),								
			  gmm::sub_matrix(AM_transp, 
			  		gmm::sub_interval(dof_t, dof_v),
					gmm::sub_interval(0, dof_t)));
	// Copying Bvv
	gmm::add(Bvv,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_t, dof_v), 
					gmm::sub_interval(dof_t, dof_v))); 

	if(ASSUMPTION==A1 || ASSUMPTION==A2){
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling mass matrices on Gamma..." << endl;
	#endif
		// Assemble Mtt, mass matrix on tissue
		if(ASSUMPTION==A1){
			getfem::asm_mass_matrix_param (Mtt, mimt, mf_Ct_Omega, mf_coeft, Kt, descr_transp.GAMMA);}
		else{
			getfem::asm_mass_matrix_param (Mtt, mimt, mf_Ct,       mf_coeft, Kt, descr_transp.GAMMA);}

		// Assemble Mtv, exchange matrix on gamma
		getfem::asm_mass_matrix_param (Mvv, mimv, mf_Cv,mf_coefv, Perimeter);
		gmm::mult(gmm::transposed(Mbar), Mvv, Mtv); // M1 * M2 ---> M3	
	
		// Copying Mtt
		gmm::add(Mtt,			 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t), 
					gmm::sub_interval(0, dof_t))); 
		// Copying Mtv

		gmm::add(scaled(Mtv,-1),
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t),
					gmm::sub_interval(dof_t, dof_v))); 

	} else { //ASSUMPTION==A3 || ASSUMPTION==A0
	// Copying Btt
		gmm::add(Btt,			 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t), 
					gmm::sub_interval(0, dof_t))); 

		// Copying Btv

		gmm::add(scaled(Btv,-1),
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t),
					gmm::sub_interval(dof_t, dof_v))); 

	}
}else if(VERSION==DUAL){
	if(ASSUMPTION==A0){
		getfem::asm_mass_matrix_param (Bvv, mimt, mf_Ct_Sigma, mf_coeft, Kt, descr_transp.GAMMA);
		gmm::add(Bvv,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_t, dof_v), 
					gmm::sub_interval(dof_t, dof_v))); 

	vector_type ones(dof_transp.Ct_Sigma());
	constant_value=1;
	interpolation_function(mf_Ct_Sigma, ones, constant_function );  


	}
	if(ASSUMPTION==A1){
		getfem::asm_mass_matrix_param (Btt, mimt, mf_Ct_Omega, mf_coeft, Kt, descr_transp.GAMMA);
		gmm::add(Btt,			 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t), 
					gmm::sub_interval(0, dof_t))); 
	}
	if(ASSUMPTION==A2){
		getfem::asm_mass_matrix_param (Btt, mimt, mf_Ct, mf_coeft, Kt, descr_transp.GAMMA);
		gmm::scale(Btt, 0.5);
		gmm::add(Btt,			 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t), 
					gmm::sub_interval(0, dof_t))); 
	}
	if(ASSUMPTION==A3){
		if( (gmm::mat_nrows(MBAR)==dof_v) && (gmm::mat_ncols(MBAR)==dof_t) 
		&&  (gmm::mat_nrows(MLIN)==dof_v) && (gmm::mat_ncols(MLIN)==dof_t)
		&&  (gmm::mat_maxnorm(MBAR)>=1e-16) && (gmm::mat_maxnorm(MLIN)>=1e-16) ) //Mbar and Mlin are already defined
			{gmm::copy(MBAR, Mbar); gmm::copy(MLIN, Mlin); }
		else{ //Mbar and Mlin are not already defined
			asm_exchange_aux_mat_transp(Mbar, Mlin, 
				mimv, mf_Ct,       mf_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);
			gmm::resize(MLIN, dof_transp.Cv(), dof_transp.Ct());
			gmm::resize(MBAR, dof_transp.Cv(), dof_transp.Ct());
			gmm::copy(Mbar,MBAR); gmm::copy(Mlin,MLIN);
		}
			bool NEWFORM = 1;//PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");	
			asm_exchange_mat(Btt, Btv, Bvt, Bvv,
				mimv, mf_Cv, mf_coefv, Mbar, Mlin, Perimeter, NEWFORM);
		gmm::add(Btt,			 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_t), 
					gmm::sub_interval(0, dof_t))); 			
		gmm::add(Bvv,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_t, dof_v), 
					gmm::sub_interval(dof_t, dof_v))); 
	}

}
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling source terms ..." << endl;
	#endif
if(VERSION==PRIMAL){
	// Update potential paramters in F and G

	kk=k;
	CC = PARAM.real_value("Cv", "value of concentration in the vessel");
	RR = param.R(0);

	//// Assemble F: source term in tissue
	vector_type F(dof.coeft());
	//Vector F is buildt with f_function via muparser
	expr = PARAM.string_value("F_EXPR", "Expression for source term in tissue (muparser)");
	interpolation_function(mf_coeft, F, f_function ); 
	//Vector F is assembled on Ct_Omega or on whole Ct?	
	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		getfem::asm_source_term(Ft, mimt, mf_Ct_Omega, mf_coeft, F);}
	else{
		getfem::asm_source_term(Ft, mimt, mf_Ct,       mf_coeft, F);}	
	

	//// Assemble G: source term in vessel
	if(ASSUMPTION==A0){	// G is defined on Sigma
		vector_type G(dof.coeft());
		//Vector G is buildt with g_function via muparser
		expr = PARAM.string_value("G_EXPR", "Expression for source term in vessel (muparser)");
		interpolation_function(mf_coeft, G, g_function );  
		//Assemble source term on Sigma
		getfem::asm_source_term(Fv, mimt, mf_Ct_Sigma, mf_coeft, G);
	vtk_export exp_Fv_A0(descr_transp.OUTPUT+"Fv_A0.vtk");
	exp_Fv_A0.exporting(mf_Ct_Sigma);
	exp_Fv_A0.write_mesh();
	exp_Fv_A0.write_point_data(mf_Ct_Sigma, Fv, "Fv_A0");
		}
	else{			// G is defined on Lambda

		vector_type G(dof_transp.Cv());gmm::clear(G);
		expr = PARAM.string_value("G_EXPR", "Expression for source term in vessel (muparser)");
		bool G_CONSTANT = PARAM.int_value("G_CONSTANT", "flag for using constant source term g in vessel");


		// If G is constant in the cross section, i just interpolate the g_function in Lambda
		if(G_CONSTANT){
			interpolation_function(mf_Cv, G, g_function );	
		}
		else{	//If G is not constant, I must build Mbarbar
			vector_type Gt(dof_transp.Ct()); gmm::clear(Gt);	
			interpolation_function(mf_Ct, Gt, g_function ); 		
		
	vtk_export exp_gt(descr_transp.OUTPUT+"Gt_A3.vtk");
	exp_gt.exporting(mf_Ct);
	exp_gt.write_mesh();
	exp_gt.write_point_data(mf_Ct, Gt, "Gt_A3");	
			// Build Mbarbar 
			if(gmm::mat_nrows(MBARBAR)==dof_transp.Cv()
			&& gmm::mat_ncols(MBARBAR)==dof_transp.Ct() 
			&& gmm::mat_maxnorm(MBARBAR)>=1e-16) //Mbarbar is already defined
				{gmm::copy(MBARBAR, Mbarbar); }
			else { 				     //Mbarbar is not already defined

				gmm::resize(MBARBAR, dof_transp.Cv(), dof_transp.Ct());	
				gmm::clear(MBARBAR);
				asm_exchange_aux_mat_bar_bar(Mbarbar,mimt, mf_Ct, mf_Cv,mf_coefv, param.R(), descr.NInt, descr_transp.NIntA, nb_branches);
				gmm::copy(Mbarbar, MBARBAR);  
				}
			gmm::mult(Mbarbar, Gt, G);		
		}
		
		//Build source term
		vector_type Area_Cv(dof_transp.Cv());gmm::clear(Area_Cv);
		getfem::interpolation(mf_coefv, mf_Cv, Area, Area_Cv);
		gmm::vscale(Area_Cv, G); //G=G*pi*R^2
		getfem::asm_source_term(Fv, mimv, mf_Cv, mf_Cv, G);

	} 


	// Add source term to RHS
	gmm::add(Ft, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(0,dof_t)));

	gmm::add(Fv, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(dof_t, dof_v)));
}else if (VERSION==DUAL){


	string functional= PARAM.string_value("FUNCTIONAL", "descriptor for the functional J of the error you wanto to estimate");

	if(functional=="MEAN_VALUE"){
	// Aux tissue source vector
	vector_type Ft(dof.coeft());
	// Aux vessels source vector
	vector_type Fv(dof.coefv());

		expr = "1";
		interpolation_function(mf_coeft, Ft, g_function );  
		interpolation_function(mf_coefv, Fv, g_function );  

	if(ASSUMPTION==A0){
		getfem::asm_source_term(Jv,mimt, mf_Ct_Sigma, mf_coeft, Ft);


}else 	if(ASSUMPTION==A1){
		getfem::asm_source_term(Jt,mimt, mf_Ct_Omega, mf_coeft, Ft);
}else	if(ASSUMPTION==A2){
		getfem::asm_source_term(Jt,mimt, mf_Ct,       mf_coeft, Ft);
}else 	if(ASSUMPTION==A3){
		getfem::asm_source_term(Jt,mimt, mf_Ct,       mf_coeft, Ft);
		getfem::asm_source_term(Jv,mimv, mf_Cv, mf_coefv, Fv );
}
	}
	else if (functional=="GAMMA_FLUX"){

	/*
vector_type ones(3,1);

	if(ASSUMPTION==A0){
		getfem::asm_homogeneous_normal_source_term (Jv, mimt, mf_Ct_Sigma, ones, descr_transp.GAMMA);
}else 	if(ASSUMPTION==A1){
		getfem::asm_homogeneous_normal_source_term (Jt, mimt, mf_Ct_Omega, ones, descr_transp.GAMMA);
}else	if(ASSUMPTION==A2){
		getfem::asm_homogeneous_normal_source_term (Jt, mimt, mf_Ct      , ones, descr_transp.GAMMA);
		gmm::scale(Jt, 0.5);
}else 	if(ASSUMPTION==A3){
		getfem::asm_homogeneous_normal_source_term (Jt, mimt, mf_Ct      , ones, descr_transp.GAMMA);
		gmm::scale(Jt, 0.5);
		for (size_type bc=0; bc < BCv_transp.size(); bc++) { 
			getfem::asm_homogeneous_normal_source_term (Jv, mimv, mf_Cv, ones, BCv_transp[bc].rg);}
} 
*/



	if(ASSUMPTION==A0){
		getfem::asm_flux(Jv, mimt, mf_Ct_Sigma, descr_transp.GAMMA);
}else 	if(ASSUMPTION==A1){
		getfem::asm_flux (Jt, mimt, mf_Ct_Omega,  descr_transp.GAMMA);
}else	if(ASSUMPTION==A2){
		getfem::asm_flux (Jt, mimt, mf_Ct      ,  descr_transp.GAMMA);
		gmm::scale(Jt, 0.5);
}else 	if(ASSUMPTION==A3){
		getfem::asm_flux (Jt, mimt, mf_Ct      ,  descr_transp.GAMMA);
		gmm::scale(Jt, 0.5);
		getfem::asm_flux (Jv, mimv, mf_Cv);
} 


/*
	if(functional=="L2NORM")
	// Assemble F: source term in tissue
	sparse_matrix_type JJt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(JJt);
		
	generic_assembly L2normt;
	L2normt.push_mi(mimt);
	L2normt.push_mf(mf_Ct);
	L2normt.push_mat(JJt);
	L2normt.set("V(#1,#1)+=comp(Base(#1).Base(#1))(:,:)");
	L2normt.assembly();
	
	//gmm::add(gmm::mat_trace(JJt),Jt);
	
	for (int i=0; i>dof_transp.Ct(); i++){
		Jt[i]=JJt(i,i);
		Jt[i]=sqrt(Jt[i]);
	};
	// Assemble G: source term in vessel
	//! /todo Extend to variable source term: use muParser and build cross section average (operator Mbarbar)
	sparse_matrix_type JJv(dof_v, dof_v);gmm::clear(JJv);
		
	generic_assembly L2normv;
	L2normv.push_mi(mimv);
	L2normv.push_mf(mf_Cv);
	L2normv.push_mat(JJv);
	L2normv.set("V(#1,#1)+=comp(Base(#1).Base(#1))(:,:)");
	L2normv.assembly();
	
	//gmm::add(gmm::mat_trace(JJv),Jv);
	for (int i=0; i>dof_v; i++){
		Jv[i]=JJv(i,i);
		Jv[i]=sqrt(Jv[i]);
		
	};
*/	


	}

	//Add to FM_temp

	if(ASSUMPTION==A0){	
			gmm::add(Jv, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(dof_t, dof_v)));
}else 	if(ASSUMPTION==A1){
			gmm::add(Jt, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(0,dof_t)));;
}else	if(ASSUMPTION==A2){
			gmm::add(Jt, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(0,dof_t)));;
}else 	if(ASSUMPTION==A3){
			gmm::add(Jt, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(0,dof_t)));

	gmm::add(Jv, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(dof_t, dof_v)));
	} 

}


	// De-allocate memory
	gmm::clear(At);	   gmm::clear(Av); 
	gmm::clear(Mtt);   gmm::clear(Mtv); 
	gmm::clear(Mvt);   gmm::clear(Mvv); 	    
	gmm::clear(Mbar);  gmm::clear(Mlin);   gmm::clear(Mbarbar);  
	gmm::clear(Btt);   gmm::clear(Btv);
	gmm::clear(Bvt);   gmm::clear(Bvv);
	gmm::clear(Ft);   gmm::clear(Fv);
	gmm::clear(Jt);   gmm::clear(Jv);

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling boundary conditions ..." << endl;
	#endif	
	gmm::resize(AM_temp, dof_tot,
			     dof_tot);	gmm::clear(AM_temp);
	gmm::resize(FM_temp, dof_tot); 	gmm::clear(FM_temp);
	gmm::copy(AM_transp,AM_temp);
	gmm::copy(FM_transp,FM_temp);

	//Boundary conditions: 
	//On vessel, we have homogeneous neumann conditions: we simply don't implement the boundary term due to diffusion.

if(VERSION==PRIMAL){
	//Homogeneous Dirichlet on all the faces of tissue	
	vector_type Dirichlet_null(dof_t,0);
	for (size_type bc=0; bc < BCt_transp.size(); ++bc) {

		if(ASSUMPTION==A1 || ASSUMPTION==A0){
		getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct_Omega, BCt_transp[bc].rg, Dirichlet_null);	}
		else{
		getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct,       BCt_transp[bc].rg, Dirichlet_null);	}
	} 
	
} else if (VERSION==DUAL){
	vector_type Dirichlet_null(dof_t,0);
	for (size_type bc=0; bc < BCt_transp.size(); ++bc) {
	if(ASSUMPTION==A1){
	getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct_Omega, BCt_transp[bc].rg, Dirichlet_null);}
	else if(ASSUMPTION!=A0){
	getfem::assembling_Dirichlet_condition(AM_temp, FM_temp, mf_Ct,       BCt_transp[bc].rg, Dirichlet_null);}
			
	} 
}


};




	//! Solve problem U1
	bool transport3d1d::solve_model(const size_type ASSUMPTION,const  size_type VERSION){

  	#ifdef M3D1D_VERBOSE_
	cout << "Solving the monolithic system ... " << endl;
	#endif
	double time = gmm::uclock_sec();
	gmm::clear(UM_transp);
	gmm::clean(AM_temp, 1E-12);
	gmm::clean(FM_temp, 1E-12);


	// define dofs
	size_type dof_t;
	size_type dof_v;

	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		dof_t  = dof_transp.Ct_Omega();	}
	else{
		dof_t  = dof_transp.Ct();	}
	
 	if(ASSUMPTION==A0){
		dof_v  = dof_transp.Ct_Sigma();	}
	else{
		dof_v  = dof_transp.Cv();	}

	size_type dof_tot= dof_t+dof_v;

	// Solve the system on AM_temp, UM_transp, FM_temp
	bool solved = solver(dof_t,dof_v);
	if (!solved) return false;
	//export solution
	#ifdef M3D1D_VERBOSE_
	std::cout<<"solved! going to export..."<<std::endl;
	#endif		



	switch(VERSION){

		case PRIMAL:{
		switch(ASSUMPTION){

			case A0: 
			{	gmm::resize(U0, dof_tot);	gmm::clear(U0);
				gmm::copy(UM_transp, U0);
				cout << endl<<"... time to solve U0: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A1: 
			{	gmm::resize(U1, dof_tot);	gmm::clear(U1);
				gmm::copy(UM_transp, U1);
				cout << endl<<"... time to solve U1: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A2: 
			{	gmm::resize(U2, dof_tot);	gmm::clear(U2);
				gmm::copy(UM_transp, U2);	
				cout << endl<<"... time to solve U2: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A3: 
			{	gmm::resize(U3, dof_tot);	gmm::clear(U3);
				gmm::copy(UM_transp, U3);	
				cout << endl<<"... time to solve U3: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
		 }
				break;} 

		case DUAL:{
		switch(ASSUMPTION){

			case A0: 
			{	gmm::resize(Z0, dof_tot);	gmm::clear(Z0);
				gmm::copy(UM_transp, Z0);	
				cout << endl<<"... time to solve Z0: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A1: 
			{	gmm::resize(Z1, dof_tot);	gmm::clear(Z1);
				gmm::copy(UM_transp, Z1);	
				cout << endl<<"... time to solve Z1: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A2: 
			{	gmm::resize(Z2, dof_tot);	gmm::clear(Z2);
				gmm::copy(UM_transp, Z2);	
				cout << endl<<"... time to solve Z2: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
			case A3: 
			{	gmm::resize(Z3, dof_tot);	gmm::clear(Z3);
				gmm::copy(UM_transp, Z3);	
				cout << endl<<"... time to solve Z3: " << gmm::uclock_sec() - time << " seconds\n";
				break;
			}
		 }
				break;} 

		}
	return true;

	};
	//! Export problem U1
	void transport3d1d::export_model(const size_type ASSUMPTION,const  size_type VERSION, const string & suff ){

	#ifdef M3D1D_VERBOSE_
	cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
	#endif
	#ifdef M3D1D_VERBOSE_
	cout << "  Saving the results from the monolithic unknown vector ... " << endl;
	#endif
	// define dofs
	size_type dof_t;
	size_type dof_v;

	if(ASSUMPTION==A1 || ASSUMPTION==A0){
		dof_t  = dof_transp.Ct_Omega();	}
	else{
		dof_t  = dof_transp.Ct();	}

 	if(ASSUMPTION==A0){
		dof_v  = dof_transp.Ct_Sigma();	}
	else{
		dof_v  = dof_transp.Cv();	}

	size_type dof_tot= dof_t+dof_v;
	
	// Array of unknown dof of the interstitial velocity
	vector_type Ct(dof_t); 

	// Array of unknown dof of the network velocity
	vector_type Cv(dof_v); 

	//Copy solution


	switch(VERSION){

		case PRIMAL:{
		switch(ASSUMPTION){

			case A0: 
			{	gmm::copy(gmm::sub_vector(U0, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(U0, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A1: 
			{	gmm::copy(gmm::sub_vector(U1, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(U1, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A2: 
			{	gmm::copy(gmm::sub_vector(U2, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(U2, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A3: 
			{	gmm::copy(gmm::sub_vector(U3, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(U3, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
		}
				break;} 

		case DUAL:{
		switch(ASSUMPTION){

			case A0: 
			{	gmm::copy(gmm::sub_vector(Z0, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(Z0, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A1: 
			{	gmm::copy(gmm::sub_vector(Z1, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(Z1, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A2: 
			{	gmm::copy(gmm::sub_vector(Z2, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(Z2, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
			case A3: 
			{	gmm::copy(gmm::sub_vector(Z3, 
					  gmm::sub_interval(0, dof_t)), Ct);
				gmm::copy(gmm::sub_vector(Z3, 
					  gmm::sub_interval(dof_t, dof_v)), Cv);
				break;
			}
		}
				break;} 

	}


	string assumption, version;	
	switch(VERSION){
		case PRIMAL:
			version="U";
			break;
		case DUAL:
			version="Z";
			break;
	}
	
	switch(ASSUMPTION){
		case A0:
			assumption="0";
			break;
		case A1:
			assumption="1";
			break;
		case A2:
			assumption="2";
			break;
		case A3:
			assumption="3";
			break;
	}

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting reduced solution in tissue: "<<version<<assumption << endl;
	#endif
	if(ASSUMPTION==A1||ASSUMPTION==A0){
		vtk_export exp_Ct(descr_transp.OUTPUT+version+"t"+assumption+suff+".vtk");
		exp_Ct.exporting(mf_Ct_Omega);
		exp_Ct.write_mesh();
		exp_Ct.write_point_data(mf_Ct_Omega, Ct, version+"t"+assumption);
	}else{
		vtk_export exp_Ct(descr_transp.OUTPUT+version+"t"+assumption+suff+".vtk");
		exp_Ct.exporting(mf_Ct);
		exp_Ct.write_mesh();
		exp_Ct.write_point_data(mf_Ct, Ct, version+"t"+assumption);
	}


	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting reduced solution in vessel: "<<version<<assumption << endl;
	#endif

	if(ASSUMPTION==A0){
		vtk_export exp_Cv(descr_transp.OUTPUT+version+"v"+assumption+suff+".vtk");
		exp_Cv.exporting(mf_Ct_Sigma);
		exp_Cv.write_mesh();
		exp_Cv.write_point_data(mf_Ct_Sigma, Cv, version+"v"+assumption);
		
	}else{
		vtk_export exp_Cv(descr_transp.OUTPUT+version+"v"+assumption+suff+".vtk");
		exp_Cv.exporting(mf_Cv);
		exp_Cv.write_mesh();
		exp_Cv.write_point_data(mf_Cv, Cv, version+"v"+assumption);
	}

	if(ASSUMPTION==A0 && VERSION==PRIMAL){
		vector_type u0(dof_transp.Ct_Sigma()); 	
		vector_type u0_1D(dof_transp.Cv());
		vector_type u0_3D(dof_transp.Ct());

		gmm::copy(gmm::sub_vector(U0, 
				  gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Ct_Sigma())), u0);

		getfem::interpolation(mf_Ct_Sigma, mf_Ct, u0, u0_3D);

		gmm::mult(MBARBAR, u0_3D,u0_1D);		
		vtk_export exp_Cv2(descr_transp.OUTPUT+"Uv0_1D.vtk");
		exp_Cv2.exporting(mf_Cv);
		exp_Cv2.write_mesh();
		exp_Cv2.write_point_data(mf_Cv, u0_1D, "Uv0_1D.vtk");
	}

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl<<endl; 
	#endif
};


	/////////////////////////////////////////////
	// Compute model errors -- final results



	//! Compute model error of A1	
	void transport3d1d::compute_model_error_1(const size_type DUALMODEL){
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  model error for assumption A1 ..." << endl;
	#endif

//////////////////////////////////////////////////////// Estimators
	// Estimators for d1
	vector_type d1_grad_ref(dof_transp.Ct()); gmm::clear(d1_grad_ref); //d1_grad_ref = (gradU, gradPhi)Sigma
	vector_type d1_grad_mean(dof_transp.Ct()); gmm::clear(d1_grad_mean); //d1_grad_mean = (grad barbarU, gradPhi)Sigma	
	vector_type d1_grad_residual(dof_transp.Ct()); gmm::clear(d1_grad_residual); //d1_grad_residual = (grad (i -barbar)U, gradPhi)Sigma
	vector_type d1_grad(dof_transp.Ct()); gmm::clear(d1_grad); //d1_grad = (grad (i -barbar)U, gradZ)Sigma

	vector_type d1_flux_ref(dof_transp.Ct()); gmm::clear(d1_flux_ref); //d1_flux_ref = (k U, Phi)Gamma
	vector_type d1_flux_mean(dof_transp.Ct()); gmm::clear(d1_flux_mean); //d1_flux_mean = (k barU, Phi)Gamma	
	vector_type d1_flux_residual(dof_transp.Ct()); gmm::clear(d1_flux_residual); //d1_flux_residual = (k (i-bar)U, Phi)Gamma
	vector_type d1_flux(dof_transp.Ct()); gmm::clear(d1_flux); //d1_flux_ = (k (i-bar)U, Z)Gamma

	vector_type d1_residual(dof_transp.Ct()); gmm::clear(d1_residual); //d1_residual = d1_flux_residual + d1_grad_residual
	vector_type d1(dof_transp.Ct()); gmm::clear(d1); //d1 = d1_grad + d1_flux

	// Estimators for l1
	vector_type l1_g_ref(dof_transp.Ct()); gmm::clear(l1_g_ref); //l1_g_ref = (g, Phi)Sigma
	vector_type l1_g_mean(dof_transp.Ct()); gmm::clear(l1_g_mean); //l1_g_mean = ( barbarg, Phi)Sigma	
	vector_type l1_g_residual(dof_transp.Ct()); gmm::clear(l1_g_residual); //l1_g_residual = ( (i -barbar)g, Phi)Sigma
	vector_type l1_g(dof_transp.Ct()); gmm::clear(l1_g); //l1_g = ((i -barbar)g, Z)Sigma

	vector_type l1_flux_ref(dof_transp.Ct()); gmm::clear(l1_flux_ref); //l1_flux_ref = (k u, Phi)Gamma
	vector_type l1_flux_mean(dof_transp.Ct()); gmm::clear(l1_flux_mean); //l1_flux_mean = (k baru, Phi)Gamma	
	vector_type l1_flux_residual(dof_transp.Ct()); gmm::clear(l1_flux_residual); //l1_flux_residual = (k (i-bar)u, Phi)Gamma
	vector_type l1_flux(dof_transp.Ct()); gmm::clear(l1_flux); //l1_flux = (k (i-bar)u, Z)Gamma

	vector_type l1_residual(dof_transp.Ct()); gmm::clear(d1_residual); //l1_residual = l1_flux_residual + l1_g_residual
	vector_type l1(dof_transp.Ct()); gmm::clear(d1); //l1 = l1_g + l1_flux	

	// Estimators for eta1
	vector_type eta1_residual(dof_transp.Ct()); gmm::clear(eta1_residual); //eta1_residual = l1_residual - d1_residual 
	vector_type eta1(dof_transp.Ct()); gmm::clear(eta1); //eta1 = l1 - d1	



//////////////////////////////////////////////////////// weight Z (dual solution) 

	vector_type Z(dof_transp.Ct()); gmm::clear(Z); // this will store the weights Z
	// reference or reduced solution?
	string dual="";
	if(DUALMODEL == REDUCED){
		vector_type Z_Lambda(dof_transp.Cv()); gmm::clear(Z_Lambda);
		vector_type Z_Sigma(dof_transp.Ct_Sigma()); gmm::clear(Z_Sigma);
		gmm::copy(gmm::sub_vector(Z3, 
			gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), Z_Lambda);
		getfem::interpolation(mf_Cv, mf_Ct_Sigma,Z_Lambda, Z_Sigma, 2);
		getfem::interpolation(mf_Ct_Sigma, mf_Ct, Z_Sigma, Z, 0);

		dual = "red";
	} 

	if(DUALMODEL == REFERENCE ){

		vector_type Z_Sigma(dof_transp.Ct_Sigma()); gmm::clear(Z_Sigma);
		gmm::copy(gmm::sub_vector(Z0, 
			gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Ct_Sigma())), Z_Sigma);
		getfem::interpolation(mf_Ct_Sigma, mf_Ct, Z_Sigma, Z, 0);
	} 	


//////////////////////////////////////////////////////// Primal solutions, u and U

	vector_type u(dof_transp.Ct()); gmm::clear(u); // u = U3_t (def su Omega)
	vector_type U(dof_transp.Cv()); gmm::clear(U); // U = U3_v (def su lambda)	
	vector_type U_Sigma(dof_transp.Ct_Sigma()); gmm::clear(U_Sigma); // U_Sigma = U (interp su Sigma, extra = 2)
	vector_type U_Omega(dof_transp.Ct()); gmm::clear(U_Omega); // U_Omega = U (interp su Omega, extra = 0)

	gmm::copy(gmm::sub_vector(U3, 
				  gmm::sub_interval(0, dof_transp.Ct())), u);
	gmm::copy(gmm::sub_vector(U3, 
				  gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv())), U);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, U, U_Sigma, 2);
	getfem::interpolation(mf_Ct_Sigma, mf_Ct, U_Sigma, U_Omega, 0);


	

//////////////////////////////////////////////////////// parameters
	scalar_type k = PARAM.real_value("kk", "value of permeability for reduced problem");
	// Permeability
	vector_type Kv(dof.coefv());
	vector_type Kt(dof.coeft());
	constant_value=k;
	interpolation_function(mf_coefv, Kv, constant_function ); 
	interpolation_function(mf_coeft, Kt, constant_function ); 
	//Perimeter = 2*pi*R*k 
	vector_type Perimeter(dof.coefv());
	gmm::copy(param.R(), Perimeter);
	gmm::scale(Perimeter, 2*pi*k);
	//Area = pi*R^2
	vector_type Area(dof.coefv());
	gmm::copy(param.R(), Area);
	gmm::vscale(param.R(), Area);
	gmm::scale(Area, pi); 


//////////////////////////////////////////////////////// l1_flux

//////////////////// Build l1_flux_ref = (k u, Phi)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l1_flux_ref..." << endl;
	#endif
	sparse_matrix_type M_Omega (dof_transp.Ct(),dof_transp.Ct()); gmm::clear(M_Omega);
	getfem::asm_mass_matrix_param(M_Omega, mimt, mf_Ct, mf_Ct, mf_coeft, Kt, descr_transp.GAMMA);
	gmm::mult(M_Omega, u, l1_flux_ref);	

//////////////////// Build l1_flux_mean = (k baru, Phi)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l1_flux_mean..." << endl;
	#endif
	vector_type u_mean(dof_transp.Cv()); gmm::clear(u_mean);
	vector_type u_sigma(dof_transp.Ct_Sigma()); gmm::clear(u_sigma);
	vector_type u_bar(dof_transp.Ct()); gmm::clear(u_bar);
	gmm::mult(MBAR, u, u_mean);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, u_mean, u_sigma, 2);
	getfem::interpolation(mf_Ct_Sigma, mf_Ct, u_sigma, u_bar, 0);
	gmm::mult(M_Omega, u_bar, l1_flux_mean);

//////////////////// Build l1_flux_residual = (k (i-bar)u, Phi)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l1_flux_residual..." << endl;
	#endif
	gmm::add(l1_flux_ref, gmm::scaled(l1_flux_mean, -1.0), l1_flux_residual);

//////////////////// Build l1_flux = (k (i-bar)u, Z)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l1_flux..." << endl;
	#endif
	gmm::copy(l1_flux_residual, l1_flux);
	gmm::vscale(Z, l1_flux);

//////////////////////////////////////////////////////// l1_g

		vector_type g(dof_transp.Ct()); gmm::clear(g);
		expr = PARAM.string_value("G_EXPR", "Expression for source term in vessel (muparser)");
		bool G_CONSTANT = PARAM.int_value("G_CONSTANT", "flag for using constant source term g in vessel");
		interpolation_function(mf_Ct, g, g_function ); 
//////////////////// Build l1_g_ref = (g, Phi)Sigma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l1_g_ref..." << endl;
	#endif
		if(G_CONSTANT){
			//do nothing! l1_g=0
		}
		else{
	sparse_matrix_type M_g (dof_transp.Ct(),dof_transp.Ct()); gmm::clear(M_g);
	getfem::asm_mass_matrix(M_g, mimt, mf_Ct, mf_Ct, descr_transp.SIGMA);
	gmm::mult(M_g, g, l1_g_ref);
	}

//////////////////// Build l1_g_mean = ( barbarg, Phi)Sigma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing l1_g_mean ..." << endl;
	#endif
		if(G_CONSTANT){
			//do nothing! l1_g=0
		}
		else{
	sparse_matrix_type M_g (dof_transp.Ct(),dof_transp.Ct()); gmm::clear(M_g);
	getfem::asm_mass_matrix(M_g, mimt, mf_Ct, mf_Ct, descr_transp.SIGMA);
	vector_type g_mean(dof_transp.Cv()); gmm::clear(g_mean);
	vector_type g_sigma(dof_transp.Ct_Sigma()); gmm::clear(g_sigma);
	vector_type g_bar(dof_transp.Ct()); gmm::clear(g_bar);
	gmm::mult(MBARBAR, g, g_mean);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, g_mean, g_sigma, 2);
	getfem::interpolation(mf_Ct_Sigma, mf_Ct, g_sigma, g_bar, 0);	
	gmm::mult(M_g, g_bar, l1_g_mean);
	}

//////////////////// Build l1_g_residual = ( (i -barbar)g, Phi)Sigma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l1_g_residual..." << endl;
	#endif
		if(G_CONSTANT){
			//do nothing! l1_g=0
		}
		else{
	gmm::add(l1_g_ref, gmm::scaled(l1_g_mean, -1.0), l1_g_residual);
	}

//////////////////// Build l1_g = ((i -barbar)g, Z)Sigma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l1_g..." << endl;
	#endif
		if(G_CONSTANT){
			//do nothing! l1_g=0
		}
		else{
	gmm::copy(l1_g_residual, l1_g);
	gmm::vscale(Z, l1_g);
	}
//////////////////////////////////////////////////////// l1

//////////////////// Build l1_residual = l1_flux_residual + l1_g_residual
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l1_residual..." << endl;
	#endif
	gmm::add(l1_flux_residual, l1_g_residual, l1_residual );

//////////////////// Build l1 = l1_g + l1_flux	
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l1..." << endl;
	#endif
	gmm::add(l1_flux, l1_g, l1 );

//////////////////////////////////////////////////////// d1_grad
			if(gmm::mat_nrows(MBARBAR)==dof_transp.Cv()
			&& gmm::mat_ncols(MBARBAR)==dof_transp.Ct() 
			&& gmm::mat_maxnorm(MBARBAR)>=1e-16) //Mbarbar is already defined
				{ /* ok! do nothing*/ }
			else { 				     //Mbarbar is not already defined: anyway, we don't need it, just resize it
				gmm::resize(MBARBAR, dof_transp.Cv(), dof_transp.Ct());	
			}

//////////////////// Build d1_grad_ref = (gradU, gradPhi)Sigma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing d1_grad_ref ..." << endl;
	#endif
	sparse_matrix_type A_Omega(dof_transp.Ct(),dof_transp.Ct()); gmm::clear(A_Omega);
	getfem::asm_stiffness_matrix_for_homogeneous_laplacian (A_Omega, mimt, mf_Ct, descr_transp.SIGMA);
	gmm::mult(A_Omega, U_Omega, d1_grad_ref);	

//////////////////// Build d1_grad_mean = (grad barbarU, gradPhi)Sigma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  d1_grad_mean..." << endl;
	#endif
	vector_type U_mean2(dof_transp.Cv()); gmm::clear(U_mean2);
	vector_type U_sigma2(dof_transp.Ct_Sigma()); gmm::clear(U_sigma2);
	vector_type U_barbar(dof_transp.Ct()); gmm::clear(U_barbar);
	gmm::mult(MBARBAR, U_Omega, U_mean2);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, U_mean2, U_sigma2, 2);
	getfem::interpolation(mf_Ct_Sigma, mf_Ct, U_sigma2, U_barbar, 0);

	gmm::mult(A_Omega, U_barbar, d1_grad_mean);

//////////////////// Build d1_grad_residual = (grad (i -barbar)U, gradPhi)Sigma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  d1_grad_residual..." << endl;
	#endif
	gmm::add(d1_grad_ref, gmm::scaled(d1_grad_mean, -1.0), d1_grad_residual);

//////////////////// Build d1_grad = (grad (i -barbar)U, gradZ)Sigma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  d1_grad..." << endl;
	#endif
	gmm::copy(d1_grad_residual, d1_grad);
	gmm::vscale(Z, d1_grad);
//////////////////////////////////////////////////////// d1_flux

//////////////////// Build d1_flux_ref = (k U, Phi)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  d1_flux_ref..." << endl;
	#endif
	gmm::mult(M_Omega, U_Omega, d1_flux_ref);	

//////////////////// Build d1_flux_mean = (k barU, Phi)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing d1_flux_mean ..." << endl;
	#endif
	vector_type U_mean(dof_transp.Cv()); gmm::clear(U_mean);
	vector_type U_sigma(dof_transp.Ct_Sigma()); gmm::clear(U_sigma);
	vector_type U_bar(dof_transp.Ct()); gmm::clear(U_bar);
	gmm::mult(MBAR, U_Omega, U_mean);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, U_mean, U_sigma, 2);
	getfem::interpolation(mf_Ct_Sigma, mf_Ct, U_sigma, U_bar, 0);
	gmm::mult(M_Omega, U_bar, d1_flux_mean);

//////////////////// Build d1_flux_residual = (k (i-bar)U, Phi)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing d1_flux_residual ..." << endl;
	#endif
	gmm::add(d1_flux_ref, gmm::scaled(d1_flux_mean, -1.0), d1_flux_residual);

//////////////////// Build d1_flux = (k (i-bar)U, Z)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing d1_flux ..." << endl;
	#endif
	gmm::copy(d1_flux_residual, d1_flux);
	gmm::vscale(Z, d1_flux);

//////////////////////////////////////////////////////// d1

//////////////////// Build d1_residual = d1_flux_residual + d1_grad_residual
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing d1_residual ..." << endl;
	#endif
	gmm::add(d1_flux_residual, d1_grad_residual, d1_residual );

//////////////////// Build d1 = d1_grad + d1_flux
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing d1 ..." << endl;
	#endif
	gmm::add(d1_flux, d1_grad, d1 );
//////////////////////////////////////////////////////// eta 1

//////////////////// Build eta1_residual = l1_residual - d1_residual
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing eta1_residual ..." << endl;
	#endif
	gmm::add(l1_residual, gmm::scaled(d1_residual, -1.0), eta1_residual );
	gmm::add(l1_residual, ETA_RESIDUAL);
//////////////////// Build eta1 = l1 - d1
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing eta1 ..." << endl;
	#endif
	gmm::add(l1, gmm::scaled(d1, -1.0), eta1 );
	gmm::add(l1, ETA);
/////////////////////////////////////////////////////// Global integrals
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  global integrals of d1 and l1..." << endl;
	#endif
	cout <<"l1_flux_ref = "<<gmm::vect_sp(Z, l1_flux_ref)<<endl;
	cout <<"l1_flux_mean = "<<gmm::vect_sp(Z, l1_flux_mean)<<endl;
	cout <<"l1_flux = "<<gmm::vect_sp(Z, l1_flux_residual)<<endl;
	cout <<"l1_g_ref = "<<gmm::vect_sp(Z, l1_g_ref)<<endl;
	cout <<"l1_g_mean = "<<gmm::vect_sp(Z, l1_g_mean)<<endl;
	cout <<"l1_g = "<<gmm::vect_sp(Z, l1_g_residual)<<endl;

	cout <<"d1_flux_ref = "<<gmm::vect_sp(Z, d1_flux_ref)<<endl;
	cout <<"d1_flux_mean = "<<gmm::vect_sp(Z, d1_flux_mean)<<endl;
	cout <<"d1_flux = "<<gmm::vect_sp(Z, d1_flux_residual)<<endl;
	cout <<"d1_grad_ref = "<<gmm::vect_sp(Z, d1_grad_ref)<<endl;
	cout <<"d1_grad_mean = "<<gmm::vect_sp(Z, d1_grad_mean)<<endl;
	cout <<"d1_grad = "<<gmm::vect_sp(Z, d1_grad_residual)<<endl;

	cout <<"l1 = "<<gmm::vect_sp(Z, l1_residual)<<endl;
	cout <<"d1 = "<<gmm::vect_sp(Z, d1_residual)<<endl;
	cout <<"eta1 = "<<gmm::vect_sp(Z, eta1_residual)<<endl;



////////////////////////////////////// E X P O R T //////////////////////////////////////
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting l1 and d1  ..." << endl;
	#endif
//////////////////////////////////////////////////////// l1_flux

//////////////////// Build l1_flux_ref
	vtk_export exp_l1_flux_ref(descr_transp.OUTPUT+"l1_flux_ref.vtk");
	exp_l1_flux_ref.exporting(mf_Ct);
	exp_l1_flux_ref.write_mesh();
	exp_l1_flux_ref.write_point_data(mf_Ct, l1_flux_ref, "l1_flux_ref");

//////////////////// Build l1_flux_mean
	vtk_export exp_l1_flux_mean(descr_transp.OUTPUT+"l1_flux_mean.vtk");
	exp_l1_flux_mean.exporting(mf_Ct);
	exp_l1_flux_mean.write_mesh();
	exp_l1_flux_mean.write_point_data(mf_Ct,l1_flux_mean , "l1_flux_mean");

//////////////////// Build l1_flux_residual
	vtk_export exp_l1_flux_residual(descr_transp.OUTPUT+"l1_flux_residual.vtk");
	exp_l1_flux_residual.exporting(mf_Ct);
	exp_l1_flux_residual.write_mesh();
	exp_l1_flux_residual.write_point_data(mf_Ct, l1_flux_residual, "l1_flux_residual");

//////////////////// Build l1_flux
	vtk_export exp_l1_flux(descr_transp.OUTPUT+"l1_flux"+dual+".vtk");
	exp_l1_flux.exporting(mf_Ct);
	exp_l1_flux.write_mesh();
	exp_l1_flux.write_point_data(mf_Ct,l1_flux , "l1_flux");
//////////////////////////////////////////////////////// l1_g

//////////////////// Build l1_g_ref
	vtk_export exp_l1_g_ref(descr_transp.OUTPUT+"l1_g_ref.vtk");
	exp_l1_g_ref.exporting(mf_Ct);
	exp_l1_g_ref.write_mesh();
	exp_l1_g_ref.write_point_data(mf_Ct, l1_g_ref, "l1_g_ref");

//////////////////// Build l1_g_mean
	vtk_export exp_l1_g_mean(descr_transp.OUTPUT+"l1_g_mean.vtk");
	exp_l1_g_mean.exporting(mf_Ct);
	exp_l1_g_mean.write_mesh();
	exp_l1_g_mean.write_point_data(mf_Ct, l1_g_mean, "l1_g_mean");

//////////////////// Build l1_g_residual
	vtk_export exp_l1_g_residual(descr_transp.OUTPUT+"l1_g_residual.vtk");
	exp_l1_g_residual.exporting(mf_Ct);
	exp_l1_g_residual.write_mesh();
	exp_l1_g_residual.write_point_data(mf_Ct, l1_g_residual, "l1_g_residual");

//////////////////// Build l1_g
	vtk_export exp_l1_g(descr_transp.OUTPUT+"l1_g"+dual+".vtk");
	exp_l1_g.exporting(mf_Ct);
	exp_l1_g.write_mesh();
	exp_l1_g.write_point_data(mf_Ct,l1_g , "l1_g");
//////////////////////////////////////////////////////// l1

//////////////////// Build l1_residual
	vtk_export exp_l1_residual(descr_transp.OUTPUT+"l1_residual.vtk");
	exp_l1_residual.exporting(mf_Ct);
	exp_l1_residual.write_mesh();
	exp_l1_residual.write_point_data(mf_Ct,l1_residual , "l1_residual");

//////////////////// Build l1
	vtk_export exp_l1(descr_transp.OUTPUT+"l1"+dual+".vtk");
	exp_l1.exporting(mf_Ct);
	exp_l1.write_mesh();
	exp_l1.write_point_data(mf_Ct,l1 , "l1");
//////////////////////////////////////////////////////// d1_grad

//////////////////// Build d1_grad_ref
	vtk_export exp_d1_grad_ref(descr_transp.OUTPUT+"d1_grad_ref.vtk");
	exp_d1_grad_ref.exporting(mf_Ct);
	exp_d1_grad_ref.write_mesh();
	exp_d1_grad_ref.write_point_data(mf_Ct, d1_grad_ref, "d1_grad_ref");

//////////////////// Build d1_grad_mean
	vtk_export exp_d1_grad_mean(descr_transp.OUTPUT+"d1_grad_mean.vtk");
	exp_d1_grad_mean.exporting(mf_Ct);
	exp_d1_grad_mean.write_mesh();
	exp_d1_grad_mean.write_point_data(mf_Ct, d1_grad_mean, "d1_grad_mean");

//////////////////// Build d1_grad_residual
	vtk_export exp_d1_grad_residual(descr_transp.OUTPUT+"d1_grad_residual.vtk");
	exp_d1_grad_residual.exporting(mf_Ct);
	exp_d1_grad_residual.write_mesh();
	exp_d1_grad_residual.write_point_data(mf_Ct,d1_grad_residual , "d1_grad_residual");

//////////////////// Build d1_grad
	vtk_export exp_d1_grad(descr_transp.OUTPUT+"d1_grad"+dual+".vtk");
	exp_d1_grad.exporting(mf_Ct);
	exp_d1_grad.write_mesh();
	exp_d1_grad.write_point_data(mf_Ct, d1_grad, "d1_grad");
//////////////////////////////////////////////////////// d1_flux

//////////////////// Build d1_flux_ref
	vtk_export exp_d1_flux_ref(descr_transp.OUTPUT+"d1_flux_ref.vtk");
	exp_d1_flux_ref.exporting(mf_Ct);
	exp_d1_flux_ref.write_mesh();
	exp_d1_flux_ref.write_point_data(mf_Ct, d1_flux_ref , "d1_flux_ref");

//////////////////// Build d1_flux_mean
	vtk_export exp_d1_flux_mean(descr_transp.OUTPUT+"d1_flux_mean.vtk");
	exp_d1_flux_mean.exporting(mf_Ct);
	exp_d1_flux_mean.write_mesh();
	exp_d1_flux_mean.write_point_data(mf_Ct, d1_flux_mean, "d1_flux_mean");

//////////////////// Build d1_flux_residual
	vtk_export exp_d1_flux_residual(descr_transp.OUTPUT+"d1_flux_residual.vtk");
	exp_d1_flux_residual.exporting(mf_Ct);
	exp_d1_flux_residual.write_mesh();
	exp_d1_flux_residual.write_point_data(mf_Ct,d1_flux_residual , "d1_flux_residual");

//////////////////// Build d1_flux
	vtk_export exp_d1_flux(descr_transp.OUTPUT+"d1_flux"+dual+".vtk");
	exp_d1_flux.exporting(mf_Ct);
	exp_d1_flux.write_mesh();
	exp_d1_flux.write_point_data(mf_Ct, d1_flux, "d1_flux");
//////////////////////////////////////////////////////// d1

//////////////////// Build d1_residual
	vtk_export exp_d1_residual(descr_transp.OUTPUT+"d1_residual.vtk");
	exp_d1_residual.exporting(mf_Ct);
	exp_d1_residual.write_mesh();
	exp_d1_residual.write_point_data(mf_Ct, d1_residual, "d1_residual");

//////////////////// Build d1
	vtk_export exp_d1(descr_transp.OUTPUT+"d1"+dual+".vtk");
	exp_d1.exporting(mf_Ct);
	exp_d1.write_mesh();
	exp_d1.write_point_data(mf_Ct, d1, "d1");
//////////////////////////////////////////////////////// eta 1

//////////////////// Build eta1_residual
	vtk_export exp_eta1_residual(descr_transp.OUTPUT+"eta1_residual.vtk");
	exp_eta1_residual.exporting(mf_Ct);
	exp_eta1_residual.write_mesh();
	exp_eta1_residual.write_point_data(mf_Ct, eta1_residual, "eta1_residual");

//////////////////// Build eta1
	vtk_export exp_eta1(descr_transp.OUTPUT+"eta1"+dual+".vtk");
	exp_eta1.exporting(mf_Ct);
	exp_eta1.write_mesh();
	exp_eta1.write_point_data(mf_Ct, eta1, "eta1");

/////////////////////////////////////////////////////// weight z1
	
	vtk_export exp_weight1(descr_transp.OUTPUT+"weight1"+dual+".vtk");
	exp_weight1.exporting(mf_Ct);
	exp_weight1.write_mesh();
	exp_weight1.write_point_data(mf_Ct, Z, "weight1");

/////////////////////////////////////////////////////// Primal solutions
	
//////////////////// Build u

	vtk_export exp_u(descr_transp.OUTPUT+"u.vtk");
	exp_u.exporting(mf_Ct);
	exp_u.write_mesh();
	exp_u.write_point_data(mf_Ct, u, "u");

//////////////////// Build U

	vtk_export exp_U(descr_transp.OUTPUT+"U.vtk");
	exp_U.exporting(mf_Cv);
	exp_U.write_mesh();
	exp_U.write_point_data(mf_Cv,U , "U");

//////////////////// Build U_Sigma

	vtk_export exp_U_Sigma(descr_transp.OUTPUT+"U_Sigma.vtk");
	exp_U_Sigma.exporting(mf_Ct_Sigma);
	exp_U_Sigma.write_mesh();
	exp_U_Sigma.write_point_data(mf_Ct_Sigma,U_Sigma , "U_Sigma");

//////////////////// Build U_Omega

	vtk_export exp_U_Omega(descr_transp.OUTPUT+"U_Omega.vtk");
	exp_U_Omega.exporting(mf_Ct);
	exp_U_Omega.write_mesh();
	exp_U_Omega.write_point_data(mf_Ct, U_Omega, "U_Omega");	
	};

	//! Compute model error of A2
	void transport3d1d::compute_model_error_2(const size_type DUALMODEL){

	#ifdef M3D1D_VERBOSE_
	cout << "  Computing model error for assumption A2 ..." << endl;
	#endif

//////////////////////////////////////////////////////// Estimators
	// Estimators for d1
	vector_type d2_residual(dof_transp.Ct()); gmm::clear(d2_residual); //d2_residual = (gradu, gradPhi)Sigma
	vector_type d2(dof_transp.Ct()); gmm::clear(d2);//d2 = (gradu, gradZ)Sigma

	vector_type l2_residual(dof_transp.Ct()); gmm::clear(l2_residual); //l2_residual = (f, Phi)Sigma
	vector_type l2(dof_transp.Ct()); gmm::clear(l2);//l2 = (f, Z)Sigma

	// Estimators for eta1
	vector_type eta2_residual(dof_transp.Ct()); gmm::clear(eta2_residual); //eta2_residual = l2_residual - d2_residual 
	vector_type eta2(dof_transp.Ct()); gmm::clear(eta2); //eta2 = l2 - d2


//////////////////////////////////////////////////////// weight Z (dual solution) 

	vector_type z(dof_transp.Ct()); gmm::clear(z); // this will store the weights Z
	// reference or reduced solution?
	vector_type Z(dof_transp.Ct()); gmm::clear(Z); // this will store the weights Z
	// reference or reduced solution?
	string dual="";

	if(DUALMODEL == REDUCED){
		gmm::copy(gmm::sub_vector(Z3, 
			gmm::sub_interval(0, dof_transp.Ct())), z);
		dual = "red";
	} 
	else if(DUALMODEL == REFERENCE ){		vector_type Z_Oplus(dof_transp.Ct_Omega()); gmm::clear(Z_Oplus);
		gmm::copy(gmm::sub_vector(Z1, 
			gmm::sub_interval(0, dof_transp.Ct_Omega())), Z_Oplus);
		getfem::interpolation(mf_Ct_Omega, mf_Ct, Z_Oplus, z, 2);
	} 	


//////////////////////////////////////////////////////// Primal solutions, u and U

	vector_type u(dof_transp.Ct()); gmm::clear(u); // u = U3_t (def su Omega)

	gmm::copy(gmm::sub_vector(U3, 
				  gmm::sub_interval(0, dof_transp.Ct())), u);


	

//////////////////////////////////////////////////////// parameters
	scalar_type k = PARAM.real_value("kk", "value of permeability for reduced problem");
	// Permeability
	vector_type Kv(dof.coefv());
	vector_type Kt(dof.coeft());
	constant_value=k;
	interpolation_function(mf_coefv, Kv, constant_function ); 
	interpolation_function(mf_coeft, Kt, constant_function ); 
	//Perimeter = 2*pi*R*k 
	vector_type Perimeter(dof.coefv());
	gmm::copy(param.R(), Perimeter);
	gmm::scale(Perimeter, 2*pi*k);
	//Area = pi*R^2
	vector_type Area(dof.coefv());
	gmm::copy(param.R(), Area);
	gmm::vscale(param.R(), Area);
	gmm::scale(Area, pi); 

//////////////////////////////////////////////////////// d2

//////////////////// Build d2_residual = (gradu, gradPhi)Sigma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  d2_residual..." << endl;
	#endif	
	sparse_matrix_type At(dof_transp.Ct(), dof_transp.Ct());gmm::clear(At);
	getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At, mimt, mf_Ct, descr_transp.SIGMA);
	gmm::mult(At, u, d2_residual);	

//////////////////// Build d2 = (gradu, gradZ)Sigma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  d2..." << endl;
	#endif	

	gmm::copy(d2_residual, d2);
	gmm::vscale(z, d2);

//////////////////////////////////////////////////////// l2

//////////////////// Build l2_residual = (f, Phi)Sigma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l2_residual..." << endl;
	#endif	
	vector_type F(dof.coeft()); gmm::clear(F);
	expr = PARAM.string_value("F_EXPR", "Expression for source term in tissue (muparser)");
	interpolation_function(mf_coeft, F, f_function ); 

	sparse_matrix_type M (dof_transp.Ct(),dof_transp.Ct()); gmm::clear(M);
	getfem::asm_mass_matrix(M, mimt, mf_Ct, mf_Ct, descr_transp.SIGMA);
	gmm::mult(M, z, l2_residual);	

//////////////////// Build l2 = (f, Z)Sigma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l2..." << endl;
	#endif	

	gmm::copy(l2_residual, l2);
	gmm::vscale(z, l2);

//////////////////////////////////////////////////////// eta2

//////////////////// Build eta2_residual = l2_residual - d2_residual 
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  eta2_residual..." << endl;
	#endif	

	gmm::add(l2_residual, gmm::scaled(d2_residual, -1.0), eta2_residual );
	gmm::add(eta2_residual, ETA_RESIDUAL);
//////////////////// Build eta2 = l2 - d2
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  eta2..." << endl;
	#endif	

	gmm::add(l2, gmm::scaled(d2, -1.0), eta2 );
	gmm::add(eta2, ETA);

/////////////////////////////////////////////////////// Global integrals
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  global integrals of d1 and l1..." << endl;
	#endif	

	cout <<"l2 = "<<gmm::vect_sp(z, l2_residual)<<endl;

	cout <<"d2 = "<<gmm::vect_sp(z, d2_residual)<<endl;

	cout <<"eta2 = "<<gmm::vect_sp(z, eta2_residual)<<endl;


	
////////////////////////////////////// E X P O R T //////////////////////////////////////
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting l2 and d2  ..." << endl;
	#endif
//////////////////////////////////////////////////////// d2

//////////////////// Build d2_residuals = (gradu, gradPhi)Sigma

	vtk_export exp_d2_residuals(descr_transp.OUTPUT+"d2_residuals.vtk");
	exp_d2_residuals.exporting(mf_Ct);
	exp_d2_residuals.write_mesh();
	exp_d2_residuals.write_point_data(mf_Ct, d2_residual , "d2_residuals");

//////////////////// Build d2 = (gradu, gradZ)Sigma

	vtk_export exp_d2(descr_transp.OUTPUT+"d2"+dual+".vtk");
	exp_d2.exporting(mf_Ct);
	exp_d2.write_mesh();
	exp_d2.write_point_data(mf_Ct,  d2, "d2");
//////////////////////////////////////////////////////// l2

//////////////////// Build l2_residuals = (f, Phi)Sigma
	vtk_export exp_l2_residuals(descr_transp.OUTPUT+"l2_residuals.vtk");
	exp_l2_residuals.exporting(mf_Ct);
	exp_l2_residuals.write_mesh();
	exp_l2_residuals.write_point_data(mf_Ct,l2_residual  , "l2_residuals");

//////////////////// Build l2 = (f, Z)Sigma
	vtk_export exp_l2(descr_transp.OUTPUT+"l2"+dual+".vtk");
	exp_l2.exporting(mf_Ct);
	exp_l2.write_mesh();
	exp_l2.write_point_data(mf_Ct, l2 , "l2");

//////////////////////////////////////////////////////// eta2

//////////////////// Build eta2_residual = l2_residual - d2_residual 
	vtk_export exp_eta2_residual(descr_transp.OUTPUT+"eta2_residual.vtk");
	exp_eta2_residual.exporting(mf_Ct);
	exp_eta2_residual.write_mesh();
	exp_eta2_residual.write_point_data(mf_Ct,  eta2_residual, "eta2_residual");
//////////////////// Build eta2 = l2 - d2
	vtk_export exp_eta2(descr_transp.OUTPUT+"eta2"+dual+".vtk");
	exp_eta2.exporting(mf_Ct);
	exp_eta2.write_mesh();
	exp_eta2.write_point_data(mf_Ct, eta2 , "eta2");

/////////////////////////////////////////////////////// weight z2
	
	vtk_export exp_weight2(descr_transp.OUTPUT+"weight2"+dual+".vtk");
	exp_weight2.exporting(mf_Ct);
	exp_weight2.write_mesh();
	exp_weight2.write_point_data(mf_Ct, z, "weight2");

	};


	//! Compute model error of A3	
	void transport3d1d::compute_model_error_3(const size_type DUALMODEL){
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing model error for assumption A3 ..." << endl;
	#endif


//////////////////////////////////////////////////////// Estimators
	// Estimators for d3

	vector_type d3_ref(dof_transp.Ct()); gmm::clear(d3_ref); //d3_ref = (k u, Phi)Gamma
	vector_type d3_mean(dof_transp.Ct()); gmm::clear(d3_mean); //d3_mean = (k baru, Phi)Gamma	
	vector_type d3_residual(dof_transp.Ct()); gmm::clear(d3_residual); //d3_residual = (k (i-bar)u, Phi)Gamma
	vector_type d3(dof_transp.Ct()); gmm::clear(d3); //d3 = (k (i-bar)u, z)Gamma

	// Estimators for l3
	vector_type l3_ref(dof_transp.Ct()); gmm::clear(l3_ref); //l3_ref = (k U, Phi)Gamma
	vector_type l3_mean(dof_transp.Ct()); gmm::clear(l3_mean); //l3_mean = (k barU, Phi)Gamma	
	vector_type l3_residual(dof_transp.Ct()); gmm::clear(l3_residual); //l3_residual = (k (i-bar)U, Phi)Gamma
	vector_type l3(dof_transp.Ct()); gmm::clear(l3); //l3 = (k (i-bar)U, z)Gamma

	// Estimators for eta3
	vector_type eta3_residual(dof_transp.Ct()); gmm::clear(eta3_residual); //eta3_residual = l3_residual - d3_residual 
	vector_type eta3(dof_transp.Ct()); gmm::clear(eta3); //eta3 = l3 - d3	


//////////////////////////////////////////////////////// weight Z (dual solution) 

	vector_type z(dof_transp.Ct()); gmm::clear(z); // this will store the weights Z
	// reference or reduced solution?
	string dual="";

	if(DUALMODEL == REDUCED){
		gmm::copy(gmm::sub_vector(Z3, 
			gmm::sub_interval(0,dof_transp.Ct())), z);
		dual = "red";
	} 
	else if(DUALMODEL == REFERENCE ){	
		gmm::copy(gmm::sub_vector(Z2, 
			gmm::sub_interval(0,dof_transp.Ct())), z);
	} 	


//////////////////////////////////////////////////////// Primal solutions, u and U

	vector_type u(dof_transp.Ct()); gmm::clear(u); // u = U3_t (def su Omega)
	vector_type U(dof_transp.Cv()); gmm::clear(U); // U = U3_v (def su lambda)	
	vector_type U_Sigma(dof_transp.Ct_Sigma()); gmm::clear(U_Sigma); // U_Sigma = U (interp su Sigma, extra = 2)
	vector_type U_Omega(dof_transp.Ct()); gmm::clear(U_Omega); // U_Omega = U (interp su Omega, extra = 0)

	gmm::copy(gmm::sub_vector(U3, 
				  gmm::sub_interval(0, dof_transp.Ct())), u);
	gmm::copy(gmm::sub_vector(U3, 
				  gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv())), U);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, U, U_Sigma, 2);
	getfem::interpolation(mf_Ct_Sigma, mf_Ct, U_Sigma, U_Omega, 0);


	

//////////////////////////////////////////////////////// parameters
	scalar_type k = PARAM.real_value("kk", "value of permeability for reduced problem");
	// Permeability
	vector_type Kv(dof.coefv());
	vector_type Kt(dof.coeft());
	constant_value=k;
	interpolation_function(mf_coefv, Kv, constant_function ); 
	interpolation_function(mf_coeft, Kt, constant_function ); 
	//Perimeter = 2*pi*R*k 
	vector_type Perimeter(dof.coefv());
	gmm::copy(param.R(), Perimeter);
	gmm::scale(Perimeter, 2*pi*k);
	//Area = pi*R^2
	vector_type Area(dof.coefv());
	gmm::copy(param.R(), Area);
	gmm::vscale(param.R(), Area);
	gmm::scale(Area, pi); 


//////////////////////////////////////////////////////// d3

//////////////////// Build d3_ref = (k u, Phi)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  d3_ref..." << endl;
	#endif
	sparse_matrix_type M_Omega (dof_transp.Ct(),dof_transp.Ct()); gmm::clear(M_Omega);
	getfem::asm_mass_matrix_param(M_Omega, mimt, mf_Ct, mf_Ct, mf_coeft, Kt, descr_transp.GAMMA);
	gmm::mult(M_Omega, u, d3_ref);	

//////////////////// Build d3_mean = (k baru, Phi)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing d3_mean ..." << endl;
	#endif
	vector_type u_mean(dof_transp.Cv()); gmm::clear(u_mean);
	vector_type u_sigma(dof_transp.Ct_Sigma()); gmm::clear(u_sigma);
	vector_type u_bar(dof_transp.Ct()); gmm::clear(u_bar);
	gmm::mult(MBAR, u, u_mean);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, u_mean, u_sigma, 2);
	getfem::interpolation(mf_Ct_Sigma, mf_Ct, u_sigma, u_bar, 0);
	gmm::mult(M_Omega, u_bar, d3_mean);

//////////////////// Build d3_residual = (k (i-bar)u, Phi)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing d3_residual ..." << endl;
	#endif
	gmm::add(d3_ref, gmm::scaled(d3_mean, -1.0), d3_residual);

//////////////////// Build d3_flux = (k (i-bar)u, Z)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing d3_flux ..." << endl;
	#endif
	gmm::copy(d3_residual, d3);
	gmm::vscale(z, d3);

//////////////////////////////////////////////////////// l3

//////////////////// Build l3_ref = (k u, Phi)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l3_ref..." << endl;
	#endif
	gmm::mult(M_Omega, U_Omega, l3_ref);	

//////////////////// Build l3_mean = (k barU, Phi)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l3_flux_mean..." << endl;
	#endif
	vector_type U_mean(dof_transp.Cv()); gmm::clear(U_mean);
	vector_type U_sigma(dof_transp.Ct_Sigma()); gmm::clear(U_sigma);
	vector_type U_bar(dof_transp.Ct()); gmm::clear(U_bar);
	gmm::mult(MBAR, U_Omega, U_mean);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, U_mean, U_sigma, 2);
	getfem::interpolation(mf_Ct_Sigma, mf_Ct, U_sigma, U_bar, 0);
	gmm::mult(M_Omega, U_bar, l3_mean);

//////////////////// Build l3_residual = (k (i-bar)U, Phi)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l3_flux_residual..." << endl;
	#endif
	gmm::add(l3_ref, gmm::scaled(l3_mean, -1.0), l3_residual);

//////////////////// Build l3_flux = (k (i-bar)u, Z)Gamma
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  l3_flux..." << endl;
	#endif
	gmm::copy(l3_residual, l3);
	gmm::vscale(z, l3);

//////////////////////////////////////////////////////// eta 3

//////////////////// Build eta3_residual = l3_residual - d3_residual
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing eta3_residual ..." << endl;
	#endif
	gmm::add(l3_residual, gmm::scaled(d3_residual, -1.0), eta3_residual );
	gmm::add(gmm::scaled(d3_residual, -1.0), ETA_RESIDUAL);
//////////////////// Build eta3 = l3 - d3
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing eta3 ..." << endl;
	#endif
	gmm::add(l3, gmm::scaled(d3, -1.0), eta3 );
	gmm::add(gmm::scaled(d3, -1.0), ETA);
/////////////////////////////////////////////////////// Global integrals
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing  global integrals of d3 and l3..." << endl;
	#endif

	cout <<"d3_ref = "<<gmm::vect_sp(z, d3_ref)<<endl;
	cout <<"d3_mean = "<<gmm::vect_sp(z, d3_mean)<<endl;
	cout <<"d3 = "<<gmm::vect_sp(z, d3_residual)<<endl;

	cout <<"l3_ref = "<<gmm::vect_sp(z, l3_ref)<<endl;
	cout <<"l3_mean = "<<gmm::vect_sp(z, l3_mean)<<endl;
	cout <<"l3 = "<<gmm::vect_sp(z, l3_residual)<<endl;

	cout <<"eta3 = "<<gmm::vect_sp(z, eta3_residual)<<endl;



////////////////////////////////////// E X P O R T //////////////////////////////////////
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting l3 and d3  ..." << endl;
	#endif
//////////////////////////////////////////////////////// l3_flux

//////////////////// Build l3_ref
	vtk_export exp_l3_ref(descr_transp.OUTPUT+"l3_ref.vtk");
	exp_l3_ref.exporting(mf_Ct);
	exp_l3_ref.write_mesh();
	exp_l3_ref.write_point_data(mf_Ct, l3_ref, "l3_ref");

//////////////////// Build l3_mean
	vtk_export exp_l3_mean(descr_transp.OUTPUT+"l3_mean.vtk");
	exp_l3_mean.exporting(mf_Ct);
	exp_l3_mean.write_mesh();
	exp_l3_mean.write_point_data(mf_Ct,l3_mean , "l3_mean");

//////////////////// Build l3_residual
	vtk_export exp_l3_residual(descr_transp.OUTPUT+"l3_residual.vtk");
	exp_l3_residual.exporting(mf_Ct);
	exp_l3_residual.write_mesh();
	exp_l3_residual.write_point_data(mf_Ct, l3_residual, "l3_residual");

//////////////////// Build l3
	vtk_export exp_l3(descr_transp.OUTPUT+"l3"+dual+".vtk");
	exp_l3.exporting(mf_Ct);
	exp_l3.write_mesh();
	exp_l3.write_point_data(mf_Ct,l3 , "l3");

//////////////////////////////////////////////////////// d3_grad

//////////////////// Build d3_ref
	vtk_export exp_d3_ref(descr_transp.OUTPUT+"d3_ref.vtk");
	exp_d3_ref.exporting(mf_Ct);
	exp_d3_ref.write_mesh();
	exp_d3_ref.write_point_data(mf_Ct, d3_ref, "d3_ref");

//////////////////// Build d3_mean
	vtk_export exp_d3_mean(descr_transp.OUTPUT+"d3_mean.vtk");
	exp_d3_mean.exporting(mf_Ct);
	exp_d3_mean.write_mesh();
	exp_d3_mean.write_point_data(mf_Ct, d3_mean, "d3_mean");

//////////////////// Build d3_residual
	vtk_export exp_d3_residual(descr_transp.OUTPUT+"d3_residual.vtk");
	exp_d3_residual.exporting(mf_Ct);
	exp_d3_residual.write_mesh();
	exp_d3_residual.write_point_data(mf_Ct, d3_residual, "d3_residual");

//////////////////// Build d3
	vtk_export exp_d3(descr_transp.OUTPUT+"d3"+dual+".vtk");
	exp_d3.exporting(mf_Ct);
	exp_d3.write_mesh();
	exp_d3.write_point_data(mf_Ct, d3, "d3");
//////////////////////////////////////////////////////// eta 3

//////////////////// Build eta3_residual
	vtk_export exp_eta3_residual(descr_transp.OUTPUT+"eta3_residual.vtk");
	exp_eta3_residual.exporting(mf_Ct);
	exp_eta3_residual.write_mesh();
	exp_eta3_residual.write_point_data(mf_Ct, eta3_residual, "eta3_residual");

//////////////////// Build eta3
	vtk_export exp_eta3(descr_transp.OUTPUT+"eta3"+dual+".vtk");
	exp_eta3.exporting(mf_Ct);
	exp_eta3.write_mesh();
	exp_eta3.write_point_data(mf_Ct, eta3, "eta3");

/////////////////////////////////////////////////////// weight z3
	
	vtk_export exp_weight3(descr_transp.OUTPUT+"weight3"+dual+".vtk");
	exp_weight3.exporting(mf_Ct);
	exp_weight3.write_mesh();
	exp_weight3.write_point_data(mf_Ct, z, "weight3");

	};



	/////////////////////////////////////////////
	//User interface

	//User function to compute all the the steps (reduced primal and Z0,Z1,Z2)
	void transport3d1d::model_error(int argc, char *argv[]){

	double timeA0=0; 
	double timeA1=0;
	double timeA2=0;
	double timeA3p=0;
	double timeA3d=0;
	double time_model_ref=0;
	double time_model_red=0;


 	double time = gmm::uclock_sec();

		//initialize 
		init_model(argc, argv);
	string DUALMODEL = PARAM.string_value("DUALMODEL", "Choose which dual solution use (reference or dual or both)");




	if(DUALMODEL == "REFERENCE"|| DUALMODEL == "BOTH"){
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"Z0: solve reference problem:"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	
		//assemble           
		assembly_model(A0, DUAL);    
		//solve     
		if (!solve_model(A0,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A0,DUAL);

	timeA0= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();
	
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"Z1: solve problem for first assumption:"<<endl;
	cout <<"    U(s,r,theta)=U(s)"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	 
		//assemble        
		assembly_model(A1, DUAL);    
		//solve     
		if (!solve_model(A1,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A1,DUAL);

	timeA1= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();

	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"Z2: solve problem for second assumption:"<<endl;
	cout <<"    Omega+ == Omega"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

		//assemble         
		assembly_model(A2, DUAL);    
		//solve     
		if (!solve_model(A2,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A2,DUAL);


	timeA2= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();
	}

	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"Z3: solve problem for third assumption:"<<endl;
	cout <<"    Neglect fluctuations"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	
	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble      
		assembly_model(A3);    
		//solve     
		if (!solve_model(A3,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A3,PRIMAL);

	timeA3p= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();

	if(DUALMODEL == "REDUCED" || DUALMODEL == "BOTH"){
	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
  

		//assemble        
		assembly_model(A3, DUAL);    
		//solve     
		if (!solve_model(A3,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A3,DUAL);

	timeA3d= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();
	}
	//Model error
	if(DUALMODEL == "REFERENCE" || DUALMODEL == "BOTH"){
	gmm::clear(ETA_RESIDUAL); gmm::clear(ETA);
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D1 AND L1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_1(REFERENCE);
	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D2 AND L2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_2(REFERENCE);


	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D3 AND L3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_3(REFERENCE);

	// rho = rho1 +rho2 + rho3 and eta = eta1 + eta2 + eta3
	vtk_export exp_ETA_RESIDUAL(descr_transp.OUTPUT+"ETA_RESIDUAL.vtk");
	exp_ETA_RESIDUAL.exporting(mf_Ct);
	exp_ETA_RESIDUAL.write_mesh();
	exp_ETA_RESIDUAL.write_point_data(mf_Ct, ETA_RESIDUAL, "ETA_RESIDUAL");

	vtk_export exp_ETA(descr_transp.OUTPUT+"ETA.vtk");
	exp_ETA.exporting(mf_Ct);
	exp_ETA.write_mesh();
	exp_ETA.write_point_data(mf_Ct, ETA, "ETA");

	time_model_ref= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();
	}

	if(DUALMODEL == "REDUCED" || DUALMODEL == "BOTH"){
	gmm::clear(ETA_RESIDUAL); gmm::clear(ETA);
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D1 AND L1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_1(REDUCED);
	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D2 AND L2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_2(REDUCED);


	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D3 AND L3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_3(REDUCED);
	vtk_export exp_ETA_RESIDUAL(descr_transp.OUTPUT+"ETA_RESIDUAL_RED.vtk");
	exp_ETA_RESIDUAL.exporting(mf_Ct);
	exp_ETA_RESIDUAL.write_mesh();
	exp_ETA_RESIDUAL.write_point_data(mf_Ct, ETA_RESIDUAL, "ETA_RESIDUAL");

	vtk_export exp_ETA(descr_transp.OUTPUT+"ETA_RED.vtk");
	exp_ETA.exporting(mf_Ct);
	exp_ETA.write_mesh();
	exp_ETA.write_point_data(mf_Ct, ETA, "ETA");

	time_model_red= -time+gmm::uclock_sec();
	time= gmm::uclock_sec();
	}




	cout << "========================================================="<<endl;
	cout << endl<<"... time to solve Z0: " << timeA0 << " seconds\n";
	cout << endl<<"... time to solve Z1: " << timeA1 << " seconds\n";
	cout << endl<<"... time to solve Z2: " << timeA2 << " seconds\n";
	cout << endl<<"... time to solve Z3: " << timeA3d << " seconds\n";
	cout << endl<<"... time to solve U3: " << timeA3p << " seconds\n";
	if(DUALMODEL == "REFERENCE" || DUALMODEL == "BOTH"){
	cout << endl<<"... time to solve error reference estimators: " << time_model_ref << " seconds\n";}
	if(DUALMODEL == "REDUCED" || DUALMODEL == "BOTH"){
	cout << endl<<"... time to solve error reduced estimators: " << time_model_red << " seconds\n";}

};




	//User method to compute model error for A1
	void transport3d1d::model_error_A0(int argc, char *argv[]){// A1: U=U(s)
	//Primal	
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A0: solve reference problem:"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

		//initialize 
		init_model(argc, argv);


	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U0 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A0);    
		//solve     
		if (!solve_model(A0,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A0,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z0 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A0, DUAL);   
		//solve     
		if (!solve_model(A0,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A0,DUAL);




};

	//User method to compute model error for A1
	void transport3d1d::model_error_A1(int argc, char *argv[]){// A1: U=U(s)
	//Primal	
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A1: solve problem for first assumption:"<<endl;
	cout <<"    U(s,r,theta)=U(s)"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;	  

		//initialize 
		init_model(argc, argv);


	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A1);    
		//solve     
		if (!solve_model(A1,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A1,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A1, DUAL);     
		//solve     
		if (!solve_model(A1,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A1,DUAL);

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D1 AND L1 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_1(REFERENCE);


};


	//User method to compute model error for A2
	void transport3d1d::model_error_A2(int argc, char *argv[]){//A2: Omega+ == Omega
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A2: solve problem for second assumption:"<<endl;
	cout <<"    Omega+ == Omega"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;

		//initialize 
		init_model(argc, argv);

	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A2);    
		//solve     
		if (!solve_model(A2,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A2,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A2, DUAL);      
		//solve     
		if (!solve_model(A2,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A2,DUAL);

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D2 AND L2 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_2(REFERENCE);
		

	};


	//User method to compute model error for A3
	void transport3d1d::model_error_A3(int argc, char *argv[]){//A3: neglect fluctuations

	cout <<endl<< "================================================"<<endl<<endl;
	cout <<"A3: solve problem for third assumption:"<<endl;
	cout <<"    Neglect fluctuations"<<endl;
	cout <<endl<< "================================================"<<endl<<endl;

		//initialize 
		init_model(argc, argv);
	//Primal
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" PRIMAL: U3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A3);    
		//solve     
		if (!solve_model(A3,PRIMAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A3,PRIMAL);

	//Dual
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" DUAL: Z3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//assemble        
		assembly_model(A3, DUAL);      
		//solve     
		if (!solve_model(A3,DUAL)) GMM_ASSERT1(false, "solve procedure has failed");  
		//export
   		export_model(A3,DUAL);

	//Model error
	cout <<endl<< "================================================"<<endl<<endl;
	cout <<" MODEL ERROR: D3 AND L3 "<<endl;
	cout <<endl<< "================================================"<<endl<<endl;
		//compute error
		compute_model_error_3(REFERENCE);

	};

/* 

	/////////////////////////////////////////////
	// Compute model errors -- testing! compute extensively many terms in order to try and validate the code

	//! Compute model error of A1	
	void transport3d1d::compute_model_error_1_test(){
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing model error for assumption A1 ..." << endl;
	#endif

	// vari integrali per d1
	vector_type d1_1(dof_transp.Ct_Sigma()); gmm::clear(d1_1); //d1_1 = int su Sigma dei gradienti di Ub(su Sigma) e z1ref(su Sigma)
	vector_type d1_2(dof_transp.Ct_Sigma()); gmm::clear(d1_2); //d1_2 = int su Lambda delle medie barbar dei gradienti (costruisci prima la matrice)
	vector_type d1_3(dof_transp.Ct_Sigma()); gmm::clear(d1_3); //d1_3 = int su Sigma del flusso (matrice di massa) con Ub(su Sigma) e z1ref(su Sigma)
	vector_type d1_4(dof_transp.Ct_Sigma()); gmm::clear(d1_4); //d1_4 = int su Lambda delle medie bar dei flussi (costruisci prima la matrice)
	vector_type d1_5(dof_transp.Ct_Sigma()); gmm::clear(d1_5); //d1_5 = int su Sigma delle medie barbar dei gradienti  (costruisci prima le medie)
	vector_type d1_6(dof_transp.Ct_Sigma()); gmm::clear(d1_6); //d1_6 = int su Gamma delle medie bar dei flussi (costruisci prima le medie)
	vector_type d1_7(dof_transp.Cv()); gmm::clear(d1_7); //d1_7 = int su Lambda delle medie barbar dei gradienti (calcola le medie e poi fai l'integrale su lambda)
	
	vector_type d1_grad(dof_transp.Ct_Sigma()); gmm::clear(d1_grad); //d1_grad = d1_1 - d1_5
	vector_type d1_flux(dof_transp.Ct_Sigma()); gmm::clear(d1_flux); //d1_flux = d1_3 - d1_6

	// vari integrali per l1_t
	vector_type l1_t1(dof_transp.Ct()); gmm::clear(l1_t1); // l1_t1 = (u, z1ref)Gamma con asm_source_term (localizzazione su u)
	vector_type l1_t2(dof_transp.Cv()); gmm::clear(l1_t2); // l1_t2 = MBAR*l1_t1
	vector_type l1_t3(dof_transp.Ct()); gmm::clear(l1_t3); // l1_t3 = (u, z1ref)Gamma con asm_mass_term_param (localizzazione su u)
	vector_type l1_t4(dof_transp.Ct()); gmm::clear(l1_t4); // l1_t4 = (u_bar, z1ref)Lambda (costruisci prima la matrice) (localizzazione su u)
	vector_type l1_t5(dof_transp.Ct_Sigma()); gmm::clear(l1_t5); // l1_t5 = (u_bar, z1ref)Lambda (costruisci prima la matrice) (localizzazione su z)
	vector_type l1_t6(dof_transp.Ct_Sigma()); gmm::clear(l1_t6); // l1_t6 = l1_t1 (localizzazione su z)
	vector_type l1_t7(dof_transp.Ct_Sigma()); gmm::clear(l1_t7); // l1_t7 = l1_t3 (localizzazione su z)
	vector_type l1_t8(dof_transp.Ct_Sigma()); gmm::clear(l1_t8); // l1_t8 = (u_bar, z1ref)Gamma (costruisci prima le medie) (localizzazione su u)
	vector_type l1_t(dof_transp.Ct()); gmm::clear(l1_t); // l1_t = L1_t1 - l1_t4

	// vari integrali per l1_g
	vector_type l1_g1(dof_transp.Ct()); gmm::clear(l1_g1); // l1_g1 = (g,z1ref)_Sigma con asm_source term (localizzato su g)
	vector_type l1_g2(dof_transp.Ct()); gmm::clear(l1_g2); // l1_g2 = (g_bar,z1ref)_lambda (costruisci prima la matrice)
	vector_type l1_g3(dof_transp.Ct_Sigma()); gmm::clear(l1_g3); // l1_g3 = (g,z1ref)_Sigma con asm_source term (localizzato su z)
	vector_type l1_g4(dof_transp.Ct_Sigma()); gmm::clear(l1_g4); // l1_g4 = (g_bar,z1ref)_sigma (costruisci prima le medie)
	vector_type l1_g(dof_transp.Ct()); gmm::clear(l1_g); // l1_g = l1_g1-l1_g2

	// Solution of dual 1D problem
	vector_type z1ref(dof_transp.Ct_Sigma()); gmm::clear(z1ref); //z1ref = Z0_v (def su Sigma)
	vector_type z1refa(dof_transp.Ct()); gmm::clear(z1refa); //z1refa = z1ref(interp su mf_Ct, extra=0)
	vector_type z1refb(dof_transp.Ct()); gmm::clear(z1refb); //z1refb = z1ref(interp su mf_Ct, extra=1)
	vector_type z1refc(dof_transp.Ct()); gmm::clear(z1refc); //z1refc = z1ref(interp su mf_Ct, extra=2)
		//importo i valori delle z
	gmm::copy(gmm::sub_vector(Z0, 
		gmm::sub_interval(dof_transp.Ct_Omega(), dof_transp.Ct_Sigma())), z1ref);
	getfem::interpolation(mf_Ct_Sigma, mf_Ct, z1ref, z1refa, 0);
	getfem::interpolation(mf_Ct_Sigma, mf_Ct, z1ref, z1refb, 1);
	getfem::interpolation(mf_Ct_Sigma, mf_Ct, z1ref, z1refc, 2);

	// Solution of primal 3D problem	
	vector_type u(dof_transp.Ct()); gmm::clear(u); // u = U3_t (def su Omega)

	vector_type U(dof_transp.Cv()); gmm::clear(U); // U = U3_v (def su lambda)
	vector_type Ua(dof_transp.Ct_Sigma()); gmm::clear(Ua); // Ua = U (interp su Sigma, extra = 0)
	vector_type Ub(dof_transp.Ct_Sigma()); gmm::clear(Ub); // Ub = U (interp su Sigma, extra = 2)
	vector_type Uc(dof_transp.Ct()); gmm::clear(Uc); // Uc = U (interp su Omega, extra = 0)
	vector_type Ud(dof_transp.Ct_Sigma()); gmm::clear(Ud); // Ud = U (interp su Sigma, extra = 1)
		// importo i valori di u e U
	gmm::copy(gmm::sub_vector(U3, 
				  gmm::sub_interval(0, dof_transp.Ct())), u);
	gmm::copy(gmm::sub_vector(U3, 
				  gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv())), U);

	getfem::interpolation(mf_Cv, mf_Ct_Sigma, U, Ua, 0);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, U, Ub, 2);
	getfem::interpolation(mf_Cv, mf_Ct, U, Uc, 2);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, U, Ud, 1);

	// valori double dei vari integrali (globali)
	scalar_type int_l1_t1=0, int_l1_t2=0, int_l1_t=0; //int_l1_t1= l1_t3, //int_l1_t2= l1_t8,   
	scalar_type int_l1_g1=0, int_l1_g2=0, int_l1_g=0;
	scalar_type int_j1=0;
	scalar_type int_d1_1=0, int_d1_2=0, int_d1_grad=0;
	scalar_type int_d1_3=0, int_d1_4=0, int_d1_flux=0;
	scalar_type int_d1=0;
	

	//parameters
	scalar_type k = PARAM.real_value("kk", "value of permeability for reduced problem");
	// Permeability
	vector_type Kv(dof.coefv());
	vector_type Kt(dof.coeft());
	constant_value=k;
	interpolation_function(mf_coefv, Kv, constant_function ); 
	interpolation_function(mf_coeft, Kt, constant_function ); 
	//Perimeter = 2*pi*R*k 
	vector_type Perimeter(dof.coefv());
	gmm::copy(param.R(), Perimeter);
	gmm::scale(Perimeter, 2*pi*k);
	//Area = pi*R^2
	vector_type Area(dof.coefv());
	gmm::copy(param.R(), Area);
	gmm::vscale(param.R(), Area);
	gmm::scale(Area, pi); 

////////////////////////////////////// l1_t

	//l1_t1 = (u, z1ref)Gamma con asm_source_term (localizzazione su u)
	asm_source_term(l1_t1, mimt, mf_Ct, mf_Ct_Sigma, z1ref,descr_transp.GAMMA); 	//manca k!!!!
	gmm::vscale(u, l1_t1);

	// l1_t2 = MBAR*l1_t1
	gmm::mult(MBAR, l1_t1, l1_t2);

	// l1_t6 = l1_t1 (localizzazione su z)
	asm_source_term(l1_t6, mimt, mf_Ct_Sigma,mf_Ct, u,descr_transp.GAMMA); 		//manca k!!!!
	gmm::vscale(z1ref, l1_t6);

	// diverse matrici di massa / exchange utili per il proseguio		
	sparse_matrix_type M_OS (dof_transp.Ct(),dof_transp.Ct_Sigma()); gmm::clear(M_OS); // mass matrix tra Ct e Ct_Sigma, su Gamma, moltiplicata per K (vedi l1_t3)
	sparse_matrix_type M_OS2 (dof_transp.Ct(),dof_transp.Ct_Sigma()); gmm::clear(M_OS2); 
	sparse_matrix_type M_OS3 (dof_transp.Ct(),dof_transp.Ct_Sigma()); gmm::clear(M_OS3); // bar matrix tra Ct e Ct_Sigma (media su entrambe le basi) (vedi l1_t4)
	sparse_matrix_type M_OS4 (dof_transp.Ct(),dof_transp.Ct_Sigma()); gmm::clear(M_OS4); // barbar matrix tra Ct e Ct_Sigma (media su entrambe le basi) (vedi l1_g2)
	sparse_matrix_type M_OS5 (dof_transp.Ct(),dof_transp.Ct_Sigma()); gmm::clear(M_OS5); 
	sparse_matrix_type M_OS6 (dof_transp.Ct_Sigma(), dof_transp.Ct()); gmm::clear(M_OS6); // mass matrix tra Ct e Ct_Sigma, su Gamma, moltiplicata per K (vedi l1_t7)
	sparse_matrix_type M_OS8 (dof_transp.Ct_Sigma(), dof_transp.Ct_Sigma()); gmm::clear(M_OS8);// mass matrix su Ct_Sigma, su Sigma (vedi l1_g4)
	sparse_matrix_type M_OS9 (dof_transp.Ct_Sigma(), dof_transp.Ct_Sigma()); gmm::clear(M_OS9);

	// l1_t3 = (u, z1ref)Gamma con asm_mass_term_param (localizzazione su u)	
	getfem::asm_mass_matrix_param(M_OS, mimt, mf_Ct, mf_Ct_Sigma, mf_coeft, Kt, descr_transp.GAMMA);
	gmm::mult(M_OS, z1ref, l1_t3);	
int_l1_t1= gmm::vect_sp(u, l1_t3); cout <<"l1_t mass= "<<int_l1_t1<<endl;
	gmm::vscale(u, l1_t3);
	
	// l1_t3 = (u, z1ref)Gamma con asm_mass_term_param (localizzazione su z)	
	getfem::asm_mass_matrix(M_OS6, mimt, mf_Ct_Sigma, mf_Ct, descr_transp.GAMMA);
	gmm::mult(M_OS6, u , l1_t7);	
	gmm::vscale(z1ref, l1_t7);

	// l1_t4 = (u_bar, z1ref)Lambda (costruisci prima la matrice)
		// Mvv = matrice di massa su Lambda (moltiplicata per 2*pi*R*K
	sparse_matrix_type Mvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Mvv);
	getfem::asm_mass_matrix_param (Mvv, mimv, mf_Cv,mf_coefv, Perimeter);
		// INT_Omega_Sigma = matrice di interpolazione tra Omega e Sigma
	sparse_matrix_type INT_Omega_Sigma(dof_transp.Ct_Sigma(), dof_transp.Ct());gmm::clear(INT_Omega_Sigma);
 	getfem::interpolation(mf_Ct, mf_Ct_Sigma, INT_Omega_Sigma, 2);
		// MBAR = media tra Omega e Lambda, Mbar_Sigma = media tra Sigma e Lambda
	sparse_matrix_type Mbar_Sigma(dof_transp.Cv(), dof_transp.Ct_Sigma());gmm::clear(Mbar_Sigma);
	gmm::mult(MBAR, gmm::transposed(INT_Omega_Sigma) , Mbar_Sigma);
		//M_OS3 = MBAR^T * Mvv * Mbar_Sigma
	sparse_matrix_type M_partial(dof_transp.Cv(), dof_transp.Ct_Sigma());gmm::clear(M_partial);
	gmm::mult(Mvv, Mbar_Sigma, M_partial);
	gmm::mult(gmm::transposed(MBAR), M_partial, M_OS3);
	gmm::mult(M_OS3, z1ref, l1_t4);	
int_l1_t2= gmm::vect_sp(u, l1_t4); cout <<"l1_t mean = "<<int_l1_t2<<endl;
	gmm::vscale(u, l1_t4);

	// l1_t = L1_t1 - l1_t4
	gmm::add(l1_t1, gmm::scaled(l1_t4,-1.0), l1_t);

	// l1_t5 = (u_bar, z1ref)Lambda (costruisci prima la matrice) (localizzazione su z)
	gmm::mult(gmm::transposed(M_OS3), u, l1_t5);	
	gmm::vscale(z1ref, l1_t5);
	

	// l1_t8 = (u_bar, z1ref)Lambda (costruisci prima le medie) (localizzazione su u)
	vector_type u_mean(dof_transp.Cv()); gmm::clear(u_mean);
	vector_type u_sigma(dof_transp.Ct_Sigma()); gmm::clear(u_sigma);
	gmm::mult(MBAR, u, u_mean);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, u_mean, u_sigma, 2);
	getfem::asm_mass_matrix_param(M_OS9, mimt, mf_Ct_Sigma, mf_coeft, Kt, descr_transp.GAMMA);
	gmm::mult(M_OS9, u_sigma, l1_t8);
	cout <<"l1_t mean alternativo = "<<gmm::vect_sp(z1ref, l1_t8)<<endl;
	gmm::vscale(z1ref, l1_t8);


///////////////////////////////// li_g

		
		vector_type G(dof_transp.Cv());gmm::clear(G);
		vector_type Gt(dof_transp.Ct()); gmm::clear(Gt);
		expr = PARAM.string_value("G_EXPR", "Expression for source term in vessel (muparser)");
		bool G_CONSTANT = PARAM.int_value("G_CONSTANT", "flag for using constant source term g in vessel");
		// If G is constant in the cross section, i just interpolate the g_function in Lambda

	// matrice di scambio (medie su D) tra lambda e sigma
	sparse_matrix_type MBARBAR_Sigma(dof_transp.Cv(), dof_transp.Ct_Sigma());gmm::clear(MBARBAR_Sigma);
		if(G_CONSTANT){
			//do nothing! l1_g=0
		}
		else{	//If G is not constant, I must build Mbarbar
			interpolation_function(mf_Cv, G, g_function );		
			interpolation_function(mf_Ct, Gt, g_function ); 		
	// l1_g1 = (g,z1ref)_Sigma con asm_source term (localizzato su g)
	asm_source_term(l1_g1, mimt, mf_Ct, mf_Ct_Sigma, z1ref,descr_transp.SIGMA);
int_l1_g1= gmm::vect_sp(Gt, l1_g1); cout <<"l1_g mass = "<<int_l1_g1<<endl;
	gmm::vscale(Gt, l1_g1);	

	// l1_g2 = (g_bar,z1ref)_lambda (costruisci prima la matrice)
		// Mvv = matrice di massa su Lambda (moltiplicata per pi*R^2	
	sparse_matrix_type Mvvg(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Mvvg);
	getfem::asm_mass_matrix_param (Mvvg, mimv, mf_Cv,mf_coefv, Area);
		// MBARBAR = media tra Omega e Lambda, Mbarbar_Sigma = media tra Sigma e Lambda
	sparse_matrix_type Mbarbar_Sigma(dof_transp.Cv(), dof_transp.Ct_Sigma());gmm::clear(Mbarbar_Sigma);
	gmm::mult(MBARBAR, gmm::transposed(INT_Omega_Sigma) , Mbarbar_Sigma);
	gmm::copy(Mbarbar_Sigma, MBARBAR_Sigma);//gmm::copy(V, W);//W <-- V
		//M_OS4 = MBARBAR^T * Mvvg * Mbarbar_Sigma
	sparse_matrix_type M_partialg(dof_transp.Cv(), dof_transp.Ct_Sigma());gmm::clear(M_partialg);
	gmm::mult(Mvvg, Mbarbar_Sigma, M_partialg);
	gmm::mult(gmm::transposed(MBARBAR), M_partialg, M_OS4);
	gmm::mult(M_OS4, z1ref, l1_g2);	
int_l1_g2= gmm::vect_sp(Gt, l1_g2); cout <<"l1_g mean = "<<int_l1_g2<<endl;
	gmm::vscale(Gt, l1_g2);

	// l1_g = l1_g1 - l1_g2
	gmm::add(l1_g1, gmm::scaled(l1_g2,-1.0), l1_g);

 	// l1_g3 = (g,z1ref)_Sigma con asm_source term (localizzato su z)
	asm_source_term(l1_g3, mimt,  mf_Ct_Sigma, mf_Ct,Gt,descr_transp.SIGMA);

	vtk_export exp_lg3p(descr_transp.OUTPUT+"l1_g3_part.vtk");
	exp_lg3p.exporting(mf_Ct_Sigma);
	exp_lg3p.write_mesh();
	exp_lg3p.write_point_data(mf_Ct_Sigma, l1_g3, "l1_g3");

	gmm::vscale(z1ref,l1_g3);

	// l1_g4 = (g_bar,z1ref)_sigma (costruisci prima le medie)
	vector_type g_mean(dof_transp.Cv()); gmm::clear(g_mean);
	vector_type g_sigma(dof_transp.Ct_Sigma()); gmm::clear(g_sigma);
	gmm::mult(MBARBAR, Gt, g_mean);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, g_mean, g_sigma, 2);
	getfem::asm_mass_matrix(M_OS8, mimt, mf_Ct_Sigma, descr_transp.SIGMA);
	gmm::mult(M_OS8, g_sigma, l1_g4);
	cout <<"l1_g mean alternativo = "<<gmm::vect_sp(z1ref, l1_g4)<<endl;
	gmm::vscale(z1ref, l1_g4);

		
	} // end l1_g
		
// end l1



////////////// d1 operators
		/////////////// Grad	
	//d1_1 = int su Sigma dei gradienti di Ub(su Sigma) e z1ref(su Sigma)
	sparse_matrix_type A_Sigma(dof_transp.Ct_Sigma(),dof_transp.Ct_Sigma()); gmm::clear(A_Sigma);
	sparse_matrix_type Bgrad_Sigma(dof_transp.Ct_Sigma(),dof_transp.Ct_Sigma()); gmm::clear(Bgrad_Sigma);
	sparse_matrix_type A_Lambda(dof_transp.Cv(),dof_transp.Cv()); gmm::clear(A_Lambda);
	getfem::asm_stiffness_matrix_for_homogeneous_laplacian (A_Sigma, mimt, mf_Ct_Sigma);
	getfem::asm_stiffness_matrix_for_laplacian (A_Lambda, mimv, mf_Cv, mf_coefv, Area);
	
	gmm::mult(A_Sigma, Ub, d1_1);	
int_d1_1= gmm::vect_sp(z1ref, d1_1); cout <<"d1 grad (Sigma)= "<<int_d1_1<<endl;
	gmm::vscale(z1ref, d1_1);	

	//d1_2 = int su Lambda delle medie barbar dei gradienti (costruisci prima la matrice)
	sparse_matrix_type A_Lambdapart(dof_transp.Cv(), dof_transp.Ct_Sigma());gmm::clear(A_Lambdapart);
	sparse_matrix_type Mbarbar_Sigma(dof_transp.Cv(), dof_transp.Ct_Sigma());gmm::clear(Mbarbar_Sigma);
	gmm::mult(MBARBAR, gmm::transposed(INT_Omega_Sigma) , Mbarbar_Sigma);

	gmm::mult(A_Lambda, Mbarbar_Sigma, A_Lambdapart);
	gmm::mult(gmm::transposed(Mbarbar_Sigma), A_Lambdapart, Bgrad_Sigma);

	gmm::mult(Bgrad_Sigma, Ub, d1_2);	
int_d1_2= gmm::vect_sp(z1ref, d1_2); cout <<"d1 grad (barbar)= "<<int_d1_2<<endl;
	gmm::vscale(z1ref, d1_2);

	//d1_5 = int su Sigma delle medie barbar dei gradienti  (costruisci prima le medie)
	vector_type z_mean2(dof_transp.Cv()); gmm::clear(z_mean2);
	vector_type z_sigma2(dof_transp.Ct_Sigma()); gmm::clear(z_sigma2);	
	gmm::mult(MBARBAR_Sigma, z1ref, z_mean2);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, z_mean2, z_sigma2, 2);
	gmm::mult(A_Sigma, z_sigma2, d1_5);
	cout <<"d1_grad mean alternativo = "<<gmm::vect_sp(Ub, d1_5)<<endl;
	gmm::vscale(Ub, d1_5);

	// d1_7 = int su Lambda delle medie barbar dei gradienti (calcola le medie e poi fai l'integrale su lambda)s
	gmm::mult(A_Lambda, U, d1_7);	
	cout <<"d1 grad mean alternativo 2= "<<gmm::vect_sp(z_mean2, d1_7)<<endl;
	gmm::vscale(z_mean2, d1_7);	

	// d1_grad = d1_1 - d1-2
	gmm::add(d1_1, gmm::scaled(d1_2,-1.0), d1_grad);
int_d1_grad=int_d1_1-int_d1_2; cout <<"d1_grad tot= "<<int_d1_grad<<endl;

		/////////////// Flux
	//d1_3 = int su Sigma del flusso (matrice di massa) con Ub(su Sigma) e z1ref(su Sigma)
	sparse_matrix_type M_Sigma(dof_transp.Ct_Sigma(),dof_transp.Ct_Sigma()); gmm::clear(M_Sigma);
	sparse_matrix_type B_Sigma(dof_transp.Ct_Sigma(),dof_transp.Ct_Sigma()); gmm::clear(B_Sigma);
	sparse_matrix_type M_Lambda(dof_transp.Cv(),dof_transp.Cv()); gmm::clear(M_Lambda);
	getfem::asm_mass_matrix_param (M_Sigma, mimt, mf_Ct_Sigma, mf_coeft, Kt);
	getfem::asm_mass_matrix_param  (M_Lambda, mimv, mf_Cv, mf_coefv, Perimeter);
	
	gmm::mult(M_Sigma, Ub, d1_3);	
int_d1_3= gmm::vect_sp(z1ref, d1_3); cout <<"d1 flux (Sigma)= "<<int_d1_3<<endl;
	gmm::vscale(z1ref, d1_3);	

	//d1_4 = int su Lambda delle medie bar dei flussi (costruisci prima la matrice)
	sparse_matrix_type M_Lambdapart(dof_transp.Cv(), dof_transp.Ct_Sigma());gmm::clear(M_Lambdapart);
	gmm::mult(A_Lambda, Mbar_Sigma, M_Lambdapart);
	gmm::mult(gmm::transposed(Mbar_Sigma), M_Lambdapart, B_Sigma);

	gmm::mult(B_Sigma, Ub, d1_4);	
int_d1_4= gmm::vect_sp(z1ref, d1_4); cout <<"d1 flux (mean)= "<<int_d1_4<<endl;
	gmm::vscale(z1ref, d1_4);

	//d1_6 = int su Gamma delle medie bar dei flussi (costruisci prima le medie)
	vector_type z_mean(dof_transp.Cv()); gmm::clear(z_mean);
	vector_type z_sigma(dof_transp.Ct_Sigma()); gmm::clear(z_sigma);	
	gmm::mult(Mbar_Sigma, z1ref, z_mean);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, z_mean, z_sigma, 2);
	gmm::mult(M_Sigma, z_sigma, d1_6);
	cout <<"d1_flux mean alternativo = "<<gmm::vect_sp(Ub, d1_6)<<endl;
	gmm::vscale(Ub, d1_6);

	// d1_flux = d1_3 -d1_4
	gmm::add(d1_3, gmm::scaled(d1_4,-1.0), d1_flux);
int_d1_flux=int_d1_3-int_d1_4; cout <<"d1_flux tot= "<<int_d1_flux<<endl;

	// integrali finali per d
int_d1 = int_d1_grad + int_d1_flux;  cout <<"d1 tot= "<<int_d1<<endl;


int_l1_t=int_l1_t1-int_l1_t2; cout <<"l1_t tot= "<<int_l1_t<<endl;
int_l1_g=int_l1_g1-int_l1_g2; cout <<"l1_g tot= "<<int_l1_g<<endl;

int_j1=int_l1_t-int_l1_g; cout <<"j1= "<<int_j1<<endl;

// controlla che gli operatori bar e barbar eseguano correttamente la media.
	vector_type check_bar(dof_transp.Cv()); gmm::clear(check_bar);
	vector_type check_barbar(dof_transp.Cv()); gmm::clear(check_barbar);
	vector_type check_bar2(dof_transp.Cv()); gmm::clear(check_bar2);
	vector_type check_barbar2(dof_transp.Cv()); gmm::clear(check_barbar2);
	gmm::mult(Mbarbar_Sigma, Ub, check_barbar);
	gmm::add(U, gmm::scaled(check_barbar, -1), check_barbar);
	gmm::mult(Mbar_Sigma, Ub, check_bar);
	gmm::add(U, gmm::scaled(check_bar, -1), check_bar);
	cout<<"check_bar (should be zero):  "<<gmm::vect_norminf(check_bar)<<endl;
	cout<<"check_barbar (should be zero):  "<<gmm::vect_norminf(check_barbar)<<endl;	

	vtk_export exp_check_bar(descr_transp.OUTPUT+"check_bar.vtk");
	exp_check_bar.exporting(mf_Cv);
	exp_check_bar.write_mesh();
	exp_check_bar.write_point_data(mf_Cv, check_bar, "check_bar");

	vtk_export exp_check_barbar(descr_transp.OUTPUT+"check_barbar.vtk");
	exp_check_barbar.exporting(mf_Cv);
	exp_check_barbar.write_mesh();
	exp_check_barbar.write_point_data(mf_Cv, check_barbar, "check_barbar");

	// export l1 and d1
	vtk_export exp_lt1(descr_transp.OUTPUT+"l1_t1.vtk");
	exp_lt1.exporting(mf_Ct);
	exp_lt1.write_mesh();
	exp_lt1.write_point_data(mf_Ct, l1_t1, "l1_t1");

	vtk_export exp_lt2(descr_transp.OUTPUT+"l1_t2.vtk");
	exp_lt2.exporting(mf_Cv);
	exp_lt2.write_mesh();
	exp_lt2.write_point_data(mf_Cv, l1_t2, "l1_t2");

	vtk_export exp_lt3(descr_transp.OUTPUT+"l1_t3.vtk");
	exp_lt3.exporting(mf_Ct);
	exp_lt3.write_mesh();
	exp_lt3.write_point_data(mf_Ct, l1_t3, "l1_t3");

	vtk_export exp_lt4(descr_transp.OUTPUT+"l1_t4.vtk");
	exp_lt4.exporting(mf_Ct);
	exp_lt4.write_mesh();
	exp_lt4.write_point_data(mf_Ct, l1_t4, "l1_t4");

	vtk_export exp_lt5(descr_transp.OUTPUT+"l1_t5.vtk");
	exp_lt5.exporting(mf_Ct_Sigma);
	exp_lt5.write_mesh();
	exp_lt5.write_point_data(mf_Ct_Sigma, l1_t5, "l1_t5");

	vtk_export exp_lt6(descr_transp.OUTPUT+"l1_t6.vtk");
	exp_lt6.exporting(mf_Ct_Sigma);
	exp_lt6.write_mesh();
	exp_lt6.write_point_data(mf_Ct_Sigma, l1_t6, "l1_t6");

	vtk_export exp_lt7(descr_transp.OUTPUT+"l1_t7.vtk");
	exp_lt7.exporting(mf_Ct_Sigma);
	exp_lt7.write_mesh();
	exp_lt7.write_point_data(mf_Ct_Sigma, l1_t7, "l1_t7");

	vtk_export exp_lt8(descr_transp.OUTPUT+"l1_t8.vtk");
	exp_lt8.exporting(mf_Ct_Sigma);
	exp_lt8.write_mesh();
	exp_lt8.write_point_data(mf_Ct_Sigma, l1_t8, "l1_t8");

	vtk_export exp_lt(descr_transp.OUTPUT+"l1_t.vtk");
	exp_lt.exporting(mf_Ct);
	exp_lt.write_mesh();
	exp_lt.write_point_data(mf_Ct, l1_t, "l1_t");


	vtk_export exp_lg1(descr_transp.OUTPUT+"l1_g1.vtk");
	exp_lg1.exporting(mf_Ct);
	exp_lg1.write_mesh();
	exp_lg1.write_point_data(mf_Ct, l1_g1, "l1_g1");

	vtk_export exp_lg2(descr_transp.OUTPUT+"l1_g2.vtk");
	exp_lg2.exporting(mf_Ct);
	exp_lg2.write_mesh();
	exp_lg2.write_point_data(mf_Ct, l1_g2, "l1_g2");

	vtk_export exp_lg3(descr_transp.OUTPUT+"l1_g3.vtk");
	exp_lg3.exporting(mf_Ct_Sigma);
	exp_lg3.write_mesh();
	exp_lg3.write_point_data(mf_Ct_Sigma, l1_g3, "l1_g3");

	vtk_export exp_lg4(descr_transp.OUTPUT+"l1_g4.vtk");
	exp_lg4.exporting(mf_Ct_Sigma);
	exp_lg4.write_mesh();
	exp_lg4.write_point_data(mf_Ct_Sigma, l1_g4, "l1_g4");

	vtk_export exp_lg(descr_transp.OUTPUT+"l1_g.vtk");
	exp_lg.exporting(mf_Ct);
	exp_lg.write_mesh();
	exp_lg.write_point_data(mf_Ct, l1_g, "l1_g");

	vtk_export exp_za(descr_transp.OUTPUT+"z1refa.vtk");
	exp_za.exporting(mf_Ct);
	exp_za.write_mesh();
	exp_za.write_point_data(mf_Ct, z1refa, "z1refa");

	vtk_export exp_zb(descr_transp.OUTPUT+"z1refb.vtk");
	exp_zb.exporting(mf_Ct);
	exp_zb.write_mesh();
	exp_zb.write_point_data(mf_Ct, z1refb, "z1refb");

	vtk_export exp_Ua(descr_transp.OUTPUT+"Ua.vtk");
	exp_Ua.exporting(mf_Ct_Sigma);
	exp_Ua.write_mesh();
	exp_Ua.write_point_data(mf_Ct_Sigma, Ua, "Ua");

	vtk_export exp_Ub(descr_transp.OUTPUT+"Ub.vtk");
	exp_Ub.exporting(mf_Ct_Sigma);
	exp_Ub.write_mesh();
	exp_Ub.write_point_data(mf_Ct_Sigma, Ub, "Ub");

	vtk_export exp_Uc(descr_transp.OUTPUT+"Uc.vtk");
	exp_Uc.exporting(mf_Ct);
	exp_Uc.write_mesh();
	exp_Uc.write_point_data(mf_Ct, Uc, "Uc");

	vtk_export exp_Ud(descr_transp.OUTPUT+"Ud.vtk");
	exp_Ud.exporting(mf_Ct_Sigma);
	exp_Ud.write_mesh();
	exp_Ud.write_point_data(mf_Ct_Sigma, Ud, "Ud");

	vtk_export exp_gt(descr_transp.OUTPUT+"Gt.vtk");
	exp_gt.exporting(mf_Ct);
	exp_gt.write_mesh();
	exp_gt.write_point_data(mf_Ct, Gt, "Gt");

	vtk_export exp_gv(descr_transp.OUTPUT+"Gv.vtk");
	exp_gv.exporting(mf_Cv);
	exp_gv.write_mesh();
	exp_gv.write_point_data(mf_Cv, G, "Gv");


	vtk_export exp_d1_1(descr_transp.OUTPUT+"d1_1.vtk");
	exp_d1_1.exporting(mf_Ct_Sigma);
	exp_d1_1.write_mesh();
	exp_d1_1.write_point_data(mf_Ct_Sigma, d1_1, "d1_1");

	vtk_export exp_d1_2(descr_transp.OUTPUT+"d1_2.vtk");
	exp_d1_2.exporting(mf_Ct_Sigma);
	exp_d1_2.write_mesh();
	exp_d1_2.write_point_data(mf_Ct_Sigma, d1_2, "d1_2");

	vtk_export exp_d1_3(descr_transp.OUTPUT+"d1_3.vtk");
	exp_d1_3.exporting(mf_Ct_Sigma);
	exp_d1_3.write_mesh();
	exp_d1_3.write_point_data(mf_Ct_Sigma, d1_3, "d1_3");

	vtk_export exp_d1_4(descr_transp.OUTPUT+"d1_4.vtk");
	exp_d1_4.exporting(mf_Ct_Sigma);
	exp_d1_4.write_mesh();
	exp_d1_4.write_point_data(mf_Ct_Sigma, d1_4, "d1_4");

	vtk_export exp_d1_5(descr_transp.OUTPUT+"d1_5.vtk");
	exp_d1_5.exporting(mf_Ct_Sigma);
	exp_d1_5.write_mesh();
	exp_d1_5.write_point_data(mf_Ct_Sigma, d1_5, "d1_5");

	vtk_export exp_d1_6(descr_transp.OUTPUT+"d1_6.vtk");
	exp_d1_6.exporting(mf_Ct_Sigma);
	exp_d1_6.write_mesh();
	exp_d1_6.write_point_data(mf_Ct_Sigma, d1_6, "d1_6");

	vtk_export exp_d1_7(descr_transp.OUTPUT+"d1_7.vtk");
	exp_d1_7.exporting(mf_Cv);
	exp_d1_7.write_mesh();
	exp_d1_7.write_point_data(mf_Cv, d1_7, "d1_7");

	vtk_export exp_d1_grad(descr_transp.OUTPUT+"d1_grad.vtk");
	exp_d1_grad.exporting(mf_Ct_Sigma);
	exp_d1_grad.write_mesh();
	exp_d1_grad.write_point_data(mf_Ct_Sigma, d1_grad, "d1_grad");

	vtk_export exp_d1_flux(descr_transp.OUTPUT+"d1_flux.vtk");
	exp_d1_flux.exporting(mf_Ct_Sigma);
	exp_d1_flux.write_mesh();
	exp_d1_flux.write_point_data(mf_Ct_Sigma, d1_flux, "d1_flux");

// De-allocate memory
	gmm::clear(l1_t);
	gmm::clear(l1_t1);
	gmm::clear(l1_t2);
	gmm::clear(l1_t3);
	gmm::clear(l1_t4);
	gmm::clear(l1_t5);

	gmm::clear(l1_g);
	gmm::clear(l1_g1);
	gmm::clear(l1_g2);
	gmm::clear(l1_g3);

	gmm::clear(U);
	gmm::clear(Ua);
	gmm::clear(Ub);
	gmm::clear(Uc);
	gmm::clear(u);
	gmm::clear(z1ref);
	gmm::clear(z1refa);
	gmm::clear(z1refb);

	gmm::clear(M_OS);
	gmm::clear(M_OS2);
	gmm::clear(M_OS3);
	gmm::clear(M_OS4);

	gmm::clear(Mvv);
	gmm::clear(M_partial);
	gmm::clear(Mbar_Sigma);
	gmm::clear(INT_Omega_Sigma);

	gmm::clear(G);
	gmm::clear(Gt);
	
	};

	//! Compute model error of A2
	void transport3d1d::compute_model_error_2_test(){

	#ifdef M3D1D_VERBOSE_
	cout << "  Computing model error for assumption A2 ..." << endl;
	#endif
	// Storing error model from bilinear form 
	vector_type d2_1(dof_transp.Ct()); gmm::clear(d2_1); // z2ref
	vector_type d2_2(dof_transp.Ct()); gmm::clear(d2_2); // z2ref interpolato su tutto Omega
	vector_type d2_3(dof_transp.Ct()); gmm::clear(d2_3); // z3ref = z2red
	// Storing error model from linear form 
	vector_type l2_1(dof_transp.Ct()); gmm::clear(l2_1); // z2ref
	vector_type l2_2(dof_transp.Ct()); gmm::clear(l2_2); // z2ref interpolato su tutto Omega
	vector_type l2_3(dof_transp.Ct()); gmm::clear(l2_3); // z3ref = z2red
	// Solution of primal 3D problem
	vector_type u(dof_transp.Ct()); gmm::clear(u);
	// Solution of dual 3D problem
	vector_type z2ref(dof_transp.Ct_Omega()); gmm::clear(z2ref);
	vector_type z3ref(dof_transp.Ct()); gmm::clear(z3ref);

	gmm::copy(gmm::sub_vector(U3, 
		gmm::sub_interval(0, dof_transp.Ct())), u);
	gmm::copy(gmm::sub_vector(Z1, 
		gmm::sub_interval(0, dof_transp.Ct_Omega())), z2ref);
	gmm::copy(gmm::sub_vector(Z2, 
		gmm::sub_interval(0, dof_transp.Ct())), z3ref);

	scalar_type int_d2_1=0, int_d2_2=0;
	scalar_type int_l2_1=0, int_l2_2=0;
	scalar_type int_j2=0;

	vector_type z2refa(dof_transp.Ct()); gmm::clear(z2refa);
	getfem::interpolation(mf_Ct_Omega, mf_Ct, z2ref, z2refa, 2);
/////////////////// d2
	sparse_matrix_type At1(dof_transp.Ct(), dof_transp.Ct_Omega());gmm::clear(At1);
	sparse_matrix_type At2(dof_transp.Ct(), dof_transp.Ct());gmm::clear(At2);
// uso z2 originale, definito solo su Omega+.

    getfem::generic_assembly
      assem("M$1(#1,#2) += comp( Grad(#1).Grad(#2) )(:, :);");
    assem.push_mi(mimt);
    assem.push_mf(mf_Ct);
    assem.push_mf(mf_Ct_Omega);
    assem.push_mat(At1);
    assem.assembly(descr_transp.SIGMA);
	gmm::mult(At1, z2ref, d2_1);	
	gmm::vscale(u,d2_1);
// uso z2ref ma definito su tutto Omega( interpola)

	getfem::asm_stiffness_matrix_for_homogeneous_laplacian(At2, mimt, mf_Ct, descr_transp.SIGMA);
	gmm::mult(At2, z2refa, d2_2);	

int_d2_1= gmm::vect_sp(u, d2_2); cout <<"d2 interp = "<<int_d2_1<<endl;
	gmm::vscale(u,d2_2);

// uso z3ref
	gmm::mult(At2, z3ref, d2_3);	

int_d2_2= gmm::vect_sp(u, d2_3); cout <<"d2 reduced = "<<int_d2_2<<endl;
	gmm::vscale(u,d2_3);
///////////// l2
	// Source for primal 3D problem
	scalar_type f = PARAM.real_value("F", "Value of source for tissue reduced problem");	
	vector_type F(dof.coeft()); gmm::clear(F);
	expr = PARAM.string_value("F_EXPR", "Expression for source term in tissue (muparser)");
	interpolation_function(mf_coeft, F, f_function ); 
// uso z2 originale, definito solo su Omega+.
	sparse_matrix_type M_OS (dof_transp.Ct(),dof_transp.Ct_Omega()); gmm::clear(M_OS);
	getfem::asm_mass_matrix(M_OS, mimt, mf_Ct, mf_Ct_Omega, descr_transp.SIGMA);
	gmm::mult(M_OS, z2ref, l2_1);	
	gmm::vscale(u,l2_1);
// uso z2ref ma definito su tutto Omega( interpola)
	sparse_matrix_type M_OS2 (dof_transp.Ct(),dof_transp.Ct()); gmm::clear(M_OS2);
	getfem::asm_mass_matrix(M_OS2, mimt, mf_Ct, mf_Ct, descr_transp.SIGMA);
	gmm::mult(M_OS2, z2refa, l2_2);	
int_l2_1= gmm::vect_sp(u, l2_2); cout <<"l2 interp = "<<int_l2_1<<endl;
	gmm::vscale(u,l2_2);

// uso z3ref
	gmm::mult(M_OS2, z3ref, l2_3);	
int_l2_2= gmm::vect_sp(u, l2_3); cout <<"l2 reduced = "<<int_l2_2<<endl;
	gmm::vscale(u,l2_3);



cout <<"j2_1= "<<int_d2_1-int_l2_1<<endl;
cout <<"j2_2= "<<int_d2_2-int_l2_2<<endl;
///////////

	
	//export d2 and l2

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting d2 and l2 ..." << endl;
	#endif
	vtk_export exp_d2_1(descr_transp.OUTPUT+"d2_1.vtk");
	exp_d2_1.exporting(mf_Ct);
	exp_d2_1.write_mesh();
	exp_d2_1.write_point_data(mf_Ct, d2_1, "d2_1");

	vtk_export exp_d2_2(descr_transp.OUTPUT+"d2_2.vtk");
	exp_d2_2.exporting(mf_Ct);
	exp_d2_2.write_mesh();
	exp_d2_2.write_point_data(mf_Ct, d2_2, "d2_2");

	vtk_export exp_d2_3(descr_transp.OUTPUT+"d2_3.vtk");
	exp_d2_3.exporting(mf_Ct);
	exp_d2_3.write_mesh();
	exp_d2_3.write_point_data(mf_Ct, d2_3, "d2_3");

	vtk_export exp_l2_1(descr_transp.OUTPUT+"l2_1.vtk");
	exp_l2_1.exporting(mf_Ct);
	exp_l2_1.write_mesh();
	exp_l2_1.write_point_data(mf_Ct, l2_1, "l2_1");

	vtk_export exp_l2_2(descr_transp.OUTPUT+"l2_2.vtk");
	exp_l2_2.exporting(mf_Ct);
	exp_l2_2.write_mesh();
	exp_l2_2.write_point_data(mf_Ct, l2_2, "l2_2");

	vtk_export exp_l2_3(descr_transp.OUTPUT+"l2_3.vtk");
	exp_l2_3.exporting(mf_Ct);
	exp_l2_3.write_mesh();
	exp_l2_3.write_point_data(mf_Ct, l2_3, "l2_3");

	vtk_export exp_z2ref(descr_transp.OUTPUT+"z2ref.vtk");
	exp_z2ref.exporting(mf_Ct_Omega);
	exp_z2ref.write_mesh();
	exp_z2ref.write_point_data(mf_Ct_Omega, z2ref, "z2ref");

	vtk_export exp_z2refa(descr_transp.OUTPUT+"z2refa.vtk");
	exp_z2refa.exporting(mf_Ct);
	exp_z2refa.write_mesh();
	exp_z2refa.write_point_data(mf_Ct, z2refa, "z2refa");

	vtk_export exp_z3ref(descr_transp.OUTPUT+"z3ref.vtk");
	exp_z3ref.exporting(mf_Ct);
	exp_z3ref.write_mesh();
	exp_z3ref.write_point_data(mf_Ct, z3ref, "z3ref");

	};


	//! Compute model error of A3	
	void transport3d1d::compute_model_error_3_test(){
	#ifdef M3D1D_VERBOSE_
	cout << "  Computing model error for assumption A3 ..." << endl;
	#endif

	// Storing error model from bilinear form 
	vector_type d3(dof_transp.Ct()); gmm::clear(d3);
	vector_type d3_1(dof_transp.Ct()); gmm::clear(d3_1); // d3_1 = (k u,   z3ref)Gamma con asm mass matrix
	vector_type d3_2(dof_transp.Ct()); gmm::clear(d3_2); // d3_2 = (k ubar,z3ref)Lambda costruendo la matrice di scambio
	vector_type d3_3(dof_transp.Ct_Sigma()); gmm::clear(d3_3);// d3_3 = (k ubar,z3ref)gamma costruendo le medie

	//Storing error model from linear form 
	vector_type l3(dof_transp.Ct()); gmm::clear(l3); 
	vector_type l3_1(dof_transp.Ct()); gmm::clear(l3_1); // l3_1 = (k U, z3ref)Gamma con asm mass matrix
	vector_type l3_2(dof_transp.Ct()); gmm::clear(l3_2); // l3_2 = (k U, z3refbar)Lambda costruendo la matrice di scambio
	vector_type l3_3(dof_transp.Ct_Sigma()); gmm::clear(l3_3);// l3_3 = (k U, z3refbar)gamma costruendo le medie
	vector_type l3_4(dof_transp.Ct()); gmm::clear(l3_4);
	//cout<<"WARNING: l3(u,z)=0!!"<<endl;

	// integrali globali
	scalar_type int_d3_1=0, int_d3_2=0, int_d3=0;
	scalar_type int_l3_1=0, int_l3_2=0, int_l3=0;
	scalar_type int_j3=0;

	// Solution of primal 3D problem
	vector_type u(dof_transp.Ct()); gmm::clear(u);

	// Solution of primal 1D problem
	vector_type U(dof_transp.Cv()); gmm::clear(U);
	vector_type Ua(dof_transp.Ct_Sigma()); gmm::clear(Ua);
	vector_type Ub(dof_transp.Ct_Sigma()); gmm::clear(Ub);
	vector_type Uc(dof_transp.Ct()); gmm::clear(Uc);

	// Solution of dual 3D problem
	vector_type z3ref(dof_transp.Ct()); gmm::clear(z3ref);

	// Import data
	gmm::copy(gmm::sub_vector(U3, 
		gmm::sub_interval(0, dof_transp.Ct())), u);
	gmm::copy(gmm::sub_vector(Z2, 
		gmm::sub_interval(0, dof_transp.Ct())), z3ref);

	gmm::copy(gmm::sub_vector(U3, 
				  gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv())), U);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, U, Ub, 2);

	//parameters
	scalar_type k = PARAM.real_value("kk", "value of permeability for reduced problem");
	// Permeability
	vector_type Kv(dof.coefv());
	vector_type Kt(dof.coeft());
	constant_value=k;
	interpolation_function(mf_coefv, Kv, constant_function ); 
	interpolation_function(mf_coeft, Kt, constant_function ); 
	//Perimeter = 2*pi*R*k 
	vector_type Perimeter(dof.coefv());
	gmm::copy(param.R(), Perimeter);
	gmm::scale(Perimeter, 2*pi*k);
	//Area = pi*R^2
	vector_type Area(dof.coefv());
	gmm::copy(param.R(), Area);
	gmm::vscale(param.R(), Area);
	gmm::scale(Area, pi); 

///////////////// d3
	// d3_1 = (k u,   z3ref)Gamma con asm mass matrix
sparse_matrix_type M_OS (dof_transp.Ct(),dof_transp.Ct()); gmm::clear(M_OS);
	getfem::asm_mass_matrix_param(M_OS, mimt, mf_Ct, mf_Ct, mf_coeft, Kt, descr_transp.GAMMA);
	gmm::scale(M_OS, 0.5);
	gmm::mult(M_OS, u, d3_1);	
int_d3_1= gmm::vect_sp(z3ref, d3_1); cout <<"d3 mass = "<<int_d3_1<<endl;
	gmm::vscale(z3ref,d3_1);

	// d3_2 = (k ubar,z3ref)Lambda costruendo la matrice di scambio
	sparse_matrix_type Btt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Btt);
	sparse_matrix_type Mvv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Mvv);
	sparse_matrix_type M_partial(dof_transp.Cv(), dof_transp.Ct());gmm::clear(M_partial);
	getfem::asm_mass_matrix_param (Mvv, mimv, mf_Cv,mf_coefv, Perimeter);
	gmm::mult(Mvv, MBAR, M_partial);
	gmm::mult(gmm::transposed(MBAR), M_partial, Btt);
	gmm::mult(Btt, u, d3_2);	
int_d3_2= gmm::vect_sp(z3ref, d3_2); cout <<"d3 mean = "<<int_d3_2<<endl;
	gmm::vscale(z3ref,d3_2);


	// d3_3 = (k ubar,z3ref)gamma costruendo le medie
	vector_type u_mean(dof_transp.Cv()); gmm::clear(u_mean);
	vector_type u_sigma(dof_transp.Ct_Sigma()); gmm::clear(u_sigma);
	gmm::mult(MBAR, u, u_mean);
	sparse_matrix_type M_omega (dof_transp.Ct_Sigma(),dof_transp.Ct()); gmm::clear(M_omega);
	getfem::asm_mass_matrix_param(M_omega, mimt, mf_Ct_Sigma, mf_Ct, mf_coeft, Kt, descr_transp.GAMMA);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, u_mean, u_sigma, 2);
	gmm::mult(M_omega,z3ref , d3_3);
	cout <<"d3 mean alternativo = "<<gmm::vect_sp(u_sigma, d3_3)<<endl;
	gmm::vscale(u_sigma, d3_3);

	gmm::add(d3_1, gmm::scaled(d3_2,-1.0), d3);
	cout <<"d3= "<<int_d3_1-int_d3_2<<endl;
//////////////////// l3
	// l3_1 = (k U, z3ref)Gamma con asm mass matrix
	sparse_matrix_type M_OS2 (dof_transp.Ct(), dof_transp.Ct_Sigma()); gmm::clear(M_OS2);
	sparse_matrix_type M_OS3 (dof_transp.Ct(), dof_transp.Ct_Sigma()); gmm::clear(M_OS3);
	getfem::asm_mass_matrix_param(M_OS2, mimt, mf_Ct, mf_Ct_Sigma, mf_coeft, Kt, descr_transp.GAMMA);
	gmm::mult(M_OS2, Ub, l3_1);	
int_l3_1= gmm::vect_sp(z3ref, l3_1); cout <<"l3 mass = "<<int_l3_1<<endl;
	gmm::vscale(z3ref, l3_1);


	// l3_2 = (k U, z3refbar)Lambda costruendo la matrice di scambio	
	sparse_matrix_type MvvS(dof_transp.Cv(), dof_transp.Cv());gmm::clear(MvvS);
	sparse_matrix_type M_partialS(dof_transp.Cv(), dof_transp.Ct_Sigma());gmm::clear(M_partialS);
	getfem::asm_mass_matrix_param (MvvS, mimv, mf_Cv,mf_coefv, Perimeter);
	sparse_matrix_type Mbar_Sigma(dof_transp.Cv(), dof_transp.Ct_Sigma());gmm::clear(Mbar_Sigma);
	sparse_matrix_type INT_Omega_Sigma(dof_transp.Ct_Sigma(), dof_transp.Ct());gmm::clear(INT_Omega_Sigma);
 	getfem::interpolation(mf_Ct, mf_Ct_Sigma, INT_Omega_Sigma, 2);
	gmm::mult(MBAR, gmm::transposed(INT_Omega_Sigma) , Mbar_Sigma);
	gmm::mult(MvvS, Mbar_Sigma, M_partialS);
	gmm::mult(gmm::transposed(MBAR), M_partialS, M_OS3);
	gmm::mult(M_OS3, Ub, l3_2);	
int_l3_2= gmm::vect_sp(z3ref, l3_2); cout <<"l3 mean = "<<int_l3_2<<endl;
	gmm::vscale(z3ref, l3_2);

	// l3_3 = (k U, z3refbar)gamma costruendo le medie
	vector_type z_mean(dof_transp.Cv()); gmm::clear(z_mean);
	vector_type z_sigma(dof_transp.Ct_Sigma()); gmm::clear(z_sigma);
	gmm::mult(MBAR, z3ref, z_mean);
	getfem::interpolation(mf_Cv, mf_Ct_Sigma, z_mean, z_sigma, 2);
	sparse_matrix_type M_sigma (dof_transp.Ct_Sigma(),dof_transp.Ct_Sigma()); gmm::clear(M_sigma);
	getfem::asm_mass_matrix_param(M_sigma, mimt, mf_Ct_Sigma, mf_Ct_Sigma, mf_coeft, Kt, descr_transp.GAMMA);
	gmm::mult(M_sigma, z_sigma, l3_3);
	cout <<"l3 mean alternativo = "<<gmm::vect_sp(Ub, l3_3)<<endl;
	gmm::vscale(Ub, l3_3);

	gmm::add(l3_1, gmm::scaled(l3_2,-1.0), l3);
	cout <<"l3= "<<int_l3_1-int_l3_2<<endl;

	// l3_4 = (k U, z3refbar)Lambda costruendo la matrice di scambio (U in Omega)	
	vector_type U_omega(dof_transp.Ct()); gmm::clear(U_omega);
	getfem::interpolation(mf_Cv, mf_Ct, U, U_omega, 2);
	sparse_matrix_type M_partialS2(dof_transp.Cv(), dof_transp.Ct());gmm::clear(M_partialS2);
	sparse_matrix_type M_OS4 (dof_transp.Ct(), dof_transp.Ct()); gmm::clear(M_OS4);
	gmm::mult(MvvS, MBAR, M_partialS2);
	gmm::mult(gmm::transposed(MBAR), M_partialS2, M_OS4);
	gmm::mult(M_OS4, U_omega, l3_4);	
	cout <<"l3 mean alternativo 2= "<< gmm::vect_sp(z3ref, l3_4) <<endl;
	gmm::vscale(z3ref, l3_4);

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting d3 and l3 ..." << endl;
	#endif
	vtk_export exp_d(descr_transp.OUTPUT+"d3.vtk");
	exp_d.exporting(mf_Ct);
	exp_d.write_mesh();
	exp_d.write_point_data(mf_Ct, d3, "d3");

	vtk_export exp_d3_1(descr_transp.OUTPUT+"d3_1.vtk");
	exp_d3_1.exporting(mf_Ct);
	exp_d3_1.write_mesh();
	exp_d3_1.write_point_data(mf_Ct, d3_1, "d3_1");

	vtk_export exp_d3_2(descr_transp.OUTPUT+"d3_2.vtk");
	exp_d3_2.exporting(mf_Ct);
	exp_d3_2.write_mesh();
	exp_d3_2.write_point_data(mf_Ct, d3_2, "d3_2");

	vtk_export exp_d3_3(descr_transp.OUTPUT+"d3_3.vtk");
	exp_d3_3.exporting(mf_Ct_Sigma);
	exp_d3_3.write_mesh();
	exp_d3_3.write_point_data(mf_Ct_Sigma, d3_3, "d3_3");

	vtk_export exp_l(descr_transp.OUTPUT+"l3.vtk");
	exp_l.exporting(mf_Ct);
	exp_l.write_mesh();
	exp_l.write_point_data(mf_Ct, l3, "l3");

	vtk_export exp_l3_1(descr_transp.OUTPUT+"l3_1.vtk");
	exp_l3_1.exporting(mf_Ct);
	exp_l3_1.write_mesh();
	exp_l3_1.write_point_data(mf_Ct, l3_1, "l3_1");

	vtk_export exp_l3_2(descr_transp.OUTPUT+"l3_2.vtk");
	exp_l3_2.exporting(mf_Ct);
	exp_l3_2.write_mesh();
	exp_l3_2.write_point_data(mf_Ct, l3_2, "l3_2");

	vtk_export exp_l3_3(descr_transp.OUTPUT+"l3_3.vtk");
	exp_l3_3.exporting(mf_Ct_Sigma);
	exp_l3_3.write_mesh();
	exp_l3_3.write_point_data(mf_Ct_Sigma, l3_3, "l3_3");

	vtk_export exp_l3_4(descr_transp.OUTPUT+"l3_4.vtk");
	exp_l3_4.exporting(mf_Ct);
	exp_l3_4.write_mesh();
	exp_l3_4.write_point_data(mf_Ct, l3_4, "l3_4");

	};

*/

} // end of namespace
