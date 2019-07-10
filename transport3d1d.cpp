/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   transport3d1d.cpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016 - May 2018
  @brief  Definition of the main class for the 3D/1D coupled transport problem.
 */
 
 #include <transport3d1d.hpp>
 #include <AMG_Interface.hpp>
 #include <cmath>
 #include "gmm/gmm_inoutput.h"
 #include "getfem/getfem_import.h"


 namespace getfem {


	// Initialize the transport problem
 	void transport3d1d::init_transp (int argc, char *argv[]) 
 	{
 	#ifdef M3D1D_VERBOSE_
 	std::cout << "initialize transport problem..."<<std::endl<<std::endl;
 	#endif

	PARAM.read_command_line(argc, argv);
	//1. Import data (algorithm specifications, boundary conditions, ...)	
	import_data_transp();
	//2. Import mesh for tissue (3D) and vessel network (1D)
	build_mesh_transp();
	//3. Set finite elements and integration methods
	set_im_and_fem_transp();
 	//4. Build problem parameters
 	build_param_transp();
	//5. Build the list of tissue boundary data
	build_tissue_boundary_transp();
	//6. Build the list of tissue boundary (and junction) data
 	build_vessel_boundary_transp();

 	}; // end of init


 	// Aux methods for init
	
	// Import algorithm specifications
	void transport3d1d::import_data_transp(void)
	{
	#ifdef M3D1D_VERBOSE_
	cout << "Importing descriptors for tissue and vessel problems ..." << endl;
	#endif

	descr.import(PARAM);
	descr_transp.import(PARAM);
	#ifdef M3D1D_VERBOSE_
	cout << descr_transp;
	#endif
	 
 
	}; //end of import_data_transp
	
	 
	
	// Import mesh for tissue (3D) and vessel (1D)  
	void transport3d1d::build_mesh_transp(void){

		//In order to have the boundary conditions for the nodes
		//we need to build again the 1D mesh from another pts file
		//! \todo write import_msh_file_transp which take in account both .pts files: modify transport3d1d::init_fluid such that it calls import_msh_file_transp.
		//! \todo branch index should start from 1, not from 0: the junctions will use the notation that +branch is inflow and -branch is outflow. Branch 0 force the user to not give in the .pts file an outflow for first branch (this is usually ok, but loses generality)
		mesht.clear();
			bool test = 0;
		test = PARAM.int_value("TEST_GEOMETRY");
		if(test==0){
			#ifdef M3D1D_VERBOSE_
			cout << "Importing the 3D mesh for the tissue ...  "   << endl;
			#endif
			string st("gmsh:"+descr.MESH_FILET);
			getfem::import_mesh(st,mesht);
				 //import_msh_file(descr.MESH_FILET, mesht);
			}
		else{
			#ifdef M3D1D_VERBOSE_
			cout << "Building the regular 3D mesh for the tissue ...  "   << endl;
			#endif
			string st("GT='" + PARAM.string_value("GT_T") + "'; " +
						   "NSUBDIV=" + PARAM.string_value("NSUBDIV_T") + "; " +  
						   "ORG=" + PARAM.string_value("ORG_T") + "; " +  
						   "SIZES=" + PARAM.string_value("SIZES_T") + "; " +  
						   "NOISED=" + PARAM.string_value("NOISED_T")); 
			#ifdef M3D1D_VERBOSE_		
			cout << "mesht description: " << st << endl;
			#endif
			regular_mesh(problem3d1d::mesht, st);
			}
	
	#ifdef M3D1D_VERBOSE_
	cout << "Importing the 1D mesh for the vessel (transport problem)... "   << endl;
	#endif
	std::ifstream ifs(descr_transp.MESH_FILEV);
	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descr_transp.MESH_FILEV);
	meshv.clear();
	bool Import=PARAM.int_value("IMPORT_CURVE");
	bool Curve=PARAM.int_value("CURVE_PROBLEM");

	if(Curve && !Import){
		import_pts_file(ifs, meshv, BCv_transp, nb_vertices, descr.MESH_TYPEV, param);
	}
	
	else if(Import && !Curve){
		GMM_ASSERT1(0,"If you want to import the curvature, you need to enable CURVE_PROBLEM=1");
	}
	else if(Import && Curve){
		std::ifstream ifc(PARAM.string_value("CURVE_FILE","curvature file location"));
		GMM_ASSERT1(ifc.good(), "impossible to read from file " << PARAM.string_value("CURVE_FILE","curvature file location"));
		
		import_pts_file(ifs,ifc, meshv, BCv_transp, nb_vertices, descr.MESH_TYPEV, param);

		ifc.close();
	} 
	else{
		import_pts_file(ifs, meshv, BCv_transp, nb_vertices, descr.MESH_TYPEV);
	}


	nb_branches = nb_vertices.size();
	ifs.close();

	
	}; // end of build_mesh_geometry

	// Set finite elements methods and integration methods 
	void transport3d1d::set_im_and_fem_transp(void)
	{

	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEMs for tissue and vessel problems ..." << endl;
	#endif
	
	
	pfem pf_Ct = fem_descriptor(descr_transp.FEM_TYPET_C);
	pfem pf_Cv = fem_descriptor(descr_transp.FEM_TYPEV_C);

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for tissue ..." << endl;
	#endif
		

	mf_Ct.set_finite_element(mesht.convex_index(), pf_Ct);

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for vessel branches ..." << endl;
	#endif

	mf_Cv.set_finite_element(meshv.convex_index(), pf_Cv);

	
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEM dimensions for tissue and vessel problems ..." << endl;
	#endif

	dof_transp.set(mf_Ct, mf_Cv);
	#ifdef M3D1D_VERBOSE_
	cout << std::scientific << dof_transp;
	#endif
	

	// I had to delete the meshes and build them again, because of issues on the boundary conditions:
	// Now i have to build again also the mesh_fem and mesh_im objects.
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

	}; // end of set_im_and_fem
	
	
	// Build problem parameters
	void transport3d1d::build_param_transp(void)
	{
	
	#ifdef M3D1D_VERBOSE_
	cout << "Building parameters for tissue and vessel problems ..." << endl;
	#endif
	param.build(PARAM, mf_coeft, mf_coefv,mf_coefvi); 
	param_transp.build(PARAM, mf_coeft, mf_coefv);
	#ifdef M3D1D_VERBOSE_
	cout << param_transp ;
	#endif

	}; // end of build_param_transp
  
  
  	//Build boundary regions on tissue
	void
	transport3d1d::build_tissue_boundary_transp (void) 
	{
	#ifdef M3D1D_VERBOSE_
	cout << "Building tissue boundary ..." << endl;
	#endif

	size_type face=descr_transp.FACE;

	BCt_transp.clear();
	BCt_transp.reserve(2*DIMT);
	// Parse BC data
	string label_in = PARAM.string_value("BClabel_transp", "Array of tissue boundary labels");
	string value_in = PARAM.string_value("BCvalue_transp", "Array of tissue boundary values");
	vector<string> labels = split(label_in, ' ');
	vector<string> values = split(value_in, ' ');
	GMM_ASSERT1(labels.size()==2*DIMT, "wrong number of BC labels");
	GMM_ASSERT1(values.size()==2*DIMT, "wrong number of BC values");
	for (unsigned f=0; f<2*DIMT; ++f) 
		{
			BCt_transp.emplace_back(labels[f], std::stof(values[f]), 0, face+f);
			#ifdef M3D1D_VERBOSE_
			cout << "  face " << f << " : " << BCt_transp.back() << endl;
			#endif
		} 
	

	// Build mesht regions
	size_type xx=0, yy=1, zz=2;

	mesh_region border_faces;
	outer_faces_of_mesh(mesht, border_faces);

	for (mr_visitor i(border_faces); !i.finished(); ++i) 
		{

		assert(i.is_face());

		// Unit outward normal : used to identify faces
		//! \todo Use getfem 5.0's function select_faces_of_normal?
		base_node un = mesht.normal_of_face_of_convex(i.cv(), i.f());
		un /= gmm::vect_norm2(un);

		     if (gmm::abs(un[xx] + 1.0) < 1.0E-7) 	// back
			mesht.region(face+0).add(i.cv(), i.f());
		else if (gmm::abs(un[xx] - 1.0) < 1.0E-7) 	// front
			mesht.region(face+1).add(i.cv(), i.f());
		else if (gmm::abs(un[yy] + 1.0) < 1.0E-7) 	// left
			mesht.region(face+2).add(i.cv(), i.f());
		else if (gmm::abs(un[yy] - 1.0) < 1.0E-7) 	// right
			mesht.region(face+3).add(i.cv(), i.f());
		else if (gmm::abs(un[zz] + 1.0) < 1.0E-7) 	// bottom
			mesht.region(face+4).add(i.cv(), i.f());
		else if (gmm::abs(un[zz] - 1.0) < 1.0E-7) 	// top
			mesht.region(face+5).add(i.cv(), i.f());
				
		}



	GMM_ASSERT1(!(PARAM.int_value("TEST_GEOMETRY") && descr_transp.CONFORMING), "Test 3D mesh cannot be conforming! Check the .param file!!!");
	
	if(descr_transp.CONFORMING){
	/*! \todo Regions SIGMA and GAMMA are described in .msh file from GMESH (using PhysicalVolume). 
                  Build directly those regions ("if(dist(element, meshv)>0){OMEGA.add(element)}") 
		  or build the indicator functions and rescale all the integrals by those functions.
	*/
	#ifdef M3D1D_VERBOSE_
	cout << "The mesh is conforming with respect to the vessel!" << endl;
	cout << "Checking the regions..." << endl;
	#endif

	size_type sigma=descr_transp.SIGMA;
	size_type omega=descr_transp.OMEGA;
	size_type gamma=descr_transp.GAMMA;	

	//check if there is a overlap in regions!
	GMM_ASSERT1(sigma!=gamma, "SIGMA=GAMMA: check the .param file");
	GMM_ASSERT1(omega!=sigma, "OMEGA=SIGMA: check the .param file");
	GMM_ASSERT1(gamma!=omega, "GAMMA=OMEGA: check the .param file");
	
	for(int i=0; i<5; i++){
	GMM_ASSERT1((face+i)!=gamma, "FACE=GAMMA: check the .param file");
	GMM_ASSERT1((face+i)!=sigma, "FACE=SIGMA: check the .param file");
	GMM_ASSERT1((face+i)!=omega, "FACE=OMEGA: check the .param file");
}

	//Check if sigma and omega are defined in the msh file
	GMM_ASSERT1(mesht.has_region(sigma), "File .msh does not contain region SIGMA! Check .msh and .param!!");
	GMM_ASSERT1(mesht.has_region(omega), "File .msh does not contain region OMEGA! Check .msh and .param!!");

	// GAMMA is the boundary of the vessel SIGMA:
	
	mesh_region gamma_region;
	outer_faces_of_mesh(mesht, mesht.region(omega), gamma_region);

	for (mr_visitor i(gamma_region); !i.finished(); ++i) {

		assert(i.is_face());

		     if (!(mesht.region(face+0).is_in(i.cv(),i.f()))&&!(mesht.region(face+1).is_in(i.cv(),i.f()))
				&&!(mesht.region(face+2).is_in(i.cv(),i.f()))&&!(mesht.region(face+3).is_in(i.cv(),i.f()))
				&&!(mesht.region(face+4).is_in(i.cv(),i.f()))&&!(mesht.region(face+5).is_in(i.cv(),i.f()))) 	// back
			mesht.region(gamma).add(i.cv(), i.f());
			mesht.region(gamma+1).add(i.cv(), i.f());
				
		} 

	mesh_region gamma_region2;
	outer_faces_of_mesh(mesht, mesht.region(sigma), gamma_region2);

	for (mr_visitor i(gamma_region2); !i.finished(); ++i) {

		assert(i.is_face());

		     if (!(mesht.region(face+0).is_in(i.cv(),i.f()))&&!(mesht.region(face+1).is_in(i.cv(),i.f()))
				&&!(mesht.region(face+2).is_in(i.cv(),i.f()))&&!(mesht.region(face+3).is_in(i.cv(),i.f()))
				&&!(mesht.region(face+4).is_in(i.cv(),i.f()))&&!(mesht.region(face+5).is_in(i.cv(),i.f()))) 	// back
			mesht.region(gamma).add(i.cv(), i.f());	
			mesht.region(gamma+2).add(i.cv(), i.f());
			//mesht.region(gamma+1).add(i.cv(), i.f());					
		}

	//outer_faces_of_mesh(mesht, mesht.region(sigma), mesht.region(gamma));
	GMM_ASSERT1(mesht.has_region(gamma), "Something wrong happened: I couldn't build the region GAMMA");

	outer_faces_of_mesh(mesht, mesht.region(omega), mesht.region(gamma+1));
	GMM_ASSERT1(mesht.has_region(gamma), "Something wrong happened: I couldn't build the region GAMMA+1");


	#ifdef M3D1D_VERBOSE_
	cout << "...Check complete! All 3D regions are correctly defined!" << endl;
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


/*	//The following are some tests for the different regions
	#ifdef M3D1D_VERBOSE_
	dal::bit_vector nn = mesht.convex_index();
	bgeot::size_type i;
	size_type cont_sigma=0;
	size_type cont_omega=0;
	size_type cont_gamma=0;
	size_type cont_gamma2=0;
	size_type cont_gamma_sigma=0;	
	size_type cont_gamma_omega=0;	
	size_type cont_sigma_omega=0;
	size_type cont_sigma_face=0;
	size_type cont_omega_face=0;
	size_type cont_gamma_sigma_face=0;	
	size_type cont_gamma_omega_face=0;
	size_type cont_gamma2_omega_face=0;
	size_type cont_gamma2_omega_face_ext=0;
	size_type cont_gamma_gamma2=0;


	for (i << nn; i != bgeot::size_type(-1); i << nn) {
		bgeot::pconvex_structure cvs = mesht.structure_of_convex(i);
		if(mesht.region(sigma).is_in(i)) cont_sigma++;
		if(mesht.region(omega).is_in(i)) cont_omega++;

		for (bgeot::short_type f = 0; f < cvs->nb_faces(); ++f) {
			if(mesht.region(gamma).is_in(i,f)) cont_gamma++;
			if(mesht.region(gamma+1).is_in(i,f)) cont_gamma2++;
			if(mesht.region(sigma).is_in(i,f)) cont_sigma_face++;
			if(mesht.region(omega).is_in(i,f)) cont_omega_face++;
			if( (mesht.region(gamma).is_in(i,f)) && (mesht.region(sigma).is_in(i,f)) ) cont_gamma_sigma++;
			if( (mesht.region(gamma).is_in(i,f)) && (mesht.region(omega).is_in(i,f)) ) cont_gamma_omega++;
			if( (mesht.region(sigma).is_in(i,f)) && (mesht.region(omega).is_in(i,f)) ) cont_sigma_omega++;
			if( (mesht.region(gamma).is_in(i,f)) && (mesht.region(sigma).is_in(i)) 
				&&!(mesht.region(face+0).is_in(i,f))&&!(mesht.region(face+1).is_in(i,f))
				&&!(mesht.region(face+2).is_in(i,f))&&!(mesht.region(face+3).is_in(i,f))
				&&!(mesht.region(face+4).is_in(i,f))&&!(mesht.region(face+5).is_in(i,f))) {				
				cont_gamma_sigma_face++;
//				cout<<"Element "<<i<<" is in Sigma with face "<<f<<" on gamma"<<endl;
}
			if( (mesht.region(gamma).is_in(i,f)) && (mesht.region(omega).is_in(i)) ) {
				cont_gamma_omega_face++;
//				cout<<"Element "<<i<<" is in Omega with face "<<f<<" on gamma"<<endl;
}
			if( (mesht.region(gamma+1).is_in(i,f)) && (mesht.region(omega).is_in(i)) ) {
				cont_gamma2_omega_face_ext++;
//				cout<<"Element "<<i<<" is in Omega2 with face "<<f<<" on gamma"<<endl;
}	
			if( (mesht.region(gamma+1).is_in(i,f)) && (mesht.region(omega).is_in(i))
			 	&&!(mesht.region(face+0).is_in(i,f))&&!(mesht.region(face+1).is_in(i,f))
				&&!(mesht.region(face+2).is_in(i,f))&&!(mesht.region(face+3).is_in(i,f))
				&&!(mesht.region(face+4).is_in(i,f))&&!(mesht.region(face+5).is_in(i,f)) ) {
				cont_gamma2_omega_face++;
//				cout<<"Element "<<i<<" is in Omega2 with face "<<f<<" on gamma (no exterior)"<<endl;
}	
			if( (mesht.region(gamma).is_in(i,f)) && (mesht.region(gamma+1).is_in(i,f)) ) cont_gamma_gamma2++;
		}
	}
	cout <<"Number of element in Sigma: "<<cont_sigma<<endl;
	cout <<"Number of element in Omega: "<<cont_omega<<endl;
	cout <<"Number of faces in Gamma: "<<cont_gamma<<endl;
	cout <<"Number of faces in Gamma2 (including exterior faces): "<<cont_gamma2<<endl;
	cout <<"Number of faces both in Gamma and Sigma: "<<cont_gamma_sigma<<endl;
	cout <<"Number of faces both in Gamma and Omega: "<<cont_gamma_omega<<endl;
	cout <<"Number of faces both in Omega and Sigma: "<<cont_sigma_omega<<endl;
	cout <<"Number of faces both in Gamma and Gamma2: "<<cont_gamma_gamma2<<endl;
	cout <<"Number of faces in Sigma: "<<cont_sigma_face<<endl;
	cout <<"Number of faces in Omega: "<<cont_omega_face<<endl;
	cout <<"Number of elements in Sigma with faces in Gamma: "<<cont_gamma_sigma_face<<endl;
	cout <<"Number of elements in Omega with faces in Gamma2: "<<cont_gamma2_omega_face_ext<<endl;
	cout <<"Number of elements in Omega with faces in Gamma2(no exterior): "<<cont_gamma2_omega_face<<endl;
	cout <<"Number of elements in Omega with faces in Gamma: "<<cont_gamma_omega_face<<endl;

	#endif
*/

}



	}; //end of build_tissue_boundary_transp

	//Build boundary regions on network

	void 
	transport3d1d::build_vessel_boundary_transp(void)
	{
	#ifdef M3D1D_VERBOSE_
	cout << "Building vessel boundary ..." << endl;
	#endif
	try {

	dal::bit_vector junctions; // global idx of junctions vertices in meshv
	dal::bit_vector extrema;   // global idx of extreme vertices in meshv

	Jv_transp.clear();
	nb_extrema=0; 
	nb_junctions=0;
	
	size_type fer = nb_branches; // first empty region
	GMM_ASSERT1(meshv.has_region(fer)==0, 
		"Overload in meshv region assembling!");

	// List all the convexes
	dal::bit_vector nn = meshv.convex_index();
	bgeot::size_type cv;
	for (cv << nn; cv != bgeot::size_type(-1); cv << nn) {

		bgeot::pconvex_structure cvs = meshv.structure_of_convex(cv);
		if (cvs->nb_points()>2) 
			cerr << "Error: convex #" << cv << "has more than 2 vertices!" << endl;
		if (cvs->nb_faces()>2)  
			cerr << "Error: convex #" << cv << "has more than 2 faces!" << endl;

		// Build regions for BCs and junctions
		// Global idx of mesh vertices
		size_type i0 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(1)[0]];
		size_type i1 = meshv.ind_points_of_convex(cv)[cvs->ind_points_of_face(0)[0]];
		// Identify vertex type
		if (meshv.convex_to_point(i0).size()==1){ /* inflow extremum */
			// Update information
			extrema.add(i0);
			nb_extrema++;
			// Build a new region made by a single face
			GMM_ASSERT1(meshv.has_region(fer)==0, 
				"Overload in meshv region assembling!");
			meshv.region(fer).add(cv, 1);
			// Store the current index and then update it
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv_transp.size())) {
				found = (i0 == BCv_transp[bc].idx);
				if (!found) bc++;
			}
			GMM_ASSERT1(found=true, "Miss a boundary node in BCv list!");

			BCv_transp[bc].rg = fer; 
			fer++;
			// Store the containing branch index
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {

				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			BCv_transp[bc].branches.emplace_back(branch); 
		}
		else if (meshv.convex_to_point(i0).size()==2){ /* trivial inflow junction */
			// DO NOTHING

		}
		else if (meshv.convex_to_point(i0).size()>=2){ /* non-trivial inflow junction */
		// Check if junction has been already stored, 
			// if not add to the junction list (J) and build a new region

	//! \todo add multiple times the junction node to the junction region. Generic_assembly looks (apparently) only at indexes of the nodes, not at his coordinates; in this way, when I build the region with the junction node from a certain branch, generic_assembly will not recognize the same node from another branch (probably i look at the basis buildt only on the first branch). In order to use the generic_assembly for junction nodes I should add all the basis to the region (e.g. the same node from all the branches)

			dal::bit_vector tmp; tmp.add(i0);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i0);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 1); // single-face region
				// Create a new junction node
				Jv_transp.emplace_back("JUN", 0, i0, fer);
				fer++;
			}
			// Search for index of containing branch (\mathcal{P}^{in}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			// Add the inflow branch (to the right junction node)
			size_type jj = 0;
			bool found = false;
			while (!found && jj < nb_junctions){
				found = (i0 == Jv_transp[jj].idx);
				if (!found) jj++;
			}
			//cout << "Branch -" << branch << " added to junction " << jj << endl;
			Jv_transp[jj].value += param.R(mimv, branch);
			Jv_transp[jj].branches.emplace_back(-branch);
			GMM_ASSERT1(branch>0, 
				"Error in network labeling: -0 makes no sense");
		}
		
		if (meshv.convex_to_point(i1).size()==1){ 
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv_transp.size())) {
				found = (i1 == BCv_transp[bc].idx);
				if (!found) bc++;
			}
			if (found){ /* outlow extremum */
				extrema.add(i1); 
				nb_extrema++; 
				// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				// Store the current index and then update it
				BCv_transp[bc].value *= +1.0;
				BCv_transp[bc].rg = fer; 
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv_transp[bc].branches.emplace_back(branch); 
			}

			else { // interior -> Mixed point 
				// "MIX" label via post-processing
			// Build a new region made by a single face
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshv.region(fer).add(cv, 0);
				BCv_transp.emplace_back("MIX", 0.0, i1, fer);
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshv.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv_transp.back().branches.emplace_back(branch); 
			}

		}
		else if (meshv.convex_to_point(i1).size()==2){ /* trivial outflow junction */

			// Search for index of first containing branch (\mathcal{P}^{out}_j)
			size_type firstbranch = 0; 
			bool contained = false;
			while (!contained && firstbranch<nb_branches ) {
				contained = meshv.region(firstbranch).is_in(cv);
				if (!contained) firstbranch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i1!");

			// Check if i1 is a trivial junction (or a INT point)
			size_type cv1 = meshv.convex_to_point(i1)[0];
			size_type cv2 = meshv.convex_to_point(i1)[1];
			bool is_junc = (meshv.region(firstbranch).is_in(cv1) < 1 ||
							meshv.region(firstbranch).is_in(cv2) < 1 );
							
			if (is_junc){
				cout << "Found a trivial junction at i1 = " << i1 << endl;
				// Check if jucntion has been already stored, 
				// if not add to the junction list (J) and build a new region
				dal::bit_vector tmp; tmp.add(i1);
				if(!junctions.contains(tmp)){
					// Store the junction vertex
					junctions.add(i1);
					nb_junctions++;
					GMM_ASSERT1(meshv.has_region(fer)==0, 
						"Overload in meshv region assembling!");
					// Build a new region with idx "first empty region"
					meshv.region(fer).add(cv, 0);
					// Create a new junction node
					Jv_transp.emplace_back("JUN", 0, i1, fer);
					fer++;
				// Search for index of second containing branch (\mathcal{P}^{out}_j)
				size_type secondbranch = 0; 
				size_type secondcv = (( cv1 == cv) ? cv2 : cv1);
				size_type firstcv = (( cv1 != cv) ? cv2 : cv1);
				contained = false;
				while (!contained && secondbranch<nb_branches ) {
					if (secondbranch!=firstbranch)
					contained = meshv.region(secondbranch).is_in(secondcv);
					if (!contained) secondbranch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				// Add the two branches
				scalar_type in;
				in=0;
				if (meshv.ind_points_of_convex(firstcv)[0]==i1) in=-1;
				else if (meshv.ind_points_of_convex(firstcv)[1]==i1) in=+1;
				GMM_ASSERT1(in!=0, "There's something wrong in firstbranch convex index");
				Jv_transp.back().branches.emplace_back(in*firstbranch);

				in=0;
				if (meshv.ind_points_of_convex(secondcv)[0]==i1) in=-1;
				else if (meshv.ind_points_of_convex(secondcv)[1]==i1) in=+1;
				GMM_ASSERT1(in!=0, "There's something wrong in secondbranch convex index");
				Jv_transp.back().branches.emplace_back(in*secondbranch);
				Jv_transp.back().value += param.R(mimv, firstbranch);
				Jv_transp.back().value += param.R(mimv, secondbranch);
				}
			}
		}
		else if (meshv.convex_to_point(i1).size()>=2){ /* non-trivial outflow junction */

			// Search for index of containing branch (\mathcal{P}^{out}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshv.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");

			// Check if jucntion has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i1);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i1);
				nb_junctions++;
				GMM_ASSERT1(meshv.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshv.region(fer).add(cv, 0);
				// Create a new junction node
				Jv_transp.emplace_back("JUN", 0, i1, fer);
				// Add the outflow branch
				Jv_transp.back().branches.emplace_back(+branch);
				Jv_transp.back().value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << i1 << endl;
				fer++;
			}
			else {
				// Add the outflow branch (to the right junction node)
				size_type jj = 0;
				bool found = false;
				while (!found && jj < nb_junctions){
					found = (i1 == Jv_transp[jj].idx);
					if (!found) jj++;
				}
				Jv_transp[jj].branches.emplace_back(+branch);
				Jv_transp[jj].value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << jj << endl;
			}
		}

	} /* end of convexes loop */  

	// Ckeck network assembly
	#ifdef M3D1D_VERBOSE_
	cout << "--- NETWORK ASSEMBLY ------------------ "   << endl;
	cout << "  Branches:   " << nb_branches << endl
		 << "  Vertices:   " << nn.size()+1 << endl;
	cout << "  Extrema:    " << extrema << endl;	  
	for (size_type i=0; i<BCv_transp.size(); ++i)
		cout << "    -  label=" << BCv_transp[i].label 
			 << ", value=" << BCv_transp[i].value << ", ind=" << BCv_transp[i].idx 
			 << ", rg=" << BCv_transp[i].rg << ", branches=" << BCv_transp[i].branches << endl; 
	cout << "  Junctions: " << junctions << endl;
	for (size_type i=0; i<Jv_transp.size(); ++i)
		cout << "    -  label=" << Jv_transp[i].label 
			 << ", value=" << Jv_transp[i].value << ", ind=" << Jv_transp[i].idx 
			 << ", rg=" << Jv_transp[i].rg << ", branches=" << Jv_transp[i].branches << endl; 
	cout << "---------------------------------------- "   << endl;
	#endif

	} 
	GMM_STANDARD_CATCH_ERROR; // catches standard errors

	} /* end of build_vessel_boundary_transp */


  
	  void transport3d1d::assembly_transp (void)
	 {

 	 //1) Build the monolithic matrix AM
	 assembly_mat_transp();

	 //2) Build the monolithic rhs FM
	 assembly_rhs_transp();
	 //CASO STAZIONARIO: NON SERVE PIÙ L'UPDATE, posso assemblare il rhs con la matrice monolitica
	 
	 }; // end of assembly
 
 
 

	void  
	transport3d1d::assembly_mat_transp(void)
	{
	#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM_transp, UM_transp, FM_transp ..." << endl;
	#endif
	gmm::resize(AM_transp, dof_transp.tot(), dof_transp.tot());	gmm::clear(AM_transp);
	gmm::resize(UM_transp, dof_transp.tot()); 			gmm::clear(UM_transp);
	gmm::resize(FM_transp, dof_transp.tot()); 			gmm::clear(FM_transp);
	
	
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM_transp ..." << endl;
	#endif
	// Mass(time derivative)  matrix for the interstitial problem
	sparse_matrix_type Mt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Mt);
	// Diffusion matrix for the interstitial problem
	sparse_matrix_type Dt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Dt);
	//Advection matrix for interstitial problem
	sparse_matrix_type Bt(dof_transp.Ct(), dof_transp.Ct());gmm::clear(Bt);
	// Linfatic drainage matrix for the interstitial problem
	sparse_matrix_type Lt(dof_transp.Ct(), dof_transp.Ct()); gmm::clear(Lt);
	// Reaction matrix for the interstitial problem
	sparse_matrix_type Rt(dof_transp.Ct(), dof_transp.Ct()); gmm::clear(Rt);

	// Mass (time derivative)  matrix for the network problem
	sparse_matrix_type Mv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Mv);		
	// Diffusion matrix for the network problem
	sparse_matrix_type Dv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Dv);
	//Advection matrix for network problem
	sparse_matrix_type Bv(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Bv);
	//Oxyhemoglobin advection for network problem
	//sparse_matrix_type Av(dof_transp.Cv(), dof_transp.Cv());gmm::clear(Av);

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
	
	
	//Descriptors for the problem solved

	#ifdef M3D1D_VERBOSE_
	cout<< "   Computing Peclet number..."<<endl;	
	#endif	
	//vectors containing the exact solution
	vector_type Ut(dof.Ut()); 	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())), 		Ut);
	vector_type Uv(dof.Uv()); 	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())), Uv);
	//Interpolate Ut on polinomials of degree 0 
	pfem pf_U = fem_descriptor("FEM_PK(3,0)");
	mesh_fem mf_U(mesht);
	mf_U.set_qdim(3);
	mf_U.set_finite_element(mesht.convex_index(), pf_U);
	vector_type Ut_(mf_U.nb_dof());
	getfem::interpolation(mf_Ut, mf_U, Ut, Ut_);
	//compute peclet
	scalar_type peclet_v= peclet(meshv, Uv, param_transp.Av(1), 1);
	scalar_type peclet_t= peclet(mesht, Ut_, param_transp.At(1), 3);
	
	#ifdef M3D1D_VERBOSE_
	cout<< "Peclet in vessels:    "<< peclet_v<<endl;
	cout<< "Peclet in tissue:    "<< peclet_t<<endl;
	#endif	
	
//***************************************inizio ciclo iterativo

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mt, Dt, Lt and Rt ..." << endl;
	#endif


	//Coefficient for mass term:
	vector_type linf_coeff(dof.Pt()); gmm::clear(linf_coeff);
	vector_type Pl(dof.Pt(),PARAM.real_value("PL"));
	gmm::scale(Pl, -1.0); 
	gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut(), dof.Pt())) ,  linf_coeff);
	gmm::add(Pl ,  linf_coeff);
	gmm::scale (linf_coeff,param_transp.Q_pl(1));

	//Coefficient for reaction term:
	vector_type consump_coeff(dof_transp.Ct()); gmm::clear(consump_coeff);

	vector_type PRESS50(dof_transp.Ct(), PARAM.real_value("Pm_50")); //riempio il vettore con lo scalare Pm_50
	gmm::scale(PRESS50, param_transp.alpha_t_); 
	vector_type ct_guess(dof_transp.Ct(), param_transp.Ct0_);
	gmm::add(PRESS50, ct_guess, consump_coeff);
	//c'è un modo migliore?
	for (size_type i=0; i<dof_transp.Ct(); i++)
	{
	consump_coeff[i] = 1.0/consump_coeff[i];
	}
	gmm::scale(consump_coeff, param_transp.m0_);
	 

	//Build Mt, Dt, Lt and Rt 
	asm_tissue_transp(Mt, Dt, Lt, Rt, mimt, mf_Ct, mf_coeft,  param_transp.At(), linf_coeff, consump_coeff );

	// Copy Mt: Time Derivative in tissue
	if(descr_transp.STATIONARY ==0)
	{ 
	gmm::scale(Mt, (1.0/param_transp.dt()));
	gmm::add(Mt,  
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 
	}
	
		
	// Check peclet number for instability
	if((descr_transp.ADVECTION==1) && (peclet_t>1))
		{ cout<<"WARNING!! Peclet > 1 in tissue: applying artificial diffusion"<<std::endl;	
	  	  gmm::scale(Dt, (1+peclet_t));}
	
	// Copy Dt: diffusion in tissue		  
	gmm::add(Dt,
			 gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct())));  
	
		
	// Copy Rt: reaction in tissue
 	gmm::add(Rt, 
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
				 	gmm::sub_interval(0, dof_transp.Ct()))); 
	// Copy Lt: linfatic drainage in tissue
	 gmm::add(Lt, 
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
				 	gmm::sub_interval(0, dof_transp.Ct()))); 

	
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mv and Dv ..." << endl;
	#endif	

	

	// Build Mv and Dv
	asm_network_transp(Mv, Dv, mimv,mf_Cv, mf_coefv, param_transp.Av(), param.R());

		
	// Copy Mv: Time Derivative in network
	if(descr_transp.STATIONARY ==0)
	{
	gmm::scale(Mv, (1.0/param_transp.dt()));
	gmm::add(Mv, 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))); 
	}

	// Check peclet number for instability
	 if((descr_transp.ADVECTION==1) && (peclet_v>1))
		{ cout<<"WARNING!! Peclet > 1 in network: applying artificial diffusion"<<endl;
   	 	  gmm::scale(Dv, (1+peclet_v)); }
		
	// Copy Dv: diffusion in network 	
	gmm::add(Dv, 	gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	
	
	 
	
	
	if(descr_transp.ADVECTION ==0)	{cout<<"No advection: only diffusion and reaction terms"<<endl;}
	else{		
	//ADVECTION	
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Bt and Bv ..." << endl;
	#endif	
	
		
	//ADVECTION IN TISSUE		
	// Build Bt
	asm_advection_tissue(Bt, mimt, mf_Ct, mf_Ut, gmm::sub_vector(UM, gmm::sub_interval(0, dof.Ut())));
	// Copy Bt: advection in tissue
	gmm::add(Bt,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 	
		
	size_type shift = 0;
	// build Bv and Av
	for(size_type i=0; i<nb_branches; ++i){

		if(i>0) shift += mf_Uvi[i-1].nb_dof();
		vector_type Uvi( mf_Uvi[i].nb_dof()); gmm::clear(Uvi);
		gmm::add(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift, mf_Uvi[i].nb_dof())) ,  Uvi);


		asm_advection_network(Bv, mimv, mf_Cv, mf_coefvi[i], mf_Uvi[i], mf_coefv, Uvi, param.lambdax(i), param.lambday(i), param.lambdaz(i),  param.R(), meshv.region(i) );
		//asm_hemoadvection_network((Av, mimv, mf_Cv, mf_coefvi[i], mf_Uvi[i], mf_coefv, Uvi, param.lambdax(i), param.lambday(i), param.lambdaz(i),  param.R(), depsi, meshv.region(i) );

	}
	gmm::scale(Bv, pi);
	//gmm::scale(Av,pi);

	// Copy Av and Bv: advection in network
	gmm::add(Bv,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	
	
	/*
	gmm::add(Av,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	*/
		}/* end of advection assembly*/ 


	bool COUPLING = PARAM.int_value("COUPLING", "flag for coupling-exchange term ");
	if(COUPLING==0)  { cout<< "Uncoupled problem: no exchange between tissue and vessels"<<endl; }
	else{
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling aux exchange matrices Mbar and Mlin ..." << endl;
	#endif

	if(PARAM.int_value("couple", "flag for coupling function (notaro 0, brambilla 1)")){
   		bool READ_ITERPOLATOR = PARAM.int_value("READ_ITERPOLATOR", "flag for read interpolator from file ");
    	if (!READ_ITERPOLATOR){
		asm_exchange_aux_mat_transp(Mbar, Mlin, 
			mimv, mf_Ct, mf_Cv, mf_coefv, param.R(), descr.NInt, nb_branches);
      	std::ostringstream mbar_name,mlin_name;
      	mbar_name << descr_transp.OUTPUT+"Mbar_transp";
      	mlin_name << descr_transp.OUTPUT+"Mlin_transp";
      	gmm::MatrixMarket_IO::write(mbar_name.str().c_str(), Mbar);
      	gmm::MatrixMarket_IO::write(mlin_name.str().c_str(), Mlin);
    	}
    	else
    	{
      	std::ostringstream mbar_name,mlin_name;
      	mbar_name << descr_transp.OUTPUT+"Mbar_transp";
      	mlin_name << descr_transp.OUTPUT+"Mlin_transp";
      	gmm::MatrixMarket_load(mbar_name.str().c_str(), Mbar);
      	gmm::MatrixMarket_load(mlin_name.str().c_str(), Mlin);
    	}
	}
	if(!PARAM.int_value("couple", "flag for coupling function (notaro 0, brambilla 1)")){
   		bool READ_ITERPOLATOR = PARAM.int_value("READ_ITERPOLATOR", "flag for read interpolator from file ");
    	if (!READ_ITERPOLATOR){
			asm_exchange_aux_mat(Mbar, Mlin, 
				mimv, mf_Ct, mf_Cv, param.R(), descr.NInt);
      	std::ostringstream mbar_name,mlin_name;
      	mbar_name << "./vtk/Mbar";
      	mlin_name << "./vtk/Mlin";
      	gmm::MatrixMarket_IO::write(mbar_name.str().c_str(), Mbar);
      	gmm::MatrixMarket_IO::write(mlin_name.str().c_str(), Mlin);
    	}
    	else
    	{
      	std::ostringstream mbar_name,mlin_name;
      	mbar_name << "./vtk/Mbar";
      	mlin_name << "./vtk/Mlin";
      	gmm::MatrixMarket_load(mbar_name.str().c_str(), Mbar);
      	gmm::MatrixMarket_load(mlin_name.str().c_str(), Mlin);
    	}
	}
	
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling exchange matrices ..." << endl;
	#endif


	bool NEWFORM = PARAM.int_value("NEW_FORMULATION", "flag for the new formulation");
	

	// bluid oncotic term fot tissue
	vector_type ONCOTIC (dof.Pv());
	gmm::copy(gmm::sub_vector(UM, 
		  		  gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())),
		  ONCOTIC);
	gmm::mult_add(gmm::scaled(Mbar,-1.0), 
		  gmm::sub_vector(UM, 
		  		  gmm::sub_interval(dof.Ut(), dof.Pt())),
		  ONCOTIC);


	scalar_type picoef=param.sigma()*(param.pi_v()-param.pi_t());
        vector_type DeltaPi(dof.Pv(),picoef);
        gmm::add(gmm::scaled(DeltaPi,-1.0), ONCOTIC);	
	gmm::scale(ONCOTIC,(1.0-param.sigma())*param.Q(0));
	
	// build permeability term for tissue
	vector_type PERM (dof.coefv());
	gmm::copy(param.R(), PERM);
	gmm::scale(PERM, 2*pi*param_transp.Y()[0]);




	//build exchange matrixes for tissue
	asm_exchange_mat_transp(Btt, Btv, Bvt, Bvv,
			mimv, mf_Cv, mf_coefv, mf_Pv, Mbar, Mlin, 
			ONCOTIC, PERM, NEWFORM);


	// Copying Btt
	gmm::add(Btt,			 
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()), 
					gmm::sub_interval(0, dof_transp.Ct()))); 
	// Copying Btv

	gmm::add(Btv,
			gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(0, dof_transp.Ct()),
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))); 
	
	// Copying Bvt

	gmm::add(Bvt,								
			  gmm::sub_matrix(AM_transp, 
			  		gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()),
					gmm::sub_interval(0, dof_transp.Ct())));
	// Copying Bvv

	gmm::add(Bvv,
			  gmm::sub_matrix(AM_transp, 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()), 
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))); 
	

	}	

	#ifdef M3D1D_VERBOSE_
	cout << "  Setting initial condition for tissue and network concentration ..." << endl;
	#endif

	vector_type C0t_vect(dof_transp.Ct(), param_transp.C0t());
	vector_type C0v_vect(dof_transp.Cv(), param_transp.C0v());

	gmm::copy(C0t_vect,
		  gmm::sub_vector(UM_transp, 
			  	gmm::sub_interval(0, dof_transp.Ct()))	);
	gmm::copy(C0v_vect,
		  gmm::sub_vector(UM_transp, 
			  	gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))	);

	gmm::clear(C0t_vect);	gmm::clear(C0v_vect);

	// De-allocate memory
	gmm::clear(Mt);    gmm::clear(Mv); 
	gmm::clear(Dt);    gmm::clear(Dv);
	gmm::clear(Bt);    gmm::clear(Bv);
	gmm::clear(Mbar);  gmm::clear(Mlin);
	gmm::clear(Btt);   gmm::clear(Btv);
	gmm::clear(Bvt);   gmm::clear(Bvv);

		
	}/* end of assembly_mat_transp */


	vector_type cv_on_branch(vector_type Cv_old, size_type nb){
	vector_type cv_i(nb);

	size_type pos=0;
		for (getfem::mr_visitor mrv(mf_Cv.linked_mesh().region(i)); !mrv.finished(); ++mrv)
		for (auto b : mf_Cv.ind_basic_dof_of_element(mrv.cv()))
			{
			cv_i[pos] = Cv_old[b]; //ottengo cv_i --> la concentrazione definita ramo per ramo
			pos++;}
		return cv_i;
	}

	void 
	transport3d1d::assembly_rhs_transp(void)
	{
 
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic rhs FM_transp ... " << endl;
	#endif
	
	#ifdef M3D1D_VERBOSE_
	cout << "  Building coupling dirichlet boundary term ..." << endl;
	#endif
	//! Re-write the cycle over the boundary conditions: now, the Dirichlet conditions are a.s. correct, but in case there is a Dirichlet point of the network in a Robin face of the tissue (or viceversa) te Dirichlet condition is not completely guaranteed. One easy solution: implement all the robin condition, then implement all the Dirichlet condition.
	asm_coupled_bc_transp (AM_temp, FM_temp, mf_Ct, mf_Cv, BCt_transp, BCv_transp);
	
	#ifdef M3D1D_VERBOSE_
	cout << "  Building tissue boundary term ..." << endl;
	#endif
	
	//Right Hand Side for tissue	
	vector_type Ft(dof_transp.Ct());
			
	sparse_matrix_type Att(dof_transp.Ct(), dof_transp.Ct());

	gmm::add(	gmm::sub_matrix(AM_temp,
			gmm::sub_interval(0,dof_transp.Ct()),
			gmm::sub_interval(0,dof_transp.Ct()))
			, Att);
	gmm::scale(	gmm::sub_matrix(AM_temp,
			gmm::sub_interval(0,dof_transp.Ct()),
			gmm::sub_interval(0,dof_transp.Ct()))
			, 0.0);	
			
	gmm::add(	 gmm::sub_vector(FM_temp, gmm::sub_interval(0,dof_transp.Ct()))
			,Ft);	 
	gmm::scale(	 gmm::sub_vector(FM_temp, gmm::sub_interval(0,dof_transp.Ct()))
			,0.0);
	
	scalar_type beta_t  = PARAM.real_value("BETAtissue_transp", "Coefficient for mixed BC for transport problem in tissue");
	asm_tissue_bc_transp(Ft, Att, mimt, mf_Ct, mf_coeft, BCt_transp,beta_t);
	gmm::add(Att, 
			gmm::sub_matrix(AM_temp,
					gmm::sub_interval(0,dof_transp.Ct()),
					gmm::sub_interval(0,dof_transp.Ct())));
	gmm::add(Ft, 
			gmm::sub_vector(FM_temp,
					gmm::sub_interval(0,dof_transp.Ct())));
	// De-allocate memory
	gmm::clear(Att);
	gmm::clear(Ft);


	#ifdef M3D1D_VERBOSE_
	cout << "  Building vessel boundary term ..." << endl;
	#endif
		
	//Right Hand Side for vessels
//RR: non lavoro più con AM_temp e FM_temp ma direttamente con il vettore monolitico FM e la matrice monolitica AM --> non ho più bisogno di aggiornare il rhs nel tempo
	vector_type cv_guess (dof_transp.Cv(), param_transp.Cv0_);
	
	sparse_matrix_type Avv(dof_transp.Cv(), dof_transp.Cv());
	vector_type Fv(dof_transp.Cv());
	vector_type Ov(dof_transp.Cv());
	
	gmm::add(	gmm::sub_matrix(AM_transp,
			gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
			gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()))
			, Avv);
	gmm::scale(	gmm::sub_matrix(AM_transp,
			gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
			gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()))
			, 0.0);	
				
	gmm::add(	 gmm::sub_vector(FM_transp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))
			,Fv);	
	gmm::scale(	 gmm::sub_vector(FM_transp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv()))
			,0.0);

	scalar_type beta_v  = PARAM.real_value("BETAvessel_transp", "Coefficient for mixed BC for transport problem in vessels");
	asm_network_bc_transp(Fv, Avv, mimv, mf_Cv, mf_coefv, BCv_transp, beta_v, param.R());
	//RR: c'è bisogno di utilizzare Fv_temp e Avv_temp? li copierei direttamente nella corretta posizione
	gmm::add(Avv, 
			gmm::sub_matrix(AM_transp,
					gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv()),
					gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv())));
	gmm::add(Fv, 
			gmm::sub_vector(FM_transp,
					gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
					
					
		size_type shift =0;
		size_type shift_h=0;
		scalar_type shift_coef=0;
		
		scalar_type k2 = pow((param_transp.Ps_50_*param_transp.alpha_t_),param_transp.delta_);
		scalar_type k1 = param_transp.N_*param_transp.MCHC_;
		
	for(size_type i=0; i<nb_branches; ++i){

		vector_type Hi(mf_Hi[i].nb_dof());
		scalar_type Ri = param.R(mimv, i);
		
		vector_type psi(mf_Hi[i].nb_dof()); gmm::clear(psi);
		
		if(i>0) shift_h += mf_Hi[i-1].nb_dof();
		if(i>0) shift_coef += mf_coefvi[i-1].nb_dof();
		if(i>0) shift += mf_Uvi[i-1].nb_dof();
		
		nbdof_mfH = mf_Hi[i].nb_dof();
		cv_i = cv_on_branch(vector_type cv_guess, size_type nbdof_mfH);
		
		gmm::copy(gmm::sub_vector(H_old, 
			gmm::sub_interval(shift_h, mf_Hi[i].nb_dof())), Hi);
			
		for (size_type i=0; i<mf_Hi[i].nb_dof(); i++)
		{
		psi[i] = Hi[i]*k1*pow(cv_i[i], param_transp.delta_)/(pow(cv_i[i], param_transp.delta_)+k2);
		}
	
	asm_hemoadvection_rhs_network(Ov, mimv, mf_Cv, mf_coefvi[i], mf_Uvi[i], mf_Hi[i], Uvi, param.lambdax(i), param.lambday(i), param.lambdaz(i),  param.R(), psi, meshv.region(i));	
	}
	
	gmm::copy(Ov, gmm::sub_vector(FM_transp, 
								gmm::sub_interval(dof_transp.Ct(),dof_transp.Cv())));
	
	
	
	// De-allocate memory
	gmm::clear(Avv);
	gmm::clear(Fv);
	gmm::clear(Ov);
	}/* end of assembly_rhs_transp */
	
/*
 //IN CASO STAZIONARIO: NON SERVE PIÙ L'UPDATE DEL RHS
	void transport3d1d::update_transp (void){

	
	/*At each time step, the right hand side is modified by the time derivative term.
	Since we must ensure in a strong way the Dirichlet conditions, by modifying the monolithic matrix and the rhs vector, we save both AM_transp and FM_transp, where are assembled the stationary terms; 	then, we work on AM_temp and FM_temp, modifying them when necessary.
	

	#ifdef M3D1D_VERBOSE_
	cout << "  Update monolithic matrix and rhs vector ..." << endl;
	#endif

	gmm::copy(AM_transp, AM_temp);
	gmm::copy(FM_transp, FM_temp);


	// update rhs (time step mass term)
	vector_type TFt(dof_transp.Ct());
	vector_type TFv(dof_transp.Cv());
	asm_source_term(TFt,mimt, mf_Ct, mf_Ct,gmm::sub_vector(UM_transp, gmm::sub_interval(0, dof_transp.Ct()))); 
	vector_type RR(dof.coefv()); gmm::clear(RR);
	gmm::add(param.R(),RR); gmm::vscale(param.R(),RR); gmm:: scale (RR, pi);
	asm_source_term(TFv,mimv, mf_Cv, mf_coefv, RR);
	gmm::vscale(gmm::sub_vector(UM_transp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), TFv);
	gmm::scale(TFt, (1.0/param_transp.dt())); // dt time step
	gmm::scale(TFv, (1.0/param_transp.dt())); // dt time step
	gmm::add(TFt, gmm::sub_vector(FM_temp, gmm::sub_interval(0, dof_transp.Ct())));
	gmm::add(TFv, gmm::sub_vector(FM_temp, gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())));
	gmm::clear(UM_transp);
	gmm::clear(TFt); gmm::clear(TFv);

	//update rhs (bundary condition terms)
	assembly_rhs_transp();


	}  end of update_transp
*/

	// Aux function for solver:
	// contains the list of different methods for solving (SuperLU, SAMG,GMRES, etc)
	// Always passes through the temporary matrixes AM_temp and FM_temp 
	bool transport3d1d::solver (const size_type dof1, 
					   const size_type dof2,
					   const size_type dof3,
					   const size_type dof4){

	gmm::csc_matrix<scalar_type> A_transp;
	gmm::clean(AM_transp, 1E-12);
	gmm::copy(AM_temp, A_transp);
	
	vector_type F_transp(gmm::vect_size(FM_transp));
	gmm::clean(FM_transp, 1E-12);
	gmm::copy(FM_temp, F_transp);
	
		
	if ( descr_transp.SOLVE_METHOD == "SuperLU" || descr_transp.SOLVE_METHOD == "SUPERLU" ) { // direct solver //
		#ifdef M3D1D_VERBOSE_
		cout << "  Applying the SuperLU method ... " << endl;
		#endif
		scalar_type cond;
		gmm::SuperLU_solve(A_transp, UM_transp, F_transp, cond);
		cout << "  Condition number (transport problem): " << cond << endl;
	}
	else if(descr_transp.SOLVE_METHOD == "SAMG"){
	#ifdef WITH_SAMG	
	#ifdef M3D1D_VERBOSE_
		cout << "Solving the monolithic system ... " << endl;
	#endif



	//////////////////////////////////////AMG INTERFACE
	#ifdef M3D1D_VERBOSE_
	std::cout<<"converting A"<<std::endl;
	#endif
	gmm::csr_matrix<scalar_type> A_csr;
	gmm::clean(AM_temp, 1E-12);


	int dim_matrix=dof1+dof2+dof3+dof4;
	gmm::copy(gmm::sub_matrix(AM_temp,
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
	gmm::copy(gmm::sub_vector(FM_temp,gmm::sub_interval(0,dim_matrix)),B);



	AMG amg("3d1d");
	amg.set_dof(dof1,dof2,dof3,dof4);

	
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

	return true;
	} /* end of solver_transp */



	 bool transport3d1d::solve_transp (void)
 	{
  	#ifdef M3D1D_VERBOSE_
	cout << "Solving the monolithic system ... " << endl;
	#endif
	gmm::resize(AM_temp, dof_transp.tot(), dof_transp.tot()); gmm::clear(AM_temp);
	gmm::resize(FM_temp, dof_transp.tot()); gmm::clear(FM_temp);
	
	double time = gmm::uclock_sec();
	double time_partial = 0;
	double time_count = 0;	



	for(double t=0;t<=(param_transp.T()+ param_transp.dt())*(!descr_transp.STATIONARY) ; t = t + param_transp.dt() + 0.1*(param_transp.dt()==0) ){ 
	time_count++; 
	time_partial=gmm::uclock_sec();

	if(descr_transp.STATIONARY){
	std::cout<<"Stationary problem... "<<std::endl;	
	}
	else{
	std::cout<<"-------------------------------------------"<<std::endl;
	std::cout<<"Iteration number: "<<time_count<<std::endl;
	std::cout<<"time = "<<t<<" s"<<std::endl<<std::endl;	
	}

	//Update rhs and boundary condition
	//update_transp();
	
	// Solve the system on AM_temp, UM_transp, FM_temp(???)
	bool solved = solver(dof_transp.Ct(),dof_transp.Cv(),0,0);
	if (!solved) return false;
	
//***** decido di portarlo dentro a fixpoint_concentration
	//export solution
	#ifdef M3D1D_VERBOSE_
	std::cout<<"solved! going to export..."<<std::endl;
	#endif	
	
	string time_suff = "";
	std::ostringstream convert;
	convert << time_count;
	time_suff = convert.str();
	export_vtk_transp(time_suff); 



	#ifdef M3D1D_VERBOSE_		
	std::cout<<"exported!"<<std::endl;
	#endif	
	
	if(!descr_transp.STATIONARY)
	cout << "... time to solve : "	<< gmm::uclock_sec() - time_partial << " seconds\n";
	
	} //end of cycle over time 
	if(!descr_transp.STATIONARY){
	cout << endl<<"... time to solve all the time steps: " << gmm::uclock_sec() - time << " seconds\n";				}
	else{
	cout << endl<<"... time to solve : " << gmm::uclock_sec() - time << " seconds\n";
	}
//********************************************************
	return true;
 	}; // end of solve_transp
	

//***************************************fine ciclo iterativo

	//Compute the residuals for mass balance at each junction 
	void transport3d1d::mass_balance(void){

	#ifdef M3D1D_VERBOSE_		
	cout << " Compute MBD and MBA "   << endl;
	#endif	

	// initialize the MBD and MBA to zero (clear at eac time step)
	for (size_type i=0; i<Jv_transp.size(); ++i){
		Jv_transp[i].MBD=0;
		Jv_transp[i].MBA=0;
	}	

	size_type shift = 0; //counter for branches
	
	for (size_type i=0; i<mf_Uvi.size(); ++i){ // branch loop 	
		if(i>0) shift += mf_Uvi[i-1].nb_dof(); 
		mesh_region &rg_branch = meshv.region(i); // branch region

		for (size_type j=0; j<Jv_transp.size(); ++j){ //junction loop
			mesh_region &rg_junction = meshv.region(Jv_transp[j].rg);  // junction region
			// Iterators for all the branches which flow in the junction j
			std::vector<long signed int>::const_iterator bb = Jv_transp[j].branches.begin();
			std::vector<long signed int>::const_iterator be = Jv_transp[j].branches.end();
		
			//Check if outflow of branch i is in junction j			
			if ((std::find(bb, be, +i) != be)){

				// find last element of the branch
				getfem::mr_visitor ii(rg_branch); int temp=0;
				for(; !ii.finished() && temp==0; ++ii){ // loop for convexes of the branch i
					getfem::mr_visitor ii_temp=ii;
					++ii_temp;
					if(ii_temp.finished()) temp=1;
				} 
				
				//import radius of the branch
				scalar_type Ri = param.R(mimv, i);
				// find the dof for last point for Cv e Uv
				size_type last_C, last_C2, last_U;
				vector_type dof_enum_C,dof_enum_U;
				int fine_C=0, fine_U=0;
				for (mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
				for (auto b : mf_Cv.ind_basic_dof_of_element(mrv.cv()))
					{dof_enum_C.emplace_back(b);
					fine_C++;}			
				for (auto b : mf_Uvi[i].ind_basic_dof_of_element(mrv.cv()))
					{dof_enum_U.emplace_back(b);
					fine_U++;}			
				}
	

				last_C = dof_enum_C[fine_C-1];
				last_C2= dof_enum_C[fine_C-2];
				last_U = dof_enum_U[fine_U-1];
				dof_enum_C.clear(); dof_enum_U.clear();
	
				scalar_type DIFF=0, ADV=0;

				//Compute the diffusive flux
				DIFF = pi* Ri * Ri*(
						UM_transp[dof_transp.Ct()+last_C]-UM_transp[dof_transp.Ct()+last_C2] )
						/estimate_h(meshv, ii.cv()) ;
				//Compute the advective fluxes
				ADV = pi* Ri * Ri*UM[dof.Ut()+dof.Pt()+shift+last_U]*UM_transp[dof_transp.Ct()+last_C];
	
				#ifdef M3D1D_VERBOSE_		
				cout << "------------------------------------------ "   << endl;
				cout <<"in branch "<< i << " and junction "<< j <<" (region number :  "<< Jv_transp[j].rg<<" )"<<endl;
				cout<<"MBD_partial = "<< DIFF<< endl;
				cout<<"MBA_partial = "<< ADV<<  endl;
				#endif	
		
				Jv_transp[j].MBD -= DIFF;
				Jv_transp[j].MBA -= ADV;
					
			}// end of check for outflow branch
	
			//Check if inflow of branch i is in junction j			
			if ( (i!=0) &&  (std::find(bb, be, -i) != be )){  //notice that, in build_vessel_transp, we assume that the branch zero cannot have the inflow in a junction. 


				// find first element of the branch
				getfem::mr_visitor ii(rg_branch);
			
				//import radius of the branch
				scalar_type Ri = param.R(mimv, i);
				// find the dof for last point for Cv e Uv
				size_type first_C, first_C2, first_U;
				vector_type dof_enum_C,dof_enum_U;
				int fine_C=0, fine_U=0;
				for (mr_visitor mrv(rg_branch); !mrv.finished(); ++mrv){
				for (auto b : mf_Cv.ind_basic_dof_of_element(mrv.cv()))
					{dof_enum_C.emplace_back(b);
					fine_C++;}			
				for (auto b : mf_Uvi[i].ind_basic_dof_of_element(mrv.cv()))
					{dof_enum_U.emplace_back(b);
					fine_U++;}			
				}
	
				first_C = dof_enum_C[0];
				first_C2= dof_enum_C[1];
				first_U = dof_enum_U[0];
				dof_enum_C.clear(); dof_enum_U.clear();
	
				scalar_type DIFF=0,ADV=0;
	
				//Compute the diffusive flux
				 DIFF = pi* Ri * Ri*(
						UM_transp[dof_transp.Ct()+first_C2]-UM_transp[dof_transp.Ct()+first_C] )
						/estimate_h(meshv, ii.cv()) ;
	
	
				//Compute the advective fluxes
				ADV = pi* Ri * Ri*UM[dof.Ut()+dof.Pt()+shift+first_U]*UM_transp[dof_transp.Ct()+first_C];
	
				#ifdef M3D1D_VERBOSE_		
				cout << "------------------------------------------ "   << endl;
				cout <<"in branch "<< i << " and junction "<< j <<" (region number :  "<< Jv_transp[j].rg<<" )"<<endl;
				cout<<"MBD_partial = "<< DIFF<< endl;
				cout<<"MBA_partial = "<< ADV<<  endl;
				#endif	
		
				Jv_transp[j].MBD += DIFF;
				Jv_transp[j].MBA += ADV;
					
			}// end of check for outflow branch

		} //end of junction loop
	} // end of branch loop	





	cout << "  Junctions: " << endl;
	for (size_type i=0; i<Jv_transp.size(); ++i){
		cout << "    -  label=" << Jv_transp[i].label 
			 << ", value=" << Jv_transp[i].value << ", ind=" << Jv_transp[i].idx 
			 << ", rg=" << Jv_transp[i].rg << ", branches=" << Jv_transp[i].branches << endl; 
		cout << " Mass balance of diffusive fluxes = " << Jv_transp[i].MBD << endl; 
		cout << " Mass balance of advective fluxes = " << Jv_transp[i].MBA << endl;
		cout << "             ------------------- "   << endl;
	} 	
		cout << "----------------------------------------------- "   << endl;


	
	}; // end of mass_balance


 

 const void transport3d1d::export_vtk_transp (const string & time_suff,const string & suff)
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
	vector_type Cv(dof_transp.Cv()); 

	//Copy solution
	gmm::copy(gmm::sub_vector(UM_transp, 
		gmm::sub_interval(0, dof_transp.Ct())),  Ct);
	gmm::copy(gmm::sub_vector(UM_transp, 
		gmm::sub_interval(dof_transp.Ct(), dof_transp.Cv())), Cv);

	

	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Ct ..." << endl;
	#endif
	vtk_export exp_Ct(descr_transp.OUTPUT+"Ct"+suff+"_t"+time_suff+".vtk");
	exp_Ct.exporting(mf_Ct);
	exp_Ct.write_mesh();
	exp_Ct.write_point_data(mf_Ct, Ct, "Ct");



	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Cv ..." << endl;
	#endif
	vtk_export exp_Cv(descr_transp.OUTPUT+"Cv"+suff+"_t"+time_suff+".vtk");
	exp_Cv.exporting(mf_Cv);
	exp_Cv.write_mesh();
	exp_Cv.write_point_data(mf_Cv, Cv, "Cv");

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview" << endl; 
	#endif

//	If you don't want to use vtk:
/*	use this method to export a vector:

	std::ofstream outFF(descr_transp.OUTPUT+"FF.txt");
		outFF << gmm::col_vector(F_transp);
		outFF.close(); 
*/ 

/*	use this method to export a matrix 

	gmm::MatrixMarket_IO::write("SM.mm" , SM);
*/

  }
 }; // end of export_transp
 
  
  

  // Interface with problem3d1d class
  	//! Initialize the problem
	void transport3d1d::init_fluid (int argc, char *argv[])
	{ 
	problem3d1d::init(argc, argv);
	};
	
	//! Assemble the problem
	void transport3d1d::assembly_fluid (void)
	{ 
	problem3d1d::assembly();
	};
	//! Solve the problem
	bool transport3d1d::solve_fluid (void)
	{ 
	return problem3d1d::solve();
	};	
	//! Export the solution
	const void transport3d1d::export_vtk_fluid (const string & suff)
	{ 
	problem3d1d::export_vtk(suff);
	};
	


 
 } // end of namespace