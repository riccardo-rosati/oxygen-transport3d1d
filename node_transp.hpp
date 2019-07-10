
 /* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Transport Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2015-2016
                  
                Copyright (C) 2016 Stefano Brambilla
======================================================================*/
/*! 
  @file   node.hpp
  @author Stefano Brambilla <s.brambilla93@gmail.com>
  @date   September 2016 - May 2018
  @brief  Definition of the class node_transp.
 */

 
#ifndef M3D1D_NODE_TRANSP_HPP_
#define M3D1D_NODE_TRANSP_HPP_

namespace getfem{

//! Class to handle the boundary and junction nodes
struct node_transp {

	//! Label ('INT','DIR','MIX','JUN')
	std::string label; 
	//! Boundary value (useless for 'JUN')
	scalar_type value;
	//! Global index
	size_type   idx;
	//! Associated mesh region
	size_type   rg;
	//! Possible list of intersecting vessel branches
	std::vector<long signed int> branches;
	//! Balance mass from diffusive fluxes
	scalar_type MBD;
	//! Balance mass from advective fluxes
	scalar_type MBA;

	//! Constructor
	node_transp(const std::string & label_="", 
		 const scalar_type & value_=0, 
		 const size_type & idx_=0, 
		 const size_type & rg_=0,
		 const scalar_type & MBD_=0,
		 const scalar_type & MBA_=0) 
		 : label(label_), value(value_), idx(idx_), rg(rg_), MBD(MBD_), MBA(MBA_)
	{}
	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const node_transp & N
		)
	{ 
		out << "('" << N.label << "'," 
			<< N.value << "," 
			<< N.idx   << ","
			<< N.rg    << ")";

		return out;            
	}

};

}

#endif
