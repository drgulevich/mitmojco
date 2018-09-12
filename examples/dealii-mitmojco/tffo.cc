/* 
	Model of T-junction flux-flow oscillator,
    see Supplementary Material to the Ref. 
    D. R. Gulevich, V. P. Koshelets, and F. V. Kusmartsev, arXiv:1709.04052. 
*/
//==============================================//
//---------- by Dmitry R. Gulevich -------------//
//--------- drgulevich@corp.ifmo.ru ------------//
//--- ITMO University, St Petersburg, Russia ---//
//==============================================//
/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2006 - 2017 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2006
 */
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/table_handler.h> // output to table
#include <deal.II/base/timer.h> // timer
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_ilu.h> // Jacobi preconditioner
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h> // mesh output
#include <deal.II/grid/grid_in.h> // load grid from file
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>
#include <iostream>
#include <numeric> // std:accumulate
#include <iomanip>      // std::setprecision
#include <sys/time.h>
#include <sys/resource.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>

//#define timer_on true

#define _USE_MATH_DEFINES // Math constants for C++
#include <cmath>

#include <cstdlib> // atof
#include <math.h>
#include <stdlib.h> // malloc
//#include <fftw3.h>
#include <string> // string

//#include <utility>  //declarations of unique_ptr

#include <omp.h> // OpenMP
#include <mitmojco/mitmojco.h> // MiTMoJCo header file
#include <mitmojco/opt_filter.h> // optimum filtration

#define OMP_NUM_THREADS 1 //  number of OpenMP threads (use 1 for small Josephson contact)
#define debug_on false


namespace FluxFlow
{
  using namespace dealii;


  template <int dim>
  class InitialValues : public Function<dim>
  {
  public:
    InitialValues (const std::string geometry_type_in, const double h_in, const double gamma_in, const double L0_in) : Function<dim>()
		{
			h=h_in;
			gamma=gamma_in;
			L0=L0_in;
			geometry_type=geometry_type_in;
		}
    double value (const Point<dim> &p,
                          const unsigned int component = 0) const;
  private:
	double h;
	double gamma;
	double L0;
	std::string geometry_type;
  };


  template <int dim>
  double InitialValues<dim>::value (const Point<dim> &p,
                                    const unsigned int /*component*/ ) const
  {
 if(geometry_type=="ffo-sharp")
	return -h*p[0];
 else if(geometry_type=="tffo")
	return -4.*h*std::pow(p[0],3)/(3.*L0*L0);
 else {
	std::cout << "# Initial values for geometry type " << geometry_type << " are not set. Starting from 0." << std::endl;
	return 0.;
	}
  } 



  template <int dim>
  class FFO
  {
  public:
    FFO (double mfield, double time_step);
    ~FFO();
    void display_info () const;
    void set_bc (const double h_field, const double gamma, const double htffo, const double hx);
    void set_ic (const double h_field, const double gamma);
    int advance_by (unsigned int n_steps);
    int advance_by (unsigned int n_steps, char output_regime);
    void shift_phase();
    void set_voltage_filters();
    void copy_state_from(const FFO<dim>& object);
    double distance_from(FFO<dim>& object);

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    DoFHandler<dim>      dof_handler;
    const unsigned int gauss_degree;
    std::set<types::boundary_id> filter_boundary_indicators;
    std::set<types::boundary_id> rad_boundary_indicators;
    std::set<types::boundary_id> boundary7; // AJTL of TFFO
    std::vector<bool> filter_boundary_dofs;
    std::vector<bool> rad_boundary_dofs;
    std::vector<bool> boundary7_dofs; // AJTL of TFFO
    unsigned int n_filter_boundary_dofs;
    unsigned int n_rad_boundary_dofs;
    unsigned int n_boundary7_dofs;

    double V_dc_f1;
    double V_dc_f3;
    double V_dc_f5;
	const std::string geometry_type;
    const std::string mesh_file;
	const std::string amp_file;
	const double a_supp;
	const double kgap;
    double lambdaJ;
    double L0, W0, S0, P0;
    double L1, W1, S1, P1;
    double Area;
    double alpha0;
    const double beta0;
    const double h0_field;
    const double time_step;
    Vector<double>  U_current, U_old;
	TunnelCurrentType *ffo_tunnel_current;
	OptFilterType *voltage_filter_1;
	OptFilterType *voltage_filter_3;
	OptFilterType *voltage_filter_5;

  private:
    void make_grid ();
    void set_bc_ffo_rect (const double h_field, const double gamma);
    void set_bc_ffo_sharp (const double h_field, const double gamma);
    void set_bc_tffo_sharp (const double h_field, const double gamma, const double htffo, const double hx);
    void setup_system ();
    void assemble_matrices ();
    void compute_sin_vector (const Vector<double> &argument,
                          Vector<double> &sin_vector) const;
    unsigned int solve ();
    void vtk_output (const double t) const;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    SparseMatrix<double> system_matrix;
    SparseMatrix<double> p_matrix;
    SparseMatrix<double> m_matrix;
    SparseMatrix<double> alpha_mass_matrix;
    SparseMatrix<double> new_matrix;
    SparseMatrix<double> c_matrix;
    SparseMatrix<double> old_matrix;

    Vector<double>       U_next;
    Vector<double>       system_rhs;
    Vector<double> 		 boundary_forces;

    const unsigned int vtk_skip;
    const double vtk_tmin;
    const double vtk_tmax;



  };


  // Constructor
    template <int dim>
    FFO<dim>::FFO (double mfield, double time_step_in)
      :
    fe (1),
    dof_handler (triangulation),
    gauss_degree(2),
	geometry_type("tffo"),
    mesh_file("tffo.msh"),
	amp_file("../../amplitudes/NbNbN_4K2_008.fit"), // tunnel current amplitudes file
	a_supp(0.7), // pair current suppression
	kgap(4.0), // normalized gap frequency (omega_g/omega_J)
    beta0(0.017), // surface damping
    h0_field(mfield),
    time_step (time_step_in),
    vtk_skip (10), // time steps to skip
    vtk_tmin(0.),
    vtk_tmax(0.) // no vtk if < tmin
    {
    	lambdaJ = 5.5; // um

    	L0 = 400./lambdaJ;
    	W0 = 16./lambdaJ;
    	S0 = 40./lambdaJ;
    	P0 = 2./lambdaJ;

    	L1 = 228./lambdaJ;
    	W1 = 4./lambdaJ;
    	S1 = 20./lambdaJ;
    	P1 = 2./lambdaJ;

    	if(geometry_type == "tffo")
    		Area = L0*W0 - S0*(W0-P0) + (L1-W0)*W1 - 0.5*S1*(W1-P1);
    	else
			std::cout << "unknown geometry type" << std::endl;

		V_dc_f1 = 0.;
		V_dc_f3 = 0.;
		V_dc_f5 = 0.;

    	make_grid ();
    	setup_system ();

    	/* Create tunnel current object (TunnelCurrentType pointer) by calling the constructor with arguments:
    			1. Tunnel current amplitudes file.
    			2. Pair current suppression parameter (1 if no suppression).
    			3. Normalized gap frequency (omega_g/omega_J).
    			4. Integration time step.
    			5. Total number of nodes (size of array phi).
    			6. Pointer to array phi.
    			7. Number of shadow nodes to be skipped.
    			8. Pointer to the array of indices of the shadow nodes.
    		*/
    	ffo_tunnel_current = mitmojco_create( amp_file.c_str(), a_supp, kgap,
    					time_step, dof_handler.n_dofs(), U_current.begin(), 0, NULL);

    	/* Throw exception if errors occurred */
    	if( ffo_tunnel_current->error )
    		throw "MiTMoJCo error";

  		alpha0 = ffo_tunnel_current->alphaN;

   		/* Create filter object (OptFilterType pointer) for voltage filtration.
    			The constructor accepts one argument - the optimum filtration level
    			(parameter n in Ref. A. A. Odintsov, V. K. Semenov and A. B. Zorin,
    			IEEE Trans. Magn. 23, 763 (1987)). n=1 corresponds to the arithmetic mean.
   		 */
   		voltage_filter_1 = opt_filter_create(1);
   		voltage_filter_3 = opt_filter_create(3);
   		voltage_filter_5 = opt_filter_create(5);

   		// Extract filter boundary dofs
   		filter_boundary_dofs.assign(dof_handler.n_dofs(),false);
   		DoFTools::extract_boundary_dofs (dof_handler, ComponentMask(),
      			filter_boundary_dofs, filter_boundary_indicators);
   		n_filter_boundary_dofs=0;
   		for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   			if (filter_boundary_dofs[i] == true)
   				n_filter_boundary_dofs++;
   		std::cout<<"# " << n_filter_boundary_dofs << " filter boundary dofs:"<<std::endl;
   		for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   			if (filter_boundary_dofs[i] == true)
   				std::cout << "# 	i=" << i <<std::endl;
   		std::cout<<"# " << std::endl;

   		// Extract radiation boundary dofs
   		rad_boundary_dofs.assign(dof_handler.n_dofs(),false);
   		DoFTools::extract_boundary_dofs (dof_handler, ComponentMask(),
   	  			rad_boundary_dofs, rad_boundary_indicators);
   		n_rad_boundary_dofs=0;
   		for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   			if (rad_boundary_dofs[i] == true)
   				n_rad_boundary_dofs++;
   		std::cout<<"# " << n_rad_boundary_dofs << " radiation boundary dofs:"<<std::endl;
   		for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   			if (rad_boundary_dofs[i] == true)
   				std::cout << "# 	i=" << i <<std::endl;
   		std::cout<<"# " << std::endl;

   		// Extract boundary 7 dofs for TFFO
   		if(geometry_type == "tffo") {
   			boundary7_dofs.assign(dof_handler.n_dofs(),false);
   			DoFTools::extract_boundary_dofs (dof_handler, ComponentMask(),
   				boundary7_dofs, boundary7);
   			n_boundary7_dofs=0;
   			for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   				if (boundary7_dofs[i] == true)
   					n_boundary7_dofs++;
   			std::cout<<"# " << n_boundary7_dofs << " AJTL end dofs:"<<std::endl;
   			for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   				if (boundary7_dofs[i] == true)
   					std::cout << "# 	i=" << i <<std::endl;
   			std::cout<<"# " << std::endl;
   		}

   		assemble_matrices ();

   		display_info();

    }

    // Destructor
    template <int dim>
    FFO<dim>::~FFO ()
    {
    /* Clear memory allocated for the objects */
    mitmojco_free( ffo_tunnel_current );
    opt_filter_free(voltage_filter_1);
    opt_filter_free(voltage_filter_3);
    opt_filter_free(voltage_filter_5);
    }


  template <int dim>
  void FFO<dim>::make_grid () {

  GridIn<dim> gridin;
  gridin.attach_triangulation(triangulation);
  std::string filename =  mesh_file;
  std::ifstream f(filename.c_str());
  gridin.read_msh(f);
  std::cout << "# Mesh loaded from " << filename << std::endl;

  // Set filter boundary indicators. Injection (left) boundary is used here
  // for average voltage readings because of a more regular dynamics
  filter_boundary_indicators.insert(1);
  std::cout << "# filter boundary indicators: ";
  std::set<types::boundary_id>::iterator it = filter_boundary_indicators.begin();
  for(unsigned int i=0;i<filter_boundary_indicators.size();i++) {
	  std::cout << (unsigned int)*it << " ";
	  std::advance(it, 1);
  }
  std::cout << std::endl;

  // Set radiation boundary indicators.
  // Radiation boundary is the one coupled to the load (e.g. SIS mixer).
  rad_boundary_indicators.insert(2);
  std::cout << "# radiation boundary indicators: ";
  it = rad_boundary_indicators.begin();
  for(unsigned int i=0;i<rad_boundary_indicators.size();i++) {
	  std::cout << (unsigned int)*it << " ";
	  std::advance(it, 1);
  }
  std::cout << std::endl;

  boundary7.insert(7);
  std::cout << "# AJTL radiation boundary indicators: ";
  it = boundary7.begin();
  for(unsigned int i=0;i<boundary7.size();i++) {
	  std::cout << (unsigned int)*it << " ";
	  std::advance(it, 1);
  }
  std::cout << std::endl;

  // --- Output mesh to .eps file ---
  //std::ofstream out ( geometry_type + ".eps");
  //GridOut grid_out;
  //grid_out.write_eps (triangulation, out);*/

}


  template <int dim>
  void FFO<dim>::set_bc (const double h_field, const double gamma, const double htffo, const double hx) {
  if(geometry_type=="tffo")
		set_bc_tffo_sharp(h_field,gamma,htffo,hx);
  else
		std::cout << "# Error: geometry_type is not set " << std::endl;
}


  // Initialize voltage filter
  template <int dim>
  void FFO<dim>::set_voltage_filters () {

  opt_filter_init(voltage_filter_1);
  opt_filter_init(voltage_filter_3);
  opt_filter_init(voltage_filter_5);

}


  template <int dim>
  void FFO<dim>::copy_state_from (const FFO<dim>& object) {

	  if(ffo_tunnel_current->Nnodes != object.ffo_tunnel_current->Nnodes ||
		 ffo_tunnel_current->Nexps != object.ffo_tunnel_current->Nexps)
		  	  throw "# Error: array sizes do not match";

	  U_current = object.U_current;
	  U_old = object.U_old;

	  memcpy( ffo_tunnel_current->memstate.sin05phi_old,
			  object.ffo_tunnel_current->memstate.sin05phi_old,
			  ffo_tunnel_current->Nnodes * sizeof(double) );

	  memcpy( ffo_tunnel_current->memstate.cos05phi_old,
			  object.ffo_tunnel_current->memstate.cos05phi_old,
			  ffo_tunnel_current->Nnodes * sizeof(double) );

	  memcpy( ffo_tunnel_current->memstate.F, object.ffo_tunnel_current->memstate.F,
			  ffo_tunnel_current->Nnodes * ffo_tunnel_current->Nexps *
			  sizeof(double _Complex) );

	  memcpy( ffo_tunnel_current->memstate.G, object.ffo_tunnel_current->memstate.G,
			  ffo_tunnel_current->Nnodes * ffo_tunnel_current->Nexps *
			  sizeof(double _Complex) );

}

//================== PRIVATE METHODS =====================

template <int dim>
void FFO<dim>::set_bc_ffo_rect (const double h_field, const double gamma)
{

  boundary_forces.reinit (dof_handler.n_dofs());
  Vector<double> tmp_vector (dof_handler.n_dofs());

  std::set<types::boundary_id> boundary_indices;

  boundary_indices.insert(1);
  VectorTools::create_boundary_right_hand_side (dof_handler,
                                                QGauss<dim-1>(gauss_degree),
                                                ConstantFunction<dim>(h_field),
                                                tmp_vector,
                                                boundary_indices);
  boundary_forces += tmp_vector;

  boundary_indices.clear();
  boundary_indices.insert(2);
  VectorTools::create_boundary_right_hand_side (dof_handler,
                                                QGauss<dim-1>(gauss_degree),
                                                ConstantFunction<dim>(-h_field),
                                                tmp_vector,
                                                boundary_indices);
  boundary_forces += tmp_vector;

  boundary_indices.clear();
  boundary_indices.insert(3);
  VectorTools::create_boundary_right_hand_side (dof_handler,
                                                QGauss<dim-1>(gauss_degree),
                                                ConstantFunction<dim>(gamma*W0/2.),
                                                tmp_vector,
                                                boundary_indices);
  boundary_forces += tmp_vector;

}



 template <int dim>
  void FFO<dim>::set_bc_ffo_sharp (const double h_field, const double gamma)
  {

    boundary_forces.reinit (dof_handler.n_dofs());
    Vector<double> tmp_vector (dof_handler.n_dofs());

    std::set<types::boundary_id> boundary_indices;

    // left boundary
    boundary_indices.insert(1);
    VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>(h_field),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

    // right boundary
    boundary_indices.clear();
    boundary_indices.insert(2);
    VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>(-h_field),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

	const double costh0=(W0-P0)/std::sqrt((W0-P0)*(W0-P0)+4.*S0*S0);
	double sinth0=2.*S0/std::sqrt((W0-P0)*(W0-P0)+4.*S0*S0);
	double hg=0.5*gamma*Area/L0;

	// other boundaries
    boundary_indices.clear();
    boundary_indices.insert(3);
    VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>(hg),
                                                  tmp_vector,
                                                  boundary_indices);

    boundary_forces += tmp_vector;

    boundary_indices.clear();
    boundary_indices.insert(4);
	VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>(hg*sinth0+h_field*costh0),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

    boundary_indices.clear();
    boundary_indices.insert(5);
	VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>(hg*sinth0-h_field*costh0),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

  }


 template <int dim>
  void FFO<dim>::set_bc_tffo_sharp (const double h_field, const double gamma, const double htffo, const double hx)
  {

    boundary_forces.reinit (dof_handler.n_dofs());
    Vector<double> tmp_vector (dof_handler.n_dofs());

    std::set<types::boundary_id> boundary_indices;

    boundary_indices.insert(0);
	VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>(htffo),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

    boundary_indices.clear();
    boundary_indices.insert(1);
    VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>(h_field),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

    boundary_indices.clear();
    boundary_indices.insert(2);
    VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>(-h_field),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

	const double costh0=(W0-P0)/std::sqrt((W0-P0)*(W0-P0)+4.*S0*S0);
	double sinth0=2.*S0/std::sqrt((W0-P0)*(W0-P0)+4.*S0*S0);
	double hg=0.5*gamma*Area/L0;

    boundary_indices.clear();
    boundary_indices.insert(3);
    VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>(hg+hx),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

    boundary_indices.clear();
   	boundary_indices.insert(4);
   	VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>((hg+hx)*sinth0+h_field*costh0),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

    boundary_indices.clear();
    boundary_indices.insert(5);
	VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>((hg+hx)*sinth0-h_field*costh0),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

    // htffo: independent TFFO side bias field
	const double costh1=(W1-P1)/std::sqrt((W1-P1)*(W1-P1)+4.*S1*S1);
	const double sinth1=2.*S1/std::sqrt((W1-P1)*(W1-P1)+4.*S1*S1);

    boundary_indices.clear();
    boundary_indices.insert(6);
	VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>((hg-hx)*costh1+htffo*sinth1),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

    boundary_indices.clear();
    boundary_indices.insert(7);
    boundary_indices.insert(8);
	VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>(hg-hx),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

    boundary_indices.clear();
   	boundary_indices.insert(9);
   	VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>((hg-hx)*sinth0+h_field*costh0),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;

    boundary_indices.clear();
    boundary_indices.insert(10);
	VectorTools::create_boundary_right_hand_side (dof_handler,
                                                  QGauss<dim-1>(gauss_degree),
                                                  ConstantFunction<dim>((hg-hx)*sinth0-h_field*costh0),
                                                  tmp_vector,
                                                  boundary_indices);
    boundary_forces += tmp_vector;
  }



  template <int dim>
    void FFO<dim>::setup_system ()
    {
      dof_handler.distribute_dofs (fe);

      sparsity_pattern.reinit (dof_handler.n_dofs(),
                               dof_handler.n_dofs(),
                               dof_handler.max_couplings_between_dofs());

      DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

      sparsity_pattern.compress ();

      mass_matrix.reinit    (sparsity_pattern);
      laplace_matrix.reinit (sparsity_pattern);
      system_matrix.reinit    (sparsity_pattern);
      new_matrix.reinit  (sparsity_pattern);
      c_matrix.reinit  (sparsity_pattern);
      old_matrix.reinit  (sparsity_pattern);

      U_old.reinit       (dof_handler.n_dofs());
      U_current.reinit       (dof_handler.n_dofs());
      U_next.reinit       (dof_handler.n_dofs());
      system_rhs.reinit     (dof_handler.n_dofs());

      MatrixCreator::create_mass_matrix (dof_handler,
                                         QGauss<dim>(gauss_degree),
                                         mass_matrix);
      MatrixCreator::create_laplace_matrix (dof_handler,
                                            QGauss<dim>(gauss_degree),
                                            laplace_matrix);

    }



  template <int dim>
  void FFO<dim>::display_info () const
  {
  std::cout << "#================================================="<< std::endl
			<< "# h_field: " << h0_field  << std::endl
			<< "#================================================="<< std::endl
			<< "# MiTMoJCo OMP_NUM_THREADS: " << OMP_NUM_THREADS << std::endl
			<< "# time step: " << time_step << std::endl
			<< "# beta: " << beta0  << std::endl
			<< "# " << std::endl
			<< "# geometry_type: " << geometry_type << std::endl
			<< "# L0: " << L0 << std::endl
			<< "# W0: " << W0 << std::endl
			<< "# S0: " << S0 << std::endl
			<< "# P0: " << P0 << std::endl
            << "# L1: " << L1 << std::endl
		    << "# W1: " << W1 << std::endl
		    << "# S1: " << S1 << std::endl
		    << "# P1: " << P1 << std::endl;
  std::cout << "# Area: " << Area << std::endl
			<< "# active cells: " << triangulation.n_active_cells() << std::endl
//			<< "# total number of cells: " << triangulation.n_cells() << std::endl
			<< "# number of dofs: " << dof_handler.n_dofs() << std::endl;

  if(vtk_tmax>vtk_tmin)
		  std::cout	<< "# vtk_tmin: " << vtk_tmin  << std::endl
					<< "# vtk_tmax: " << vtk_tmax  << std::endl;
  }



  template <int dim>
  void FFO<dim>::assemble_matrices ()
  {
    new_matrix.copy_from (mass_matrix);
    new_matrix.add (0.5*alpha0*time_step, mass_matrix);
    new_matrix.add (0.5*beta0*time_step, laplace_matrix);

    c_matrix.copy_from (mass_matrix);
    c_matrix.add (1.0, mass_matrix);
    c_matrix.add (-time_step*time_step, laplace_matrix);

    old_matrix.copy_from (mass_matrix);
    old_matrix*=-1;
    old_matrix.add (0.5*alpha0*time_step, mass_matrix);
    old_matrix.add (0.5*beta0*time_step, laplace_matrix);
  }



  template <int dim>
  void FFO<dim>::set_ic (const double h_field, const double gamma)
  {
	std::cout << "# Starting from initial conditions"  << std::endl;
    const InitialValues<dim> initial_values (geometry_type, h_field, gamma, L0);
    ConstraintMatrix constraints;
    constraints.close();
    VectorTools::project (dof_handler,
                            constraints,
                            QGauss<dim>(gauss_degree),
                            initial_values,
                            U_old);
    U_current=U_old;

    mitmojco_init( ffo_tunnel_current );  // Initialize the tunnel current object

  }


  // Used for sine-Gordon test
  template <int dim>
  void FFO<dim>::compute_sin_vector (const Vector<double> &argument,
                                                Vector<double> &sin_vector) const
  {
    const QGauss<dim> quadrature_formula (gauss_degree);
    FEValues<dim>     fe_values (fe, quadrature_formula,
                                 update_values |
                                 update_JxW_values |
                                 update_quadrature_points);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    Vector<double> local_sin_vector (dofs_per_cell); // initialized with zeros
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
    std::vector<double> argument_values (n_q_points);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    sin_vector=0; // set to zero
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        fe_values.get_function_values (argument, argument_values);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            local_sin_vector(i) += std::sin( argument_values[q_point] ) *
                                 fe_values.shape_value (i, q_point) *
                                 fe_values.JxW (q_point);

        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          sin_vector(local_dof_indices[i]) += local_sin_vector(i);

        local_sin_vector = 0;
      }
  }



  template <int dim>
  unsigned int
  FFO<dim>::solve ()
  {
    SolverControl solver_control (1000, 1e-12*system_rhs.l2_norm());

    SolverCG<> cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(new_matrix, 1.0);

    cg.solve (new_matrix, U_next,
              system_rhs,
              preconditioner);

    return solver_control.last_step();
  }


  template <int dim>
  void FFO<dim>::vtk_output (const double t) const
  {
    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (U_current, "u");
    data_out.build_patches ();

    const std::string filename =  "t" + Utilities::int_to_string ( (int) round(t*10), 5) + ".vtk";

    std::ofstream output (filename.c_str());
    data_out.write_vtk (output);
  }


  template <int dim>
  int FFO<dim>::advance_by (unsigned int n_steps)
  {
  #ifdef timer_on
    Timer timer;
  #endif

  Vector<double>  tmp_vector (dof_handler.n_dofs());
  Vector<double>  jbar_vector (dof_handler.n_dofs());

  /*double sum;
  double U_mean_start, U_mean_finish;
  sum=0.;
  for(unsigned int i=0; i<dof_handler.n_dofs(); i++)
	  sum += U_current[i];
  U_mean_start=sum/dof_handler.n_dofs();*/

  for (unsigned int step=0; step < n_steps; step++) {

    	  // system_rhs
    	  c_matrix.vmult (system_rhs, U_current);
    	  old_matrix.vmult (tmp_vector, U_old);
    	  system_rhs.add (1.0, tmp_vector);
    	  system_rhs.add (time_step*time_step, boundary_forces);

    	  // Update the tunnel current
    	  mitmojco_update( ffo_tunnel_current );
    	  memcpy(tmp_vector.begin(),ffo_tunnel_current->jbar, dof_handler.n_dofs()*sizeof(double));
    	  mass_matrix.vmult(jbar_vector,tmp_vector);
    	  system_rhs.add (-time_step*time_step, jbar_vector);

    	  // sine-Gordon
    	  /*compute_sin_vector(U_current, tmp_vector);
    	  system_rhs.add (-time_step*time_step, tmp_vector);*/

    	  solve(); // U_next is calculated

   		  double U_boundary_old=0.;
   		  double U_boundary_new=0.;
   		  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   			  if (filter_boundary_dofs[i] == true) {
   				  U_boundary_old+=U_old[i];
   				  U_boundary_new+=U_next[i];
   			  }
   		  // Make a record to the filter object
   		  opt_filter_update(voltage_filter_1, U_boundary_new - U_boundary_old);
   		  opt_filter_update(voltage_filter_3, U_boundary_new - U_boundary_old);
   		  opt_filter_update(voltage_filter_5, U_boundary_new - U_boundary_old);

    	  U_old=U_current;
    	  U_current=U_next;

  } // end of time loop

  /*sum=0.;
  for(unsigned int i=0; i<dof_handler.n_dofs(); i++)
	  sum += U_current[i];
  U_mean_finish=sum/dof_handler.n_dofs();
  V_dc = 0.5*(U_mean_finish-U_mean_start)/(n_steps*time_step);*/

	/* DC voltage in units hbar*omega_J/e.	Factor 0.25 comes from the Josephson relation
		Vdc = 0.5*dphi/dt, and 2nd order discretization for the derivative dphi/dt */
  V_dc_f1 = 0.25*opt_filter_result(voltage_filter_1)/(time_step*n_filter_boundary_dofs); // equivalent to the arithmetic mean
  V_dc_f3 = 0.25*opt_filter_result(voltage_filter_3)/(time_step*n_filter_boundary_dofs);
  V_dc_f5 = 0.25*opt_filter_result(voltage_filter_5)/(time_step*n_filter_boundary_dofs);

  return 0;
  }


  // Input: 
  // 	'V': output voltages
  // 	'P': output voltages and phase space variables
  template <int dim>
  int FFO<dim>::advance_by (unsigned int n_steps, char output_regime)
  {
  #ifdef timer_on
    Timer timer;
  #endif

  Vector<double>  tmp_vector (dof_handler.n_dofs());
  Vector<double>  jbar_vector (dof_handler.n_dofs());

  /*double sum;
  double U_mean_start, U_mean_finish;
  sum=0.;
  for(unsigned int i=0; i<dof_handler.n_dofs(); i++)
	  sum += U_current[i];
  U_mean_start=sum/dof_handler.n_dofs();*/

  for (unsigned int step=0; step < n_steps; step++) {

    	  // system_rhs
    	  c_matrix.vmult (system_rhs, U_current);
    	  old_matrix.vmult (tmp_vector, U_old);
    	  system_rhs.add (1.0, tmp_vector);
    	  system_rhs.add (time_step*time_step, boundary_forces);

    	  // Update the tunnel current
    	  mitmojco_update( ffo_tunnel_current );
    	  memcpy(tmp_vector.begin(),ffo_tunnel_current->jbar, dof_handler.n_dofs()*sizeof(double));
    	  mass_matrix.vmult(jbar_vector,tmp_vector);
    	  system_rhs.add (-time_step*time_step, jbar_vector);

    	  // sine-Gordon
    	  /*compute_sin_vector(U_current, tmp_vector);
    	  system_rhs.add (-time_step*time_step, tmp_vector);*/

    	  solve(); // U_next is calculated

   		  double U_boundary_old=0.;
   		  double U_boundary_new=0.;
   		  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   			  if (filter_boundary_dofs[i] == true) {
   				  U_boundary_old+=U_old[i];
   				  U_boundary_new+=U_next[i];
   			  }
   		  // Make a record to the filter object
   		  opt_filter_update(voltage_filter_1, U_boundary_new - U_boundary_old);
   		  opt_filter_update(voltage_filter_3, U_boundary_new - U_boundary_old);
   		  opt_filter_update(voltage_filter_5, U_boundary_new - U_boundary_old);

   		  if(output_regime=='V') {

   			  U_boundary_old=0.;
   			  U_boundary_new=0.;
   			  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   				  if (rad_boundary_dofs[i] == true) {
   					  U_boundary_old+=U_old[i];
   					  U_boundary_new+=U_next[i];
   				  }
   			  double Voltage_current = 0.25*(U_boundary_new - U_boundary_old)/(time_step*n_rad_boundary_dofs);
   			  std::cout << Voltage_current;

   	   		  if(geometry_type == "tffo") {
   	   			  U_boundary_old=0.;
   	   			  U_boundary_new=0.;
   	   			  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   	   				  if (boundary7_dofs[i] == true) {
   	   					  U_boundary_old+=U_old[i];
   	   					  U_boundary_new+=U_next[i];
   	   				  }
   	   			  double Voltage_current = 0.25*(U_boundary_new - U_boundary_old)/(time_step*n_boundary7_dofs);
   	   			  std::cout << "  " << Voltage_current;
   	   		  }

   	   		  std::cout << std::endl;
   		  }
   		  else if(output_regime=='P') {

   			  U_boundary_old=0.;
   			  double U_boundary_current=0.;
   			  U_boundary_new=0.;
   			  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   				  if (filter_boundary_dofs[i] == true) {
   					  U_boundary_old+=U_old[i];
                      U_boundary_current+=U_current[i];
   					  U_boundary_new+=U_next[i];
   				  }
              double phi1 = U_boundary_current/n_filter_boundary_dofs;
   			  double phi1_dot = 0.5*(U_boundary_new - U_boundary_old)/(time_step*n_filter_boundary_dofs);
   			  //double Voltage_at_1 = 0.5*phi1_dot;

   			  U_boundary_old=0.;
   			  U_boundary_current=0.;
   			  U_boundary_new=0.;
   			  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
   				  if (rad_boundary_dofs[i] == true) {
   					  U_boundary_old+=U_old[i];
                      U_boundary_current+=U_current[i];
   					  U_boundary_new+=U_next[i];
   				  }
              double phi2 = U_boundary_current/n_rad_boundary_dofs;
   			  double phi2_dot = 0.5*(U_boundary_new - U_boundary_old)/(time_step*n_rad_boundary_dofs);
   			  double Voltage_at_2 = 0.5*phi2_dot;

   	   		  if(geometry_type == "tffo") {
       			  U_boundary_old=0.;
       			  U_boundary_current=0.;
       			  U_boundary_new=0.;
       			  for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
       				  if (boundary7_dofs[i] == true) {
       					  U_boundary_old+=U_old[i];
                          U_boundary_current+=U_current[i];
       					  U_boundary_new+=U_next[i];
       				  }
                  double phi3 = U_boundary_current/n_boundary7_dofs;
       			  double phi3_dot = 0.5*(U_boundary_new - U_boundary_old)/(time_step*n_boundary7_dofs);
   	   			  double Voltage_at_3 = 0.5*phi3_dot;
   	   			  std::cout << Voltage_at_2 << "  " << Voltage_at_3
							<< "  " << phi2 - phi1 << "  " << phi2_dot - phi1_dot
       			  			<< "  " << phi3 - phi1 << "  " << phi3_dot - phi1_dot 
							<< "  " << phi3 - phi2 << "  " << phi3_dot - phi2_dot << std::endl;
   	   		  } 
              else
       			  std::cout << Voltage_at_2 
							<< "  " << phi2 - phi1 << "  " << phi2_dot - phi1_dot << std::endl;

   		  }

   		  /*if(time>=vtk_tmin && time<=vtk_tmax)
       		  if (count % vtk_skip == 0)
       			  vtk_output (gamma,time);*/

    	  U_old=U_current;
    	  U_current=U_next;

  } // end of time loop

  /*sum=0.;
  for(unsigned int i=0; i<dof_handler.n_dofs(); i++)
	  sum += U_current[i];
  U_mean_finish=sum/dof_handler.n_dofs();
  V_dc = 0.5*(U_mean_finish-U_mean_start)/(n_steps*time_step);*/

	/* DC voltage in units hbar*omega_J/e.	Factor 0.25 comes from the Josephson relation
		Vdc = 0.5*dphi/dt, and 2nd order discretization for the derivative dphi/dt */
  V_dc_f1 = 0.25*opt_filter_result(voltage_filter_1)/(time_step*n_filter_boundary_dofs); // equivalent to the arithmetic mean
  V_dc_f3 = 0.25*opt_filter_result(voltage_filter_3)/(time_step*n_filter_boundary_dofs);
  V_dc_f5 = 0.25*opt_filter_result(voltage_filter_5)/(time_step*n_filter_boundary_dofs);

  return 0;
  }


  template <int dim>
  void FFO<dim>::shift_phase()
  {
      const double shift = 4.*M_PI * int(U_current[0] /(4.*M_PI));
      for(unsigned int i=0; i<dof_handler.n_dofs(); i++) {
    	  U_old[i] = U_old[i] - shift;
    	  U_current[i] = U_current[i] - shift;
      }
  }


  // Distance between two states: root mean square deviation at different points
  template <int dim>
  double FFO<dim>::distance_from(FFO<dim>& object) {

  Vector<double> dU = U_current;
  dU -= object.U_current;

  Vector<double> dV = U_current;
  dV -= U_old;
  dV -= object.U_current;
  dV += object.U_old;
  dV *= 1./time_step;

  return sqrt( (dU.norm_sqr()+dV.norm_sqr())/(2.*dof_handler.n_dofs()) );
  }


}   //// close namespace

// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------



#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

using namespace std;
using namespace dealii;
using namespace FluxFlow;

void run_lyapunov(double mfield, double settling_time, double integration_time, double time_step,
		double gamma_start, double gamma_finish, double gamma_step, double htffo, double hx, bool distout);

void run_ffo(double mfield, double settling_time, double integration_time, double time_step,
		double gamma_start, double gamma_finish, double gamma_step, double htffo, double hx);

int main (int argc, char* argv[])
{
  try
    {
    //int optind=1;
    int c;

	opterr = 0;
	while ((c = getopt (argc, argv, "f:")) != -1)
    	switch (c)
    	{
    	case '?':
                if (isprint (optopt))
                  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                  fprintf (stderr,
                           "Unknown option character `\\x%x'.\n",
                           optopt);
                return 1;
    	default:
                abort ();
	    }

	  double mfield;
	  double settling_time, integration_time, time_step;
	  double gamma_start, gamma_finish, gamma_step;
      double htffo, hx;

	  if(argc-optind==5) {
		  mfield = atof(argv[optind]);
		  settling_time = atof(argv[optind+1]);
		  integration_time = atof(argv[optind+2]);
		  time_step = atof(argv[optind+3]);
		  gamma_start = atof(argv[optind+4]);
		  gamma_finish = gamma_start;
		  gamma_step = 0.01; // arbitrary
          htffo=0.;
          hx=0.;
		  std::cout << "# --- Input arguments ---" << std::endl;
		  std::cout << "# magnetic field: " << mfield << std::endl;
		  std::cout << "# settling_time: " << settling_time << std::endl;
		  std::cout << "# integration_time: " << integration_time << std::endl;
		  std::cout << "# time_step: " << time_step << std::endl;
		  std::cout << "# gamma: " << gamma_start << std::endl;
		  std::cout << "# htffo: " << htffo << std::endl;
		  std::cout << "# hx: " << hx << std::endl;
	  } else if(argc-optind==9) {
		  mfield = atof(argv[optind]);
		  settling_time = atof(argv[optind+1]);
		  integration_time = atof(argv[optind+2]);
		  time_step = atof(argv[optind+3]);
		  gamma_start = atof(argv[optind+4]);
		  gamma_finish = atof(argv[optind+5]);
		  gamma_step = atof(argv[optind+6]);
		  htffo = atof(argv[optind+7]);
          hx = atof(argv[optind+8]);
		  std::cout << "# --- Input arguments ---" << std::endl;
		  std::cout << "# magnetic field: " << mfield << std::endl;
		  std::cout << "# settling_time: " << settling_time << std::endl;
		  std::cout << "# integration_time: " << integration_time << std::endl;
		  std::cout << "# time_step: " << time_step << std::endl;
		  std::cout << "# gamma_start: " << gamma_start << std::endl;
		  std::cout << "# gamma_finish: " << gamma_finish << std::endl;
		  std::cout << "# gamma_step: " << gamma_step << std::endl;
		  std::cout << "# htffo: " << htffo << std::endl;
		  std::cout << "# hx: " << hx << std::endl;
	  } else {
		  std::cout << "# Incorrect number of arguments." << std::endl;
		  std::cout << "# Please, supply 5 or 9 arguments in the following order:" << std::endl;
		  std::cout << "# 1. magnetic field" << std::endl;
		  std::cout << "# 2. settling_time" << std::endl;
		  std::cout << "# 3. integration_time" << std::endl;
		  std::cout << "# 4. time_step"  << std::endl;
		  std::cout << "# 5. gamma_start" << std::endl;
		  std::cout << "# 6. gamma_finish" << std::endl;
		  std::cout << "# 7. gamma_step"  << std::endl;
		  std::cout << "# 8. htffo" << std::endl;
		  std::cout << "# 9. hx" << std::endl;
		  return -1;
	  }

      deallog.depth_console (0);

      omp_set_num_threads(OMP_NUM_THREADS); // number of OpenMP threads

      gamma_step = (gamma_start<= gamma_finish)?gamma_step:(-gamma_step);

      Timer timer;
      timer.start ();

      /* Calculation of IVC of T-junction FFO */  
      run_ffo( mfield, settling_time, integration_time, time_step,
   		  	  	  	  	  	  gamma_start, gamma_finish, gamma_step, htffo, hx);


      /* Calculation of Lyapunov exponents for T-junction FFO */  
//      run_lyapunov( mfield, settling_time, integration_time, time_step,
//  		  	  	  	  	  	  gamma_start, gamma_finish, gamma_step, htffo, hx, 0);



      timer.stop ();

      std::cout << "# Wall time: " << std::setprecision(3)
				<< timer.wall_time() << " seconds."
				<< "  User CPU time: " << std::setprecision(3)
				<< timer() << " seconds." << std::endl;

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "#----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "#----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------Unknown exception------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------Unknown exception------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}




void run_lyapunov(double mfield, double settling_time, double integration_time, double time_step,
		double gamma_start, double gamma_finish, double gamma_step, double htffo, double hx, bool distout) {

    FFO<2> ffo1(mfield, time_step);
    FFO<2> ffo2(mfield, time_step);

    double gamma;
    double V_dc_max=0.8*ffo1.kgap;

    std::cout << "#-------------------------------------------------"<< std::endl;
    std::cout << "# 0 gamma"<< std::endl;
    std::cout << "# 1 V_dc_f1"<< std::endl;
    std::cout << "# 2 V_dc_f3"<< std::endl;
    std::cout << "# 3 V_dc_f5"<< std::endl;

    double rescale_dist = 1.e-10;
    srand(time(NULL));
    //srand(1234);
    unsigned int n_steps = 100;
    unsigned int step_count;
    unsigned int n_settling_steps = round(settling_time/time_step);
    unsigned int n_lyapunov_settling_steps = round(50./time_step);
    unsigned int n_integration_steps = round(integration_time/time_step);

    std::cout << "# n_steps: " << n_steps << std::endl;
    std::cout << "# n_settling_steps: " << n_settling_steps << std::endl;
    std::cout << "# n_lyapunov_settling_steps: " << n_lyapunov_settling_steps << std::endl;
    std::cout << "# n_integration_steps: " << n_integration_steps << std::endl;

    ffo1.set_ic(mfield, gamma_start);

    for(gamma=gamma_start; (((gamma < gamma_finish + 0.5*gamma_step) && (gamma_step>0)) ||
    		((gamma > gamma_finish + 0.5*gamma_step) && (gamma_step<0))) && ffo1.V_dc_f1 <= V_dc_max; gamma+=gamma_step) {

   	  // set boundary conditions
      ffo1.set_bc(mfield, gamma, htffo, hx);
      ffo2.set_bc(mfield, gamma, htffo, hx);

      // preliminary run of ffo1 for settling time interval
      ffo1.advance_by(n_settling_steps);

      // start recording the voltage
  	  ffo1.set_voltage_filters();

        // synchronize ffo2 to ffo1 and introduce a small distortion to ffo2
  	  ffo2.copy_state_from(ffo1);
  	  double noise_amp = 1.e-7;
        for (unsigned int i=0; i<ffo2.dof_handler.n_dofs(); ++i) {
      	  double r = 2.*((double)rand()/(double)RAND_MAX - 0.5); // random number in (-1,1)
      	  ffo2.U_current[i] += 2.*noise_amp*r;
        }

  	  double t=0.;

        // initial deviation between ffo1 and ffo2
  	  double dist = ffo1.distance_from(ffo2);
  	  double dist_initial=dist;

  	  // initial run: estimate the direction of Lyapunov vector
  	  step_count = 0;
      do {

      	  if(distout) {
      		  std::cout << t << "  " << dist << std::endl;
      		  t += n_steps*time_step;
      	  }

      	  step_count += n_steps;

      	  // run ffo1 (without distortion)
      	  ffo1.advance_by(n_steps);

      	  // run ffo2 (with distortion)
      	  ffo2.advance_by(n_steps);

      	  // calculate deviation between ffo1 and ffo2
      	  dist = ffo1.distance_from(ffo2);

      } while(dist<0.1 && step_count < n_lyapunov_settling_steps);

        // performance run: calculate max Lyapunov exponent
        unsigned int count=0;
        unsigned int integration_step=0;
        double lyapunov_average=0.;
        do {

      	  dist = ffo1.distance_from(ffo2);
      	  double epsilon = rescale_dist/dist;

      	  // shift phases of ffo1 and ffo2 by the same amount
            const double shift = 4.*M_PI * int(ffo1.U_current[0] /(4.*M_PI));
            for(unsigned int i=0; i<ffo1.dof_handler.n_dofs(); i++) {
          	  ffo1.U_old[i] -= shift;
          	  ffo1.U_current[i] -= shift;
          	  ffo2.U_old[i] -= shift;
          	  ffo2.U_current[i] -= shift;
            }

            // rescale perturbation in the direction of the Lyapunov vector by epsilon
      	  ffo2.U_current *= epsilon;
      	  ffo2.U_current.add(1.-epsilon, ffo1.U_current);
      	  ffo2.U_old *= epsilon;
      	  ffo2.U_old.add(1.-epsilon, ffo1.U_old);
      	  unsigned int Nexps = ffo2.ffo_tunnel_current->Nexps;
      	  for(unsigned int i=0; i<ffo2.dof_handler.n_dofs(); i++) {
      		  ffo2.ffo_tunnel_current->memstate.sin05phi_old[i] *= epsilon;
      		  ffo2.ffo_tunnel_current->memstate.sin05phi_old[i] += (1.-epsilon)*ffo1.ffo_tunnel_current->memstate.sin05phi_old[i];
        		  ffo2.ffo_tunnel_current->memstate.cos05phi_old[i] *= epsilon;
        		  ffo2.ffo_tunnel_current->memstate.cos05phi_old[i] += (1.-epsilon)*ffo1.ffo_tunnel_current->memstate.cos05phi_old[i];
        		  for(unsigned int n=0;n<Nexps;n++) {
        			  ffo2.ffo_tunnel_current->memstate.F[i*Nexps+n] *= epsilon;
        			  ffo2.ffo_tunnel_current->memstate.F[i*Nexps+n] += (1.-epsilon)*ffo1.ffo_tunnel_current->memstate.F[i*Nexps+n];
        			  ffo2.ffo_tunnel_current->memstate.G[i*Nexps+n] *= epsilon;
        			  ffo2.ffo_tunnel_current->memstate.G[i*Nexps+n] += (1.-epsilon)*ffo1.ffo_tunnel_current->memstate.G[i*Nexps+n];
        		  }
      	  }

  		  // check the distance
      	  dist = ffo1.distance_from(ffo2);
      	  dist_initial=dist;

      	  step_count = 0;
      	  do {
          	  if(distout) {
          		  std::cout << t << "  " << dist << std::endl;
      		  	  t += n_steps*time_step;
          	  }

      		  step_count += n_steps;
      		  integration_step += n_steps;

      		  ffo1.advance_by(n_steps);
      		  ffo2.advance_by(n_steps);
      		  dist = ffo1.distance_from(ffo2);

      	  } while(dist<0.1 && integration_step < n_integration_steps);

      	  double lyapunov_exp = log(dist/dist_initial)/(step_count*time_step);
      	  std::cout << "# count: " << count << "  lyapunov_exp: " << lyapunov_exp << std::endl;

      	  if(integration_step < n_integration_steps) {
      		  lyapunov_average += lyapunov_exp;
      	  	  count++;
      	  }
      	  else if(count==0) {
      		  lyapunov_average = lyapunov_exp;
      	  	  count = 1;
      	  }

        } while(integration_step < n_integration_steps);

        std::cout << gamma << "  " << ffo1.V_dc_f1 << "  "  << ffo1.V_dc_f3 << "  " << ffo1.V_dc_f5 << "  " << lyapunov_average/count << std::endl;

    } // end of gamma loop

}


void run_ffo(double mfield, double settling_time, double integration_time, double time_step,
		double gamma_start, double gamma_finish, double gamma_step, double htffo, double hx) {

    FFO<2> ffo(mfield, time_step);

    double gamma;
    double V_dc_max=0.8*ffo.kgap;

    std::cout << "#-------------------------------------------------"<< std::endl;
    std::cout << "# 0 gamma"<< std::endl;
    std::cout << "# 1 V_dc_f1"<< std::endl;
    std::cout << "# 2 V_dc_f3"<< std::endl;
    std::cout << "# 3 V_dc_f5"<< std::endl;

    unsigned int n_settling_steps = round(settling_time/time_step);
    unsigned int n_integration_steps = round(integration_time/time_step);

    std::cout << "# n_settling_steps: " << n_settling_steps << std::endl;
    std::cout << "# n_integration_steps: " << n_integration_steps << std::endl;

    ffo.set_ic(mfield, gamma_start);

    for(gamma=gamma_start; (((gamma < gamma_finish + 0.5*gamma_step) && (gamma_step>0)) ||
    		((gamma > gamma_finish + 0.5*gamma_step) && (gamma_step<0))) && ffo.V_dc_f1 <= V_dc_max; gamma+=gamma_step) {

  	  // set boundary conditions
      ffo.set_bc(mfield, gamma, htffo, hx);

      // preliminary run of ffo for settling time interval
      ffo.advance_by(n_settling_steps);

      // start recording the voltage
  	  ffo.set_voltage_filters();

  	  ffo.advance_by(n_integration_steps, 'V');

      std::cout << "# " << gamma << "  " << ffo.V_dc_f1 << "  " << ffo.V_dc_f3 << "  " << ffo.V_dc_f5 << std::endl;

      ffo.shift_phase();

    } // end of gamma loop

}

