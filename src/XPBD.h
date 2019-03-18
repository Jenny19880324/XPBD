#pragma once
#include <igl/igl_inline.h>
#include "XPBDEnergyType.h"
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unordered_set>

namespace xpbd
{
	struct XPBDData
	{
		// n               #V
		// energy          type of energy to use
		// with_dynamics   whether using dynamics
		// f_ext           #V by dim list of external forces
		// vel             #V by dim list of velocities
		// h               dynamics time step
		// ym              Young's modulus smaller is softer, larger is more rigid/stiff
		// max_iter        maximum inner iterations
		// M               mass matrix
		// b               list of boundary indices into V
		// dim             dimension being used for solving

		int n;
		int dim;
		int max_iter;
		XPBDEnergyType energy;
		double h;
		double ym;
		
		Eigen::MatrixXd f_ext, vel, V;
		Eigen::VectorXd lambda_iter;
		Eigen::SparseMatrix<double> M, M_inv;
		std::unordered_set<int> b;
		
		XPBDData() :
			n(0),
			max_iter(20),
			energy(XPBD_ENERGY_TYPE_DEFAULT),
			f_ext(),
			h(0.01),
			ym(1),
			b(),
			dim(2)
		{
		};
	};

	// Compute necessary information to start using xpbd 
	//
	// Inputs:
	//    V     #V by dim list of mesh positions
	//    F     #F by simplex-size list of triangle|tet indices into V
	//    dim   dimension being used at solve time.
	//          For deformation usually dim = V.cols(), 
	//          for surface parameterization V.cols() = 3 and dim = 2
	//    b     #b list of "boundary" fixed vertex indices into V
	// Outputs:
	//   data   struct containing necessary precomputation
	template<
		typename DerivedV,
		typename DerivedF,
		typename DerivedE>
		IGL_INLINE bool xpbd_precomputation(
			const Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F,
			const Eigen::PlainObjectBase<DerivedE> & E,
			const int dim,
			const std::unordered_set<int> & b,
			XPBDData & data);


	// update vertex positions with a solver in Gauss-Seidel fashion
	// Algorithm 1 XPBD simulation loop described in
	// XPBD: Position-Based Simulation of Compliant Constrained Dynamics
	// http://mmacklin.com/xpbd.pdf
	//
	// Inputs:
	//    V     #V by dim initial guess
	//    F     #F by simplex-size list of triangle|tet indices into V
	//    dim   dimension being used at solve time.
	//          For deformation usually dim = V.cols(),
	//          for surface parameterization V.cols() = 3 and dim 2
	//    b     #b list of "boundary" fixed vertex indices into V
	//    bc    #b by dim list of boundary conditions
	template <
		typename DerivedV,
		typename DerivedE>
		IGL_INLINE bool xpbd_solve(
			Eigen::PlainObjectBase<DerivedV> & U,
			const Eigen::PlainObjectBase<DerivedE> & E,
			XPBDData & data
			);

} // namespace XPBD

#ifndef IGL_STATIC_LIBRARY
#include "xpbd.cpp"
#endif