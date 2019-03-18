#pragma once
#include <igl/igl_inline.h>
#include "XPBDEnergyType.h"
#include <Eigen/Core>
#include <Eigen/Sparse>

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

		bool with_dynamics;
		int n;
		int dim;
		int max_iter;
		XPBDEnergyType energy;
		double h;
		double ym;
		
		Eigen::MatrixXd f_ext, vel;
		Eigen::VectorXi b;
		
		XPBDData() :
			n(0),
			energy(XPBD_ENERGY_TYPE_DEFAULT),
			with_dynamics(false),
			f_ext(),
			h(0.033),
			ym(1.0),
			b(),
			dim(2)
		{
		};
	};

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
		typename DerivedF,
		typename Derivedb,
		typename Derivedbc >
		IGL_INLINE bool xpbd_solve(
			Eigen::PlainObjectBase<DerivedV> & V,
			const Eigen::PlainObjectBase<DerivedF> & F,
			const int dim,
			const Eigen::PlainObjectBase<Derivedb> & b,
			const Eigen::PlainObjectBase<Derivedbc> & bc
			);

} // namespace XPBD

#ifndef IGL_STATIC_LIBRARY
#include "xpbd.cpp"
#endif