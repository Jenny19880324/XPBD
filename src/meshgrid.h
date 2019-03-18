#pragma once

#include <igl/igl_inline.h>
#include <Eigen/Dense>

namespace igl
{
	// like matlab's meshgrid function
	// meshgrid(x, y, V, F) V is the #x by dim grid coordinates
	// usage:
	// x = Eigen::VectorXi::LinSpaced(n, 0, n-1);
	// y = Eigen::vectorXi::LinSpaced(n, 0, n-1);
	// meshgrid(x, y, V, F, E);
	// 
	//  ----------
	//  |\/|\/|\/|
	//  |/\|/\|/\|
	//  ----------
	//  |\/|\/|\/|
	//  |/\|/\|/\|
	//  ----------
	//  E contains 2 crossing springs within one grid
	//
	// Inputs:
	//   x  x-coordinates of points, specified as a 1D vector
	//   y  y-coordinates of points, specified as a 1D vector
	// Outputs:
	//   V #V by dim mesh positions
	//   F #F by simplex size indices into V, for visualization
	//   E #E by simplex size indices into V, for simulation
	template <
		typename Derivedx,
		typename DerivedV,
		typename DerivedF,
		typename DerivedE>
		IGL_INLINE void meshgrid(
			const Eigen::PlainObjectBase<Derivedx> & x,
			const Eigen::PlainObjectBase<Derivedx> & y,
			Eigen::PlainObjectBase<DerivedV> & V,
			Eigen::PlainObjectBase<DerivedF> & F,
			Eigen::PlainObjectBase<DerivedE> & E);
} // namespace igl

#ifndef IGL_STATIC_LIBRARY
#	include "meshgrid.cpp"
#endif
