#include "meshgrid.h"
#include <igl/LinSpaced.h>

template <
	typename Derivedx,
	typename DerivedV,
	typename DerivedF,
	typename DerivedE>
	IGL_INLINE void igl::meshgrid(
		const Eigen::PlainObjectBase<Derivedx> & x,
		const Eigen::PlainObjectBase<Derivedx> & y,
		Eigen::PlainObjectBase<DerivedV> & V,
		Eigen::PlainObjectBase<DerivedF> & F,
		Eigen::PlainObjectBase<DerivedE> & E
	)
{
	const int m = y.rows();
	const int n = x.rows();
	
	V.resize(m * n, 3);
	F.resize((m - 1)*(n-1) * 2, 3);
	//E.resize((m - 1)*n + (n - 1)*m, 2);
	E.resize((m - 1)*n + (n - 1)*m + 2 * (m - 1)*(n - 1), 2);

	// generate the vertices
	for (int x_i = 0; x_i < n; x_i++) {
		for (int y_i = 0; y_i < m; y_i++) {
			V(n * y_i + x_i, 0) = x(x_i);
			V(n * y_i + x_i, 1) = y(y_i);
			V(n * y_i + x_i, 2) = 0.;
		}
	}

	// generate the triangles for visualization
	for (int x_i = 0; x_i < n - 1; x_i++) {
		for (int y_i = 0; y_i < m - 1; y_i++) {
			F(2 * ((n - 1) * y_i + x_i) + 0, 0) = n * y_i + x_i;
			F(2 * ((n - 1) * y_i + x_i) + 0, 1) = n * y_i + x_i + 1;
			F(2 * ((n - 1) * y_i + x_i) + 0, 2) = n * (y_i + 1) + x_i + 1;
			F(2 * ((n - 1) * y_i + x_i) + 1, 0) = n * y_i + x_i;
			F(2 * ((n - 1) * y_i + x_i) + 1, 1) = n * (y_i + 1) + x_i + 1;
			F(2 * ((n - 1) * y_i + x_i) + 1, 2) = n * (y_i + 1) + x_i;
		}
	}

	// generate the edges for simulation
	// generate all the horizontal edges
	int edge_idx = 0;
	for (int x_i = 0; x_i < n - 1; x_i++) {
		for (int y_i = 0; y_i < m; y_i++) {
			E(edge_idx, 0) = n * y_i + x_i;
			E(edge_idx, 1) = n * y_i + x_i + 1;
			edge_idx++;
		}
	}

	// generate all the vertical edges
	for (int x_i = 0; x_i < n; x_i++) {
		for (int y_i = 0; y_i < m - 1; y_i++) {
			E(edge_idx, 0) = n * y_i + x_i;
			E(edge_idx, 1) = n * (y_i + 1) + x_i;
			edge_idx++;
		}
	}

	// generate all the crossing edges
	for (int x_i = 0; x_i < n - 1; x_i++) {
		for (int y_i = 0; y_i < m - 1; y_i++) {
			E(edge_idx, 0) = n * y_i + x_i;
			E(edge_idx, 1) = n * (y_i + 1) + x_i + 1;
			edge_idx++;
			E(edge_idx, 0) = n * y_i + x_i + 1;
			E(edge_idx, 1) = n * (y_i + 1) + x_i;
			edge_idx++;
		}
	}
}

#ifdef IGL_STATIC_LIBRARY
template void igl::meshgrid<class Eigen::Matrix<int, -1, 1, 0, -1, 1>, class Eigen::Matrix<double, -1, -1, 0, -1, -1>, class Eigen::Matrix<int, -1, -1, 0, -1, -1>, class Eigen::Matrix<int, -1, -1, 0, -1, -1> >(class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, 1, 0, -1, 1> > const &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, 1, 0, -1, 1> > const &, class Eigen::PlainObjectBase<class Eigen::Matrix<double, -1, -1, 0, -1, -1> > &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1> > &, class Eigen::PlainObjectBase<class Eigen::Matrix<int, -1, -1, 0, -1, -1> > &);
#endif