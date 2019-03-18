#include "xpbd.h"
#include <igl/massmatrix.h>
#include <unordered_map>


template <
	typename DerivedV,
	typename DerivedF,
	typename DerivedE,
	typename Derivedb>
	IGL_INLINE bool xpbd::xpbd_precomputation(
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F,
		const Eigen::PlainObjectBase<DerivedE> & E,
		const int dim,
		const Eigen::PlainObjectBase<Derivedb> & b,
		XPBDData & data
	)
{
	using namespace std;
	using namespace Eigen;

	data.V = V;
	// number of vertices
	const int n = V.rows();
	data.n = n;
	assert((b.size() == 0 || b.maxCoeff() < n) && "b out of bounds");
	assert((b.size() == 0 || b.minCoeff() >= 0) && "b out of bounds");
	// remember b
	data.b = b;
	assert((dim == 3 || dim == 2) && "dim should be 2 or 3");
	data.dim = dim;
	data.f_ext = MatrixXd::Zero(n, data.dim);
	data.vel = MatrixXd::Zero(n, data.dim);
	data.lambda_iter = VectorXd::Zero(E.rows());

	SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, data.M);
	data.M_inv = data.M.cwiseInverse();

	return true;
}


template <
	typename DerivedV,
	typename DerivedF,
	typename Derivedbc>
	IGL_INLINE bool xpbd::xpbd_solve(
		Eigen::PlainObjectBase<DerivedV> & U,
		const Eigen::PlainObjectBase<DerivedF> & E,
		const Eigen::PlainObjectBase<Derivedbc> & bc,
		XPBDData & data)
{
	using namespace Eigen;

	static std::vector<double> d(E.rows()); // rest_length
	const Eigen::MatrixXd V = data.V;
	for (int j = 0; j < E.rows(); j++) {
		d[j] = (V.row(E(j, 0)) - V.row(E(j, 1))).norm();
	}


	auto C = [&](int j) {
		return (U.row(E(j, 0)) - U.row(E(j, 1))).norm() - d[j];
	};

	auto dC = [&](int j) -> Eigen::RowVector3d {
		Eigen::RowVector3d l = U.row(E(j, 0)) - U.row(E(j, 1));
		return (l.norm() - d[j]) * l.normalized();
	};

	auto dC_norm = [&](int j) {
		Eigen::RowVector3d l = U.row(E(j, 0)) - U.row(E(j, 1));
		return l.norm() - d[j];
	};

	// predict position
	const double h = data.h;
	MatrixXd U_prev = U;
	MatrixXd U0 = U + h * data.vel + h * h * data.M_inv * data.f_ext;

	int iter = 0;
	double tilde_alpha = data.ym / (h * h);
	while (iter < data.max_iter) {

		for (int j = 0; j < E.rows(); j++) {

			// compute dlambda using Eq(18)
			double Cj = C(j);
			double dCj_norm = dC_norm(j);
			Eigen::RowVector3d dCj = dC(j);
			double m_inv = data.M_inv.coeff(E(j, 0), E(j, 0));
			double dlambda_j = (-Cj - tilde_alpha * data.lambda_iter(j)) / (m_inv * dCj_norm + tilde_alpha);

			// compute dx using Eq(17)
			Eigen::RowVector3d dxj = m_inv * dCj * dlambda_j;

			// update lambda_{i+1} <= lambda_{i} + \nabla lambda
			data.lambda_iter(j) += dlambda_j;

			// update x_{i+1} <= x_{i} + \nabla x
			U0.row(E(j, 0)) += dxj;
			U0.row(E(j, 1)) -= dxj;
		}

		iter++;
	}

	U = U0;
	data.vel = (U - U_prev) / h;

	return true;
}


#ifdef IGL_STATIC_LIBRARY
template bool xpbd::xpbd_precomputation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, 1, 0, -1, 1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, int, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, 1, 0, -1, 1> > const &, struct xpbd::XPBDData &);
template bool xpbd::xpbd_solve<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<double, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, struct xpbd::XPBDData &);
#endif