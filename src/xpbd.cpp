#include "xpbd.h"
#include <igl/massmatrix.h>
#include <unordered_map>


template <
	typename DerivedV,
	typename DerivedF,
	typename DerivedE>
	IGL_INLINE bool xpbd::xpbd_precomputation(
		const Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F,
		const Eigen::PlainObjectBase<DerivedE> & E,
		const int dim,
		const std::unordered_set<int> & b,
		XPBDData & data)
{
	using namespace std;
	using namespace Eigen;

	data.V = V;
	// number of vertices
	const int n = V.rows();
	data.n = n;
	// remember b
	data.b = b;
	assert((dim == 3 || dim == 2) && "dim should be 2 or 3");
	data.dim = dim;
	data.f_ext = MatrixXd::Zero(n, data.dim);
	data.vel = MatrixXd::Zero(n, data.dim);
	data.lambda_iter = VectorXd::Zero(E.rows());

	SparseMatrix<double> M;
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, data.M);
	data.M = 1e-5 * data.M;
	data.M_inv = data.M.cwiseInverse();

	return true;
}


template <
	typename DerivedV,
	typename DerivedF>
	IGL_INLINE bool xpbd::xpbd_solve(
		Eigen::PlainObjectBase<DerivedV> & U,
		const Eigen::PlainObjectBase<DerivedF> & E,
		XPBDData & data)
{
	using namespace Eigen;

	// predict position
	const double h = data.h;
	const int n = U.rows();
	Eigen::MatrixXd U_prev = U;
	Eigen::MatrixXd U0 = U + h * data.vel + h * h * data.M_inv * data.f_ext;

	// The Position of the vertex is simply set to the static target position
	// or updated at every time step to coincide with the position of the kinematic object.
	for (auto it = data.b.begin(); it != data.b.end(); it++) {
		U0.row(*it) = U.row(*it);
	}


	static std::vector<double> d(E.rows()); // rest_length
	const Eigen::MatrixXd V = data.V;
	for (int j = 0; j < E.rows(); j++) {
		d[j] = (V.row(E(j, 0)) - V.row(E(j, 1))).norm();
	}


	auto C = [&](int j) {
		return (U0.row(E(j, 0)) - U0.row(E(j, 1))).norm() - d[j];
	};

	auto dC = [&](int j) -> Eigen::RowVector3d {
		Eigen::RowVector3d l = U0.row(E(j, 0)) - U0.row(E(j, 1));
		return l.normalized();
	};

	auto dC_squaredNorm = [&](int j) {
		Eigen::RowVector3d l = U0.row(E(j, 0)) - U0.row(E(j, 1));
		double norm = l.norm() - d[j];
		//return norm * norm;
		return 1.;
	};

	int iter = 0;
	double tilde_alpha = (1./ data.ym) / (h * h);
	while (iter < data.max_iter) {

		for (int j = 0; j < E.rows(); j++) {

			// compute dlambda using Eq(18)
			double Cj = C(j);
			double dCj_squaredNorm = dC_squaredNorm(j);
			Eigen::RowVector3d dCj = dC(j);
			double m_inv_0 = data.M_inv.coeff(E(j, 0), E(j, 0));
			double m_inv_1 = data.M_inv.coeff(E(j, 1), E(j, 1));
			double m_inv_sum = m_inv_0 + m_inv_1;
			double dlambda_j = (-Cj - tilde_alpha * data.lambda_iter(j)) / (m_inv_sum + tilde_alpha);

			// compute dx using Eq(17)
			Eigen::RowVector3d dxj = dCj * dlambda_j;

			// update lambda_{i+1} <= lambda_{i} + \nabla lambda
			data.lambda_iter(j) += dlambda_j;

			// update x_{i+1} <= x_{i} + \nabla x
			if (data.b.find(E(j, 0)) == data.b.end()) {
				U0.row(E(j, 0)) += m_inv_0 * dxj;
			}
			if (data.b.find(E(j, 1)) == data.b.end()) {
				U0.row(E(j, 1)) -= m_inv_1 * dxj;
			}
		}

		iter++;
	}

	U = U0;
	data.vel = (U - U_prev) / h;

	return true;
}


#ifdef IGL_STATIC_LIBRARY
template bool xpbd::xpbd_precomputation<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, int, std::unordered_set<int, struct std::hash<int>, struct std::equal_to<int>, std::allocator<int> > const &, struct xpbd::XPBDData &);
template bool xpbd::xpbd_solve<Eigen::Matrix<double, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1> >(Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1> > &, Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const &, struct xpbd::XPBDData &);
#endif