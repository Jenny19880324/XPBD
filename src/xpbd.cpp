#include "xpbd.h"

template <
	typename DerivedV,
	typename DerivedF,
	typename Derivedb,
	typename Derivedbc>
	IGL_INLINE bool xpbd::xpbd_solve(
		Eigen::PlainObjectBase<DerivedV> & V,
		const Eigen::PlainObjectBase<DerivedF> & F,
		const int dim,
		const Eigen::PlainObjectBase<Derivedb> & b,
		const Eigen::PlainObjectBase<Derivedbc> & bc)
{
	using namespace std;
	using namespace Eigen;

	// number of vertices
	const int n = V.rows();

	assert((b.size() == 0 || b.maxCoeff() < n) && "b out of bounds");
	assert((b.size() == 0 || b.minCoeff() >= 0) && "b out of bounds");

	assert((dim == 3 || dim == 2) && "dim should be 2 or 3");


}


#ifdef IGL_STATIC_LIBRARY
#endif