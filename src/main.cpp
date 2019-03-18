#include <igl/opengl/glfw/Viewer.h>
#include "meshgrid.h"
#include "xpbd.h"
#include <unordered_set>

Eigen::MatrixXd V, U;
Eigen::MatrixXi F, E;
xpbd::XPBDData xpbd_data;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
	std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
	if ((key == 'e') || (key == 'E')) {
		static bool show_edges = false;

		show_edges = !show_edges;

		if (show_edges) {
			viewer.data().set_edges(V, E, Eigen::RowVector3d(1., 0., 0.));
		}
		else {
			viewer.data().set_edges(Eigen::MatrixXd(), Eigen::MatrixXi(), Eigen::RowVector3d(1., 0., 0.));
		}
	}

	if (key == ' ') {
		viewer.core.is_animating = !viewer.core.is_animating;
	}

	return false;
}

int main(int argc, char *argv[])
{
	const int rows = 8;
	const int cols = 12;
	Eigen::VectorXi x = Eigen::VectorXi::LinSpaced(cols, 0, cols - 1);
	Eigen::VectorXi y = Eigen::VectorXi::LinSpaced(rows, 0, rows - 1);


	igl::meshgrid(x, y, V, F, E);

	U = V;

	std::unordered_set<int> b = { (rows - 1) * cols, rows * cols - 1 };
	if (!xpbd::xpbd_precomputation(V, F, E, 3, b, xpbd_data)) {
		std::cerr << "xpbd_precomputation failed." << std::endl;
		return EXIT_FAILURE;
	}

	const int n = V.rows();
	xpbd_data.f_ext = xpbd_data.M * Eigen::RowVector3d(0, -9.8, 0).replicate(n, 1);
	for (int i = 0; i < n; i++) {
		xpbd_data.vel.row(i) = V(i, 1) * Eigen::RowVector3d(0., 0., 2.);
	}

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.callback_key_down = &key_down;
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);

	viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer & viewer) ->bool
	{
		if (viewer.core.is_animating)
		{
			xpbd::xpbd_solve(U, E, xpbd_data);
			viewer.data().set_vertices(U);
			viewer.data().compute_normals();
		}
		return false;
	};

	viewer.launch();
}
