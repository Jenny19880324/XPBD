#include <igl/opengl/glfw/Viewer.h>
#include "meshgrid.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F, E;

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

	return false;
}

int main(int argc, char *argv[])
{
	const int rows = 16;
	const int cols = 24;
	Eigen::VectorXi x = Eigen::VectorXi::LinSpaced(cols, 0, cols - 1);
	Eigen::VectorXi y = Eigen::VectorXi::LinSpaced(rows, 0, rows - 1);


	igl::meshgrid(x, y, V, F, E);

	// Plot the mesh
	igl::opengl::glfw::Viewer viewer;
	viewer.callback_key_down = &key_down;
	viewer.data().set_mesh(V, F);
	viewer.data().set_face_based(true);
	viewer.launch();
}
