cmake_minimum_required(VERSION 3.12)
project(XPBD)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library"  ON)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                    ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                      ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"             ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI      "Use ImGui"                     ON)
option(LIBIGL_WITH_PNG               "Use PNG"                       ON)

find_package(LIBIGL REQUIRED QUIET)

# Add your project files
file(GLOB SRCFILES 
src/*.cpp
src/*.h)

add_executable(${PROJECT_NAME} ${SRCFILES})
target_link_libraries(${PROJECT_NAME} 
	igl::core 
	igl::opengl
	igl::opengl_glfw
	igl::opengl_glfw_imgui
	igl::png)
target_include_directories(${PROJECT_NAME} PUBLIC ${LIBIGL_EXTERNAL}/stb)
