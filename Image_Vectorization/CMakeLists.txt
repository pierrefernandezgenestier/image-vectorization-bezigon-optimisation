
cmake_minimum_required(VERSION 3.1)
PROJECT(code_Image_Vectorization)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Enable C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

#OpenCV
SET(OpenCV_STATIC OFF)
FIND_PACKAGE(OpenCV REQUIRED QUIET)
INCLUDE_DIRECTORIES(${OpenCV_INCLUDE_DIRS})

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
find_package(LIBIGL REQUIRED QUIET)

file(GLOB SRCFILES *.cpp)
file(GLOB HEADERFILES *.h)
add_executable(image_vect ${SRCFILES} ${HEADERFILES})

target_link_libraries(image_vect igl::core igl::opengl_glfw)
TARGET_LINK_LIBRARIES(image_vect ${OpenCV_LIBS})






