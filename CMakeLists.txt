cmake_minimum_required(VERSION 3.3)
project(main)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DINFOLOG -O3 -std=c++17")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG -O0 -g -std=c++11")

set(SOURCE_FILES
	pch.cpp
    draw.cpp
    draw.h
    main.cpp
	Logger.cpp
	src/Matrix.cpp
	src/Vector.cpp
	src/MulticlassSVM.cpp
	src/MulticlassSVMAdapter.cpp
	src/RenderingResources.cpp
	src/InputHandler.cpp)

include_directories(.)

add_executable(main ${SOURCE_FILES})

# LINK_DIRECTORIES(/usr/lib/x86_64-linux-gnu/)
TARGET_LINK_LIBRARIES(main glfw GLEW GLU GL pthread)
# TARGET_LINK_LIBRARIES(main glfw GLEW GLU GL pthread png12)

