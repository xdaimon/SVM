This code is for a paper I wrote in my linear algebra class. The program demonstrates the gradient descent method for the task of classification. While writing the code, I secretly joked with myself about there being more linear algebra involved in implementing the graphics than in coding the maths for the paper.

The program depends on the following libraries

```
glfw
GLEW
GLU
GL
pthread
```

You will need CMake installed to build the code.
To build and run just the following commands

```
mkdir build
cd build/
cmake ..
make
./main
```

![Here is a preview of the app](demo.gif)
