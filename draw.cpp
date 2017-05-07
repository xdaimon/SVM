#include "pch.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <vector>
#include <iostream>
using std::cout;
using std::endl;
#include <math.h>

#include <unistd.h>
#include <sys/time.h>

#include "src/RenderingResources.h"
#include "Logger.h"

#include "draw.h"

static GLFWwindow* window;

static GLuint vbo;
static GLuint vao;
static GLuint program;
// static GLint time_name;
static GLint mouse_name;
static GLint resolution_name;
static GLint position_attrib_location;

void log_gl_error() {
	const GLubyte* err = gluErrorString(glGetError());
	while (*err) {
		std::cout << *err;
		err++;
	}
	std::cout << std::endl;
}

static bool compile_shader(const char* s, GLuint& sn, GLenum stype) {
	sn = glCreateShader(stype);
	glShaderSource(sn, 1, &s, NULL);
	glCompileShader(sn);
	GLint isCompiled = 0;
	glGetShaderiv(sn, GL_COMPILE_STATUS, &isCompiled);
	if(isCompiled == GL_FALSE) {
		GLint maxLength = 0;
		glGetShaderiv(sn, GL_INFO_LOG_LENGTH, &maxLength);
		std::vector<GLchar> errorLog(maxLength);
		glGetShaderInfoLog(sn, maxLength, &maxLength, &errorLog[0]);
		for (auto c : errorLog)
			cout << c;
		cout << endl;
		glDeleteShader(sn);
		return false;
	}
	return true;
}
static bool compile_shaders(SVM::RenderingResources& rr) {

	if (program)
		glDeleteProgram(program);

	const char* vertex_shader = R"(
		#version 300 es
		in vec4 position;
		void main() {
			gl_Position = position;
		}
	)";
	const char* fragment_shader = rr.fragmentShader.data();
	GLuint vs, fs;
	if (!compile_shader(vertex_shader, vs, GL_VERTEX_SHADER)) return false;
	if (!compile_shader(fragment_shader, fs, GL_FRAGMENT_SHADER)) return false;
	program = glCreateProgram();
	glAttachShader(program, fs);
	glAttachShader(program, vs);
	glLinkProgram(program);
	GLint isLinked = 0;
	glGetProgramiv(program, GL_LINK_STATUS, (int*)&isLinked);
	if (isLinked == GL_FALSE) {
		GLint maxLength = 0;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);
		std::vector<GLchar> infoLog(maxLength);
		glGetProgramInfoLog(program, maxLength, &maxLength, &infoLog[0]);
		for (auto c : infoLog)
			std::cout << c;
		glDeleteProgram(program);
		glDeleteShader(vs);
		glDeleteShader(fs);
		return false;
	}

	glUseProgram(program);

	for (size_t i = 0; i < rr.data.size(); ++i)
		rr.data[i].initializer(rr.data[i]);

	// time_name = glGetUniformLocation(program, "time");
	//mouse_name = glGetUniformLocation(program, "mouse");
	position_attrib_location = glGetAttribLocation(program, "position");
	resolution_name = glGetUniformLocation(program, "resolution");

	return true;
}

void init_desktop() {
	glfwInit();

	window = glfwCreateWindow(INITIAL_WIDTH, INITIAL_HEIGHT, "MSVM", NULL, NULL);
	glfwMakeContextCurrent(window);
	glfwSetWindowSizeCallback(window, window_size_callback);
	glfwSwapInterval(1);

	glewExperimental = GL_TRUE;
	glewInit();

	// const GLubyte* renderer = glGetString(GL_RENDERER);
	// const GLubyte* version = glGetString(GL_VERSION);
	// const GLubyte* ext = glGetString(GL_EXTENSIONS);
	// printf ("Renderer: %s\n", renderer);
	// printf ("OpenGL version supported %s\n", version);
	// printf ("OpenGL extensions supported %s\n", ext);

	glViewport(0, 0, INITIAL_WIDTH, INITIAL_HEIGHT);

	glfwSetKeyCallback(window, key_callback);
	glfwSetCursorPosCallback(window, cursor_position_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);
}

static bool once = false;
bool initialize_gl(SVM::RenderingResources& rr) {
	if (!once)
		init_desktop();
	once = true;

	glEnable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);

	bool b = compile_shaders(rr);

	glDisable(GL_DEPTH_TEST);
	glClearColor(0.0, 0.0, 0.0, 1.0);

	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	// clang-format off
	GLfloat vertexPositions[] = {
		-1.0f,-1.0f, 0.0,
		-1.0f, 1.0f, 0.0,
		 1.0f, 1.0f, 0.0,
		 1.0f,-1.0f, 0.0
	};

	GLuint indices[] = {
		0,1,2,
		0,2,3
	}; // clang-format on

	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertexPositions), vertexPositions, GL_STATIC_DRAW);

	GLuint indexBuffer;
	glGenBuffers(1, &indexBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, indexBuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertexPositions), indices, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glVertexAttribPointer(
	    position_attrib_location, // attrib location in shader
	    3,                        // 3 coords per vertex
	    GL_FLOAT,                 // type of each vertex
	    GL_FALSE,                 // Normalize?
	    sizeof(GLfloat) * 3,      // stride in bytes
	    0);
	glEnableVertexAttribArray(position_attrib_location);

	glUseProgram(program);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBindVertexArray(vao);

	return b;
}

void deinit_drawing() {
	glfwDestroyWindow(window);
}

void draw(SVM::RenderingResources& rr, float width, float height) {
	glClear(GL_COLOR_BUFFER_BIT);
	glUniform2f(resolution_name, (GLfloat)width, (GLfloat)height);
	for (size_t i = 0; i < rr.data.size(); ++i)
		rr.data[i].uploader(rr.data[i]);
	glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
	glfwSwapBuffers(window);
}

bool should_loop() {
	return !glfwWindowShouldClose(window);
}
