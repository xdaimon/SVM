#include "pch.h"

#include "src/ClassifierAdapter.h"
#include "src/MutliclassSVMAdapter.h"
#include "src/InputHandler.h"
#include "src/RenderingResources.h"
#include "Logger.h"

#include "draw.h"

static struct mouse m = {0.f, 0.f, -1.f, -1.f};
static double step_size = 0.001;
static double regularization = 0.001;
static int step_itters = 1;
static bool animate = false;

static float WindowWidth = INITIAL_WIDTH;
static float WindowHeight = INITIAL_HEIGHT;

static int currentClass;

using namespace SVM;

static RenderState state = RenderState::Render2D;
static ClassifierAdapter* mClassifierAdapter;
static InputHandler* mInputHandler;

static double start_time = 0.0;
static double getTime() {
	struct timespec res;
	clock_gettime(CLOCK_MONOTONIC, &res);
	return res.tv_sec + ((double)res.tv_nsec) / 1e9f - start_time;
}

int main() {

	// Setup timing stuff
	double current_time = 0.0;
	double last_time = 0.0;
	struct timespec res;
	clock_gettime(CLOCK_MONOTONIC, &res);
	start_time = res.tv_sec + ((float)res.tv_nsec) / 1e9f;

	int frame_count = 0;

	mClassifierAdapter = new MSVMAdapter();
	mClassifierAdapter->Setup(state);
	mInputHandler = new InputHandler(mClassifierAdapter);

	if (!initialize_gl(mClassifierAdapter->GetRenderingResources()))
		return 0;

	glfwSetTime(0.0);
	while (1) {
		if (!should_loop())
			break;

		if (animate) {
			mClassifierAdapter->TrainStep(step_size, step_itters, regularization);
		}

		mInputHandler->SetPointerCoordinateScale(WindowWidth, WindowHeight);
		draw(mClassifierAdapter->GetRenderingResources(), WindowWidth, WindowHeight);

		clock_gettime(CLOCK_MONOTONIC, &res);
		current_time = res.tv_sec + ((float)res.tv_nsec) / 1e9f - start_time;
		if (current_time - last_time > 1) {
			last_time = current_time;
			// printf("fps == %d\n", frame_count);
			frame_count = 0;
			LOG("reg_stgth: \t %f", regularization);
			LOG("step_size: \t %f", step_size);
			LOG("step_itters: \t %d", step_itters);
			LOG("\n");
		}
		frame_count++;
		glfwPollEvents();
	}
	deinit_drawing();
	glfwTerminate();
	return 0;
}

void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
	(void)window;
	m.mousey = ypos;
	m.mousex = xpos;
	if (m.mouse1_down > 0.0f)
		mInputHandler->OnPointerMove(m.mousex, m.mousey, getTime());
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
	(void)window;
	(void)mods;

	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		m.mouse1_down = 1.f;
		mInputHandler->OnPointerDown(m.mousex, m.mousey, getTime());
	}

	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
		m.mouse1_down = -1.f;
		mInputHandler->OnPointerUp(m.mousex, m.mousey, getTime(), currentClass + 1, state == RenderState::Render3D);
	}

	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
		m.mouse2_down = 1.f;

	if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_RELEASE)
		m.mouse2_down = -1.f;
}

static int currentDataSet = 0;
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	(void)scancode;
	(void)mods;

	if (action == GLFW_PRESS) {
		if (key == GLFW_KEY_ESCAPE)
			glfwSetWindowShouldClose(window, GL_TRUE);

		if (key == GLFW_KEY_C) {
			LOG("KEY C");
			currentClass++;
			currentClass %= SVM::MaxNumClasses;
			LOG("%d", currentClass);
		}

		if (key == GLFW_KEY_D) {
			LOG("KEY D");
			RenderState s;
			switch (state) {
			case RenderState::Render2D:
				s = RenderState::Render3D;
				break;
			case RenderState::Render3D:
				s = RenderState::Render2D;
				break;
			}

			mClassifierAdapter->Clear();
			mInputHandler->Clear();

			mClassifierAdapter->Setup(s);

			currentDataSet = 0;
			mClassifierAdapter->SetDataSet(currentDataSet);

			RenderingResources& rr = mClassifierAdapter->GetRenderingResources();

			if (!initialize_gl(rr)) {
				LOG("failed to initialize gl on RenderState switch");
				exit(-1);
			}

			state = s;
		}

		if (key == GLFW_KEY_W)
			mClassifierAdapter->PrintInformation();

		// if (key == GLFW_KEY_F)
		// 	mClassifierAdapter->WriteInformation();

		if (key == GLFW_KEY_T)
			mClassifierAdapter->Train();

		if (key == GLFW_KEY_K) {
			regularization *= 2.;
			if (regularization > 1000.)
				regularization = 1000.;
		}

		if (key == GLFW_KEY_J) {
			regularization /= 2.;
			if (regularization > 1000.)
				regularization = 1000.;
		}

		if (key == GLFW_KEY_UP) {
			step_size *= 2.;
			if (step_size > 1000.)
				step_size = 1000.;
		}

		if (key == GLFW_KEY_DOWN) {
			step_size /= 2.;
			if (step_size < 0.000001)
				step_size = 0.000001;
		}

		if (key == GLFW_KEY_RIGHT) {
			step_itters *= 2;
			if (step_itters > 2048)
				step_itters = 2048;
		}

		if (key == GLFW_KEY_LEFT) {
			step_itters /= 2;
			if (step_itters == 0)
				step_itters = 1;
		}

		if (key == GLFW_KEY_S)
			mClassifierAdapter->TrainStep(step_size, step_itters, regularization);

		if (key == GLFW_KEY_X)
			mClassifierAdapter->Clear();

		if (key == GLFW_KEY_R) {
			mClassifierAdapter->Clear();
			mClassifierAdapter->Setup(state);

			// hack
			mInputHandler->OnPointerDown(m.mousex, m.mousey, getTime());
			mInputHandler->OnPointerMove(m.mousex, m.mousey, getTime());

			mClassifierAdapter->Update();
			RenderingResources& rr = mClassifierAdapter->GetRenderingResources();
			for (size_t i = 0; i < rr.data.size(); ++i)
				rr.data[i].initializer(rr.data[i]);
			if (!initialize_gl(rr))
				exit(-1);
		}

		if (key == GLFW_KEY_Z) {
			mClassifierAdapter->Clear();

			currentDataSet++;
			currentDataSet %= SVM::NumDataSets;
			mClassifierAdapter->SetDataSet(currentDataSet);

			mClassifierAdapter->Setup(state);

			// hack
			mInputHandler->OnPointerDown(m.mousex, m.mousey, getTime());
			mInputHandler->OnPointerMove(m.mousex, m.mousey, getTime());

			mClassifierAdapter->Update();
			RenderingResources& rr = mClassifierAdapter->GetRenderingResources();
			for (size_t i = 0; i < rr.data.size(); ++i)
				rr.data[i].initializer(rr.data[i]);
			if (!initialize_gl(rr))
				exit(-1);
		}

		if (key == GLFW_KEY_P) {
			step_itters = 1;
			// step_size = 0.001;
			animate = !animate;
		}

		// if (key == GLFW_KEY_SPACE && action == GLFW_PRESS)
	}
}

void window_size_callback(GLFWwindow* window, int width, int height) {
	(void)window;
	WindowWidth = (float)width;
	WindowHeight = (float)height;
	glViewport(0, 0, width, height);
}
