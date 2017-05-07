#include "src/RenderingResources.h"

struct mouse
{
    float mousex;
    float mousey;
    float mouse1_down;
    float mouse2_down;
};

bool initialize_gl(SVM::RenderingResources & rs);
void deinit_drawing();

void draw(SVM::RenderingResources & rs, float width, float height);

bool should_loop();

void window_size_callback(GLFWwindow* window, int width, int height);
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);

#define INITIAL_WIDTH 1366
#define INITIAL_HEIGHT 768
