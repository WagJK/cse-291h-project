#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <cstdlib>
#include <cmath>
#include "shader.h"
#include "camera.h"
#include "sph.h"
#include "sphintegrator.h"

using namespace std;
using namespace glm;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);

void mouse_callback(GLFWwindow* window, double xpos, double ypos);

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);

void processInput(GLFWwindow* window);

void checkGLerrors() {
    GLenum err;
    while ((err = glGetError()) != GL_NO_ERROR)
        std::cout << err;
}

// settings
const unsigned int SCR_WIDTH = 1920;
const unsigned int SCR_HEIGHT = 1080;

// camera
Camera camera(vec3(0.0f, 3.0f, 23.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;	// time between current frame and last frame
float totalTime = 0.0f;
float lastFrame = 0.0f;

void display_sph(Shader myShader, float* vertices, int len_vertices) {
    // activate shader
    myShader.use();

    // pass projection matrix to shader (note that in this case it could change every frame)
    mat4 projection = perspective(radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    mat4 view = camera.GetViewMatrix();
    mat4 model = mat4(1.0f);
    vec3 fragcolor(1.0f, 1.0f, 1.0f);

    myShader.setMat4("projection", projection);
    myShader.setMat4("model", model);
    myShader.setMat4("view", view);
    myShader.setVec3("fragColor", fragcolor);
    myShader.setFloat("pointScale", 1000.0);
    myShader.setFloat("pointRadius", 1.0);

    // render boxes
    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, len_vertices * sizeof(float), vertices, GL_STATIC_DRAW);

    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    // glBufferData(GL_ELEMENT_ARRAY_BUFFER, len_indices * sizeof(unsigned int), indices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(VAO);
    glDrawArrays(GL_POINTS, 0, len_vertices);
    glBindVertexArray(0);

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
}

void display_acc(Shader myShader, float* accs, int len_accs) {
    // activate shader
    myShader.use();

    // pass projection matrix to shader (note that in this case it could change every frame)
    mat4 projection = perspective(radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    myShader.setMat4("projection", projection);

    // camera/view transformation
    mat4 view = camera.GetViewMatrix();
    myShader.setMat4("view", view);

    // calculate the model matrix for each object and pass it to shader before drawing
    mat4 model = mat4(1.0f); // make sure to initialize matrix to identity matrix first
    myShader.setMat4("model", model);

    vec3 fragcolor(1.0f, 1.0f, 0.0f);
    myShader.setVec3("fragColor", fragcolor);

    // render boxes
    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, len_accs * sizeof(float), accs, GL_STATIC_DRAW);

    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    // glBufferData(GL_ELEMENT_ARRAY_BUFFER, len_indices * sizeof(unsigned int), indices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(VAO);
    glDrawArrays(GL_LINES, 0, len_accs);
    // glDrawElements(GL_LINES, len_indices, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
}

void display_ctn(Shader myShader, vec3 lb, vec3 ub) {
    vec3 half = {0.5, 0.5, 0.5};
    lb -= half;
    ub += half;
    float vertices[24] = {
        lb.x, lb.y, lb.z,
        lb.x, lb.y, ub.z,
        lb.x, ub.y, lb.z,
        lb.x, ub.y, ub.z,
        ub.x, lb.y, lb.z,
        ub.x, lb.y, ub.z,
        ub.x, ub.y, lb.z,
        ub.x, ub.y, ub.z,
    };

    unsigned int indices[24] = {
        0, 1, 0, 2, 0, 4, 1, 3, 1, 5, 2, 3,
        2, 6, 3, 7, 4, 5, 4, 6, 5, 7, 6, 7,
    };

    // activate shader
    myShader.use();

    // pass projection matrix to shader (note that in this case it could change every frame)
    mat4 projection = perspective(radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    myShader.setMat4("projection", projection);

    // camera/view transformation
    mat4 view = camera.GetViewMatrix();
    myShader.setMat4("view", view);

    // calculate the model matrix for each object and pass it to shader before drawing
    mat4 model = mat4(1.0f); // make sure to initialize matrix to identity matrix first
    myShader.setMat4("model", model);

    vec3 fragcolor(1.0f, 1.0f, 0.0f);
    myShader.setVec3("fragColor", fragcolor);

    // render boxes
    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, 24 * sizeof(float), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 24 * sizeof(unsigned int), indices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    glBindVertexArray(VAO);
    glDrawElements(GL_LINES, 24, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
}

bool pause = true;
int main()
{
    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this statement to fix compilation on OS X
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);
    glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);

    // build and compile our shader zprogram
    // ------------------------------------
    Shader sphereShader ("project2/shader/sphere_vertex.shader",    "project2/shader/sphere_fragment.shader");
    Shader lineShader   ("project2/shader/line_vertex.shader",      "project2/shader/line_fragment.shader");

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------    
    const float k               = 5e-4;
    const float density0        = 1e3;
    const float supportRadius   = 6.0;
    const float smoothingRadius = 3.0;
    const float penalty         = 200000000;
    const vec3  pos(0, 0, 0);
    const vec3  size(4, 8, 4);
    const vec3  gap(1, 0.2, 1);
    const vec3  m_d(supportRadius, supportRadius, supportRadius);
    const vec3  container_lb(-0.1, 0, -0.1);
    const vec3  container_ub(3.1, 7.1, 3.1);

    SPHSystem sph(pos, size, gap, m_d, container_lb, container_ub, penalty, k, density0, supportRadius, smoothingRadius);
    SPHIntegrator itg(penalty);
    
    float* vertices = new float[3 * sph.getSize()];
    float* accs     = new float[6 * sph.getSize()];

    sph.getPositions(vertices); 
    sph.applyGForces();
    display_sph(sphereShader, vertices, 3 * sph.getSize());

    // render loop
    // -----------
    int cnt = 0;
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        // --------------------
        float currentFrame = glfwGetTime();
        // deltaTime = currentFrame- lastFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (!pause) {
            deltaTime = 0.01f;
            totalTime += deltaTime;
            lastFrame = currentFrame;

            sph.build();
            sph.clearTempForces();
            sph.applySPHForces();

            vec3 shrink(1, 0, 0);
            shrink = shrink * (sin(totalTime * 3 - PI / 2) + 1) * 1 / 2;
            sph.setContainer(container_lb, container_ub - shrink);
            sph.applyPenaltyForces();

            itg.Integrate(sph, deltaTime);

            sph.getPositions(vertices);
            sph.getPosAcc(accs);
            
            if (cnt % 40 == 0) {
                display_ctn(lineShader, container_lb, container_ub - shrink);
                display_sph(sphereShader, vertices, 3 * sph.getSize());
                // display_acc(myShader, accs,  6 * sph.getSize());
            }
        }
        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
        checkGLerrors();
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS)
        pause ^= true;
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }
    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top
    lastX = xpos;
    lastY = ypos;
    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}