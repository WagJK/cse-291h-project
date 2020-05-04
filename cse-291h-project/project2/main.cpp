#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <cstdlib>
#include <cmath>
#include "shader.h"
#include "camera.h"
#include "input.h"
#include "sph.h"
#include "sphintegrator.h"

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
Camera camera(glm::vec3(0.0f, 3.0f, 23.0f));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;	// time between current frame and last frame
float lastFrame = 0.0f;


void display(Shader myShader, float* vertices, int len_vertices) {

    // activate shader
    myShader.use();

    // pass projection matrix to shader (note that in this case it could change every frame)
    glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    myShader.setMat4("projection", projection);

    // camera/view transformation
    glm::mat4 view = camera.GetViewMatrix();
    myShader.setMat4("view", view);

    // calculate the model matrix for each object and pass it to shader before drawing
    glm::mat4 model = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    myShader.setMat4("model", model);

    glm::vec3 fragcolor(1.0f, 1.0f, 1.0f);
    myShader.setVec3("fragColor", fragcolor);


    // render boxes
    unsigned int VBO, VAO, EBO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, len_vertices * sizeof(float), vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, len_indices * sizeof(unsigned int), indices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);


    glBindVertexArray(VAO);
    glDrawElements(GL_LINES, len_indices, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
}


void display_acc(Shader myShader, float* accs, int len_accs) {
    // activate shader
    myShader.use();

    // pass projection matrix to shader (note that in this case it could change every frame)
    glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    myShader.setMat4("projection", projection);

    // camera/view transformation
    glm::mat4 view = camera.GetViewMatrix();
    myShader.setMat4("view", view);

    // calculate the model matrix for each object and pass it to shader before drawing
    glm::mat4 model = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    myShader.setMat4("model", model);

    glm::vec3 fragcolor(1.0f, 1.0f, 0.0f);
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

    // build and compile our shader zprogram
    // ------------------------------------
    Shader myShader("project2/shader/vertex.shader", "project2/shader/fragment.shader");

    // set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------    
    const float k               = 1.0;
    const float density0        = 1.0;
    const float supportRadius   = 10.0;
    const float smoothingRadius = 10.0;
    const float penalty         = 200;
    const vec3 pos(0, 0, 0), size(10, 10, 10), gap(1, 1, 1), vel(0, 0, 0);
    SPHSystem sph(pos, size, gap, k, density0, supportRadius, smoothingRadius);
    SPHIntegrator itg(penalty);
    
    float* vertices = new float[3 * sph.getSize()];
    float* accs     = new float[6 * sph.getSize()];
    unsigned int* indices = new unsigned int[2 * sph.getNumLines()];

    sph.getPositions(vertices); 
    sph.applyGForces();
    
    display(myShader, vertices, 3 * sph.getSize());

    // render loop
    // -----------
    int cnt = 0;
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        // --------------------
        float currentFrame = glfwGetTime();
        // deltaTime = currentFrame- lastFrame;
        deltaTime = 0.01f;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        sph.clearTempForces();
        sph.applySPHForces();

        if (!pause) itg.Integrate(sph, deltaTime);

        sph.getPositions(vertices);
        // sph.getPosAcc(accs);
        // sph.PrintPositions();
        
        if (cnt % 5 == 0) {
            display(myShader, vertices, 3 * sph.getSize());
            //display_acc(myShader, accs, 6 * sph.getSize());
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


