#include <Windows.h>
#include <iostream>
#include <limits>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

#include <GL/glew.h>
#include <GL/GL.h>
#include <GL/freeglut.h>

#define GLFW_INCLUDE_GLU
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <vector>

#define GLM_SWIZZLE
#include <glm/glm.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/string_cast.hpp>

using namespace glm;

// -------------------------------------------------
// Global Variables
// -------------------------------------------------
int Width = 512;
int Height = 512;
std::vector<float> OutputImage;
// -------------------------------------------------
//resterization
//--------------------------------------------------
int gNumVertices = 0;    // Number of 3D vertices.
int gNumTriangles = 0;    // Number of triangles.
int* gIndexBuffer = NULL; // Vertex indices for the triangles.
vec3* gVertexBuffer = NULL; // Vertex positions for the vertices.
// -------------------------------------------------

// -------------------------------------------------
// Classes
// -------------------------------------------------

class Surface {
public:
    vec3 color;
    Surface(const vec3& c) : color(c) {};
};

class Sphere : public Surface {
public:
    vec3 center;
    float radius;

    Sphere(const vec3& c, float r, const vec3& col)
        : Surface(col), center(c), radius(r) {
    }
};

//class Plane : public Surface {
//public:
//    vec3 normal;
//    float d;
//
//    Plane(const vec3& n, float d, const vec3& col)
//        : Surface(col), normal(normalize(n)), d(d) {
//    }
//};


class Camera {
public:
    vec3 position, direction;
    Camera(const vec3& pos, const vec3& dir)
        : position(pos), direction(normalize(dir)) {
    }
};

class Scene {
public:
    std::vector<Surface*> surfaces;
    Camera camera;

    Scene(const std::vector<Surface*>& sur, const Camera& c)
        : surfaces(sur), camera(c) {
    }
};

class Image {
public:
    void set(int x, int y, vec3 color) {
        OutputImage[(y * Width + x) * 3] = color.x;
        OutputImage[(y * Width + x) * 3 + 1] = color.y;
        OutputImage[(y * Width + x) * 3 + 2] = color.z;
    }
};

void create_scene()
{
    int width = 32;
    int height = 16;

    float theta, phi;
    int t;

    gNumVertices = (height - 2) * width + 2;
    gNumTriangles = (height - 2) * (width - 1) * 2;

    //TODO: Allocate an array for gNumVertices vertices.
    gVertexBuffer = new vec3[gNumVertices];

    gIndexBuffer = new int[3 * gNumTriangles];

    t = 0;
    for (int j = 1; j < height - 1; ++j)//1~15
    {
        for (int i = 0; i < width; ++i)//0~31
        {
            theta = (float)j / (height - 1) * M_PI;
            phi = (float)i / (width - 1) * M_PI * 2;

            float   x = sinf(theta) * cosf(phi);
            float   y = cosf(theta);
            float   z = -sinf(theta) * sinf(phi);

            //Set vertex t in the vertex array to {x, y, z}.
            gVertexBuffer[t] = vec3(x, y, z);

            t++;
        }
    }

    //Set vertex t in the vertex array to {0, 1, 0}.
    gVertexBuffer[t] = vec3(0, 1, 0);
    t++;

    //Set vertex t in the vertex array to {0, -1, 0}.
    gVertexBuffer[t] = vec3(0, -1, 0);
    t++;

    t = 0;
    for (int j = 0; j < height - 3; ++j)
    {
        for (int i = 0; i < width - 1; ++i)
        {
            gIndexBuffer[t++] = j * width + i;
            gIndexBuffer[t++] = (j + 1) * width + (i + 1);
            gIndexBuffer[t++] = j * width + (i + 1);
            gIndexBuffer[t++] = j * width + i;
            gIndexBuffer[t++] = (j + 1) * width + i;
            gIndexBuffer[t++] = (j + 1) * width + (i + 1);
        }
    }
    for (int i = 0; i < width - 1; ++i)
    {
        gIndexBuffer[t++] = (height - 2) * width;
        gIndexBuffer[t++] = i;
        gIndexBuffer[t++] = i + 1;
        gIndexBuffer[t++] = (height - 2) * width + 1;
        gIndexBuffer[t++] = (height - 3) * width + (i + 1);
        gIndexBuffer[t++] = (height - 3) * width + i;
    }
}

mat4 cal_modelTransform(Sphere sphere) {
    //1. Modeling Transform(scale → translate)
    mat4 S(1.0f), T(1.0f);
    S[0][0] = sphere.radius;
    S[1][1] = sphere.radius;
    S[2][2] = sphere.radius;
    T[3][0] = sphere.center.x;
    T[3][1] = sphere.center.y;
    T[3][2] = sphere.center.z;
    return T * S;
}

mat4 cal_mvp(Camera camera, mat4 Mm) {
    // 2. Camera Transform
    mat4 Mcam(1.0f);
    Mcam[3][0] = -1.0f * camera.position.x;
    Mcam[3][1] = -1.0f * camera.position.y;
    Mcam[3][2] = -1.0f * camera.position.z;

    // 3. Perspective Projection
    float l = -0.1f, r = 0.1f, b = -0.1f, t = 0.1f, n = -0.1f, f = -1000.0f;
    mat4 P(0.0f);
    P[0][0] = (2.0f * n) / (r - l);
    P[1][1] = (2.0f * n) / (t - b);
    P[2][0] = (l + r) / (l - r);
    P[2][1] = (b + t) / (b - t);
    P[2][2] = (f + n) / (n - f);
    P[3][2] = (2.0f * f * n) / (f - n);
    P[2][3] = 1.0f;

    // 4. Viewport Transform
    mat4 Vp(1.0f);
    Vp[0][0] = Width / 2.0f;
    Vp[1][1] = Height / 2.0f;
    Vp[3][0] = (Width - 1.0f) / 2.0f;
    Vp[3][1] = (Height - 1.0f) / 2.0f;

    // Final Matrix
    return Vp * P * Mcam * Mm;
}

//baycentric algorithm
void cal_baycentric(mat4 MVP, Image image, Sphere sphere, std::vector<std::vector<float>> depthBuffer) {
    for (int i = 0; i < gNumTriangles; ++i) {
        int k0 = gIndexBuffer[3 * i + 0];
        int k1 = gIndexBuffer[3 * i + 1];
        int k2 = gIndexBuffer[3 * i + 2];

        vec4 v0 = MVP * vec4(gVertexBuffer[k0], 1.0f);
        vec4 v1 = MVP * vec4(gVertexBuffer[k1], 1.0f);
        vec4 v2 = MVP * vec4(gVertexBuffer[k2], 1.0f);

        v0 /= v0.w;
        v1 /= v1.w;
        v2 /= v2.w;

        vec2 p0 = vec2(v0.x, v0.y);
        vec2 p1 = vec2(v1.x, v1.y);
        vec2 p2 = vec2(v2.x, v2.y);

        float z0 = v0.z;
        float z1 = v1.z;
        float z2 = v2.z;

        int minX = (int)floor(std::min({ p0.x, p1.x, p2.x }));
        int maxX = (int)ceil(std::max({ p0.x, p1.x, p2.x }));
        int minY = (int)floor(std::min({ p0.y, p1.y, p2.y }));
        int maxY = (int)ceil(std::max({ p0.y, p1.y, p2.y }));

        float beta_n = ((p0.y - p2.y) * p1.x + (p2.x - p0.x) * p1.y + p0.x * p2.y - p2.x * p0.y);
        float gamma_n = ((p0.y - p1.y) * p2.x + (p1.x - p0.x) * p2.y + p0.x * p1.y - p1.x * p0.y);
        float beta_x = (p0.y - p2.y) / beta_n;
        float beta_y = (p2.x - p0.x) / beta_n;
        float gamma_x = (p0.y - p1.y) / gamma_n;
        float gamma_y = (p1.x - p0.x) / gamma_n;
        float beta = ((p0.y - p2.y) * minX + (p2.x - p0.x) * minY + p0.x * p2.y - p2.x * p0.y) / beta_n;
        float gamma = ((p0.y - p1.y) * minX + (p1.x - p0.x) * minY + p0.x * p1.y - p1.x * p0.y) / gamma_n;
        int k = (maxX - minX) + 1;
        for (int y = minY; y <= maxY; y++) {
            for (int x = minX; x <= maxX; x++) {
                if (beta > 0.0f && gamma > 0.0f && (beta + gamma) < 1.0f) {
                    float alpha = 1.0f - beta - gamma;
                    float z = alpha * z0 + beta * z1 + gamma * z2;

                    if (z < depthBuffer[y][x]) {
                        depthBuffer[y][x] = z;
                        image.set(x, y, sphere.color);
                    }
                }
                beta += beta_x;
                gamma += gamma_x;
            }
            beta += beta_y - k * beta_x;
            gamma += gamma_y - k * gamma_x;
        }
    }
}

void render() {
    OutputImage.resize(Width * Height * 3, 0.0f);
    std::vector<std::vector<float>> depthBuffer(Height, std::vector<float>(Width, std::numeric_limits<float>::infinity()));

    Sphere sphere(vec3(0.0f, 0.0f, -7.0f), 2.0f, vec3(1.0f)); // sphere define
    Camera camera(vec3(0.0f, 0.0f, 0.0f), vec3(0.0f, 0.0f, -1.0f));//camera define
    std::vector<Surface*> surfaces = { &sphere };
    Scene scene(surfaces, camera);
    Image image;
    create_scene();  // 버퍼는 전역 변수로 생성됨

    //1 cal Mm
    mat4 Mm = cal_modelTransform(sphere);

    //2,3,4 calculate MVP
    mat4 MVP = cal_mvp(camera, Mm);

    //Baycentric algoritm
    cal_baycentric(MVP, image, sphere, depthBuffer);

    //delete buffer
    delete[] gVertexBuffer;
    delete[] gIndexBuffer;
}


void resize_callback(GLFWwindow*, int nw, int nh)
{
    //This is called in response to the window resizing.
    //The new width and height are passed in so we make 
    //any necessary changes:
    Width = nw;
    Height = nh;
    //Tell the viewport to use all of our screen estate
    glViewport(0, 0, nw, nh);

    //This is not necessary, we're just working in 2d so
    //why not let our spaces reflect it?
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(0.0, static_cast<double>(Width)
        , 0.0, static_cast<double>(Height)
        , 1.0, -1.0);

    //Reserve memory for our render so that we don't do 
    //excessive allocations and render the image
    OutputImage.reserve(Width * Height * 3);
    render();
}

int main(int argc, char* argv[])
{
    // -------------------------------------------------
    // Initialize Window
    // -------------------------------------------------

    GLFWwindow* window;

    /* Initialize the library */
    if (!glfwInit())
        return -1;

    /* Create a windowed mode window and its OpenGL context */
    window = glfwCreateWindow(Width, Height, "OpenGL Viewer", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    //We have an opengl context now. Everything from here on out 
    //is just managing our window or opengl directly.

    //Tell the opengl state machine we don't want it to make 
    //any assumptions about how pixels are aligned in memory 
    //during transfers between host and device (like glDrawPixels(...) )
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);

    //We call our resize function once to set everything up initially
    //after registering it as a callback with glfw
    glfwSetFramebufferSizeCallback(window, resize_callback);
    resize_callback(NULL, Width, Height);

    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window))
    {
        //Clear the screen
        glClear(GL_COLOR_BUFFER_BIT);

        // -------------------------------------------------------------
        //Rendering begins!
        glDrawPixels(Width, Height, GL_RGB, GL_FLOAT, &OutputImage[0]);
        //and ends.
        // -------------------------------------------------------------

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();

        //Close when the user hits 'q' or escape
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS
            || glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
        {
            glfwSetWindowShouldClose(window, GL_TRUE);
        }
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

