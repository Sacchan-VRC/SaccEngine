//
#define _USE_MATH_DEFINES
#include <math.h>
#define GLAD_GL_IMPLEMENTATION
#include <glad/glad.c>
#include <glad/glad.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <unistd.h>
#include <ctime>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <vector>
#include <learnopengl/shader.h>
#include <iostream>
#include <fstream>
#include <linmath.h>
#include <random>
#include <glm/gtx/euler_angles.hpp>
#include <glm/gtx/quaternion.hpp>
#include "SaccEngine.h"

#define CGLTF_IMPLEMENTATION
#include <cgltf.h>

bool fullScreen = false;
float FoV = 90;
double FPS_Limit = 200;
double FPS_Limit_FrameLen = (double)1 / FPS_Limit;
double deltaTime;

glm::vec3 gravity = glm::vec3(0, -9.80665, 0);

unsigned int GroundObject;
unsigned int colliderShip;

std::vector<myObject> *allObjects_Ptr;
std::vector<unsigned int> physicsObjs;
std::vector<unsigned int> staticColliderObjs;

glm::vec3 *camPos_Ptr;
glm::vec3 *camRot_Ptr;

typedef struct loadedTex
{
    unsigned int texID;
    std::string filename;
} loadedTex;

std::vector<loadedTex> allLoadedTextures;

typedef struct collisionTriangle
{
    glm::vec3 pos[3];
} collisionTriangle;

typedef struct Vertex
{
    glm::vec3 pos;
    glm::vec4 col;
    glm::vec3 norm;
    glm::vec2 uv;
} Vertex;

typedef struct convexMesh
{
    glm::vec3 verts;
    std::vector<unsigned short> indices;
} convexMesh;

typedef struct myMesh
{
    std::vector<Vertex> verts;
    std::vector<unsigned short> indices;
    GLuint vertex_array;
    glm::vec3 bounds;
    glm::vec3 boundsCenter;
    unsigned int submeshIndex;
    unsigned int numVertices;
    unsigned int numIndices;
    bool isConvex;
    convexMesh *convexVersion;
    myMesh() : indices(0), numVertices(0), numIndices(0) {}
} myMesh;

typedef struct DrawCall
{
    myMesh *mesh;
    Shader shader;
    loadedTex textures[16];
    unsigned int numTextures = 0;
    DrawCall() : shader("", ""), numTextures{0} {}
} DrawCall;

typedef struct myFullMesh
{
    std::vector<myMesh> subMeshes;
    std::string filename;
    unsigned int meshID;
} myFullMesh;

typedef struct myObject
{
    myFullMesh *mesh;
    physicsProperties *physics = nullptr;
    glm::vec3 position;
    glm::quat rotation;
    glm::vec3 scale;
    std::vector<DrawCall> dcs;
    unsigned int ID;
    myObject() : mesh(nullptr), physics(nullptr), position{0.0f, 0.0f, 0.0f}, rotation{0.0f, 0.0f, 0.0f, 0.0f}, scale{1.0f, 1.0f, 1.0f} {}
} myObject;

struct inertiaTensor
{
    glm::mat3 values;
    inertiaTensor(
        double xx, double xy, double xz,
        double yx, double yy, double yz,
        double zx, double zy, double zz)
        : values{{xx, xy, xz}, {yx, yy, yz}, {zx, zy, zz}} {}
};

typedef struct physicsProperties
{
    unsigned int myObj;
    int Collidertype = -1; // negative means static object
    myMesh *colliderMesh;
    glm::vec3 colliderScale = glm::vec3(1, 1, 1);
    float mass = 1.0f;
    float friction = 0.05f;
    float restitution = 0.5f;
    float Drag = 0.1f;
    float angDrag = 0.1f;
    glm::vec3 physPosition;
    glm::quat physRotation;
    glm::vec3 velocity;
    glm::vec3 AngVel;
    inertiaTensor IT = inertiaTensor(0, 0, 0, 0, 0, 0, 0, 0, 0);
} physicsProperties;

loadedTex loadTexture(const char *filePath)
{
    for (size_t i = 0; i < allLoadedTextures.size(); i++)
    {
        if (filePath == allLoadedTextures[i].filename)
        {
            return allLoadedTextures[i];
        }
    }

    loadedTex result;
    result.filename = filePath;
    stbi_set_flip_vertically_on_load(true);
    unsigned int texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    // set the texture wrapping/filtering options (on the currently bound texture object)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    // load and generate the texture
    int width2, height2, nrChannels2;
    unsigned char *data = stbi_load(filePath, &width2, &height2, &nrChannels2, 0);
    if (data)
    {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width2, height2, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);
    }
    else
    {
        std::cout << "Failed to load texture" << std::endl;
    }
    stbi_image_free(data);
    result.texID = texture;
    allLoadedTextures.push_back(result);
    return result;
}

DrawCall createDrawCall(
    const char *vertLoc, const char *fragLoc,
    myMesh *mesh,
    std::vector<std::string> texLoc)
{
    DrawCall dc;
    Shader ourShader(vertLoc, fragLoc);
    dc.shader = ourShader;
    dc.mesh = mesh;

    for (size_t i = 0; i < texLoc.size(); i++)
    {
        dc.textures[i] = loadTexture(texLoc[i].c_str());
        dc.numTextures++; // should check load success first
    }

    return dc;
}
static void error_callback(int error, const char *description)
{
    fprintf(stderr, "Error: %s\n", description);
}
void reloadShaders() { std::cout << "Reloading Shaders" << std::endl; }

glm::quat getQuaternionFromEuler(const glm::vec3 &eulerAngles)
{
    glm::quat qPitch = glm::angleAxis(eulerAngles.x, glm::vec3(1.0f, 0.0f, 0.0f));
    glm::quat qYaw = glm::angleAxis(eulerAngles.y, glm::vec3(0.0f, 1.0f, 0.0f));
    glm::quat qRoll = glm::angleAxis(eulerAngles.z, glm::vec3(0.0f, 0.0f, 1.0f));
    return qPitch * qYaw * qRoll;
}
bool check = false;
static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);

    if (key == GLFW_KEY_X && action == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, GL_TRUE);
    }
    if (key == GLFW_KEY_N && action == GLFW_PRESS)
    {
        reloadShaders();
    }
    if (key == GLFW_KEY_F && action == GLFW_PRESS)
    {
        if (fullScreen)
        {
            glfwSetWindowMonitor(window, NULL, 500, 500, 640, 480, 0);
        }
        else
        {
            glfwSetWindowMonitor(window, glfwGetPrimaryMonitor(), 0, 0, 2560, 1440, -1);
        }
        fullScreen = !fullScreen;
    }

    if (key == GLFW_KEY_H && action == GLFW_PRESS)
    {
        unsigned int newBowl = createObj("mesh/testbowl.glb", allObjects_Ptr);
        (*allObjects_Ptr)[newBowl].physics = new physicsProperties();
        (*allObjects_Ptr)[newBowl].physics->Collidertype = -1;
        (*allObjects_Ptr)[newBowl].physics->physPosition = *camPos_Ptr + glm::vec3(0.0f, -1.0f, 0.0f);
        (*allObjects_Ptr)[newBowl].position = (*allObjects_Ptr)[newBowl].physics->physPosition;
        staticColliderObjs.push_back(newBowl); // this needs to be dynamic pointers or something
        std::cout << "created mesh submeshes = " << (*allObjects_Ptr)[newBowl].mesh->subMeshes.size() << std ::endl;
    }
    if (key == GLFW_KEY_C && action == GLFW_PRESS)
    {
        unsigned int newCube = createObj("mesh/testcube.glb", allObjects_Ptr);
        (*allObjects_Ptr)[newCube].physics = new physicsProperties();
        (*allObjects_Ptr)[newCube].physics->Collidertype = -1;
        (*allObjects_Ptr)[newCube].physics->physPosition = *camPos_Ptr + glm::vec3(0.0f, -1.0f, 0.0f);
        (*allObjects_Ptr)[newCube].position = (*allObjects_Ptr)[newCube].physics->physPosition;
        staticColliderObjs.push_back(newCube); // this needs to be dynamic pointers or something
    }
    if (key == GLFW_KEY_N && action == GLFW_PRESS)
    {
        unsigned int newPlane = createObj("mesh/testplane.glb", allObjects_Ptr);
        (*allObjects_Ptr)[newPlane].physics = new physicsProperties();
        (*allObjects_Ptr)[newPlane].physics->Collidertype = -1;
        staticColliderObjs.push_back(newPlane); // this needs to be dynamic pointers or something

        (*allObjects_Ptr)[newPlane].physics->physPosition = *camPos_Ptr + glm::vec3(0.0f, -1.0f, 0.0f);
        (*allObjects_Ptr)[newPlane].scale = (glm::vec3(24.0f, 24.0f, 24.0f));
        (*allObjects_Ptr)[newPlane].physics->colliderScale = (*allObjects_Ptr)[newPlane].scale;
        (*allObjects_Ptr)[newPlane].position = (*allObjects_Ptr)[newPlane].physics->physPosition;
        glm::quat camQuat = glm::quat(glm::vec3(0.0f, (*camRot_Ptr).y, 0.0f));
        camQuat = glm::quat(glm::vec3((*camRot_Ptr).x, 0.0f, 0.0f)) * camQuat;
        const glm::vec3 forward = glm::vec3(0.0f, 0.0f, -1.0f) * camQuat;
        (*allObjects_Ptr)[newPlane].physics->velocity = forward * 10.0f;
    }
    if (key == GLFW_KEY_M && action == GLFW_PRESS)
    {
        unsigned int newShip = createObj("mesh/shipmulti.glb", allObjects_Ptr);
        (*allObjects_Ptr)[newShip].physics = new physicsProperties();
        (*allObjects_Ptr)[newShip].physics->Collidertype = -1;
        staticColliderObjs.push_back(newShip); // this needs to be dynamic pointers or something

        (*allObjects_Ptr)[newShip].physics->physPosition = *camPos_Ptr + glm::vec3(0.0f, -1.0f, 0.0f);
        (*allObjects_Ptr)[newShip].position = (*allObjects_Ptr)[newShip].physics->physPosition;
        glm::quat camQuat = getQuaternionFromEuler(*camRot_Ptr);
        const glm::vec3 forward = glm::vec3(0.0f, 0.0f, -1.0f) * camQuat;
        (*allObjects_Ptr)[newShip].scale = glm::vec3(10.0f, 10.0f, 10.0f);
        (*allObjects_Ptr)[newShip].physics->colliderScale = (*allObjects_Ptr)[newShip].scale;
        (*allObjects_Ptr)[newShip].physics->physRotation = glm::inverse(camQuat);
        (*allObjects_Ptr)[newShip].rotation = glm::inverse(camQuat);
        (*allObjects_Ptr)[newShip].physics->velocity = forward * 10.0f;
    }
    if (key == GLFW_KEY_R && action == GLFW_PRESS)
    {
        check = !check;
    }
    if (key == GLFW_KEY_T && action == GLFW_PRESS)
    {
        // debug messeges for testing vector push back causing corruption
        if (physicsObjs.size() > 0)
            std::cout << "checking mesh submeshes0 = " << (*allObjects_Ptr)[physicsObjs[0]].mesh->subMeshes.size() << std ::endl;
        if (physicsObjs.size() > 1)
            std::cout << "checking mesh submeshes1 = " << (*allObjects_Ptr)[physicsObjs[1]].mesh->subMeshes.size() << std ::endl;
        if (physicsObjs.size() > 2)
            std::cout << "checking mesh submeshes2 = " << (*allObjects_Ptr)[physicsObjs[2]].mesh->subMeshes.size() << std ::endl;
        if (physicsObjs.size() > 3)
            std::cout << "checking mesh submeshes3 = " << (*allObjects_Ptr)[physicsObjs[3]].mesh->subMeshes.size() << std ::endl;
    }
}

void mouseRelative(GLFWwindow *window, double *xpos, double *ypos)
{
    double x, y;
    glfwGetCursorPos(window, &x, &y);
    int windowWidth, windowHeight;
    glfwGetWindowSize(window, &windowWidth, &windowHeight);
    *xpos = (double)x / windowWidth;
    *ypos = (double)y / windowHeight;
}

char *read_file(const std::string &filename)
{
    std::ifstream theFile(filename.c_str(), std::ios::binary);

    if (!theFile.is_open())
    {
        std::cerr << filename + " : Error opening the file." << std::endl;
        return nullptr;
    }

    theFile.seekg(0, std::ios::end);
    size_t fileSize = theFile.tellg();
    theFile.seekg(0, std::ios::beg);

    char *content = new char[fileSize + 1];
    theFile.read(content, fileSize);
    content[fileSize] = '\0';

    char *result = new char[fileSize + 1];
    strcpy(result, content);

    delete[] content;

    return result;
}

glm::vec3 projectToPoint(glm::vec3 vec, glm::vec3 point)
{
    glm::vec3 inVec = glm::normalize(vec);
    float distToPoint = glm::dot(inVec, point);
    return inVec * distToPoint;
}

glm::vec3 getPointVelocity(const myObject *myObj, glm::vec3 pos)
{
    if (glm::length2(myObj->physics->AngVel) == 0)
        return myObj->physics->velocity;
    pos -= myObj->physics->physPosition;
    float distFromCOM = glm::length(pos);

    float circum = distFromCOM * M_PI * 2;
    glm::vec3 rotAxis = glm::normalize(myObj->physics->AngVel);
    float RPS = glm::degrees(glm::dot(rotAxis, myObj->physics->AngVel)) / 360.0f;
    float surfaceSpeed = RPS * circum;
    // Calculate the velocity of a point on the object relative to its center of mass
    glm::vec3 pointVelocity = myObj->physics->velocity + (glm::cross(glm::normalize(pos), rotAxis) * surfaceSpeed);

    return pointVelocity;
}

bool Blast = false;
void addForceAtPosition(myObject *physObj, glm::vec3 forceVector, glm::vec3 position)
{
    glm::vec3 localPosition = position - physObj->physics->physPosition;
    glm::vec3 torque = glm::cross(forceVector, localPosition); // momentArm = localPosition
    glm::vec3 angularAcceleration = (glm::inverse(physObj->physics->IT.values) * (torque));
    physObj->physics->AngVel += angularAcceleration;
    physObj->physics->velocity += forceVector;
}

void Inputs_MoveObject(GLFWwindow *window, float dt, glm::vec3 *campos, glm::vec3 *camrot, myObject *physObj)
{
    float moveStrength = glm::length(gravity);

    int state = glfwGetKey(window, GLFW_KEY_LEFT_SHIFT);
    if (state == GLFW_PRESS)
    {
        moveStrength *= 2;
    }
    state = glfwGetKey(window, GLFW_KEY_I);
    if (state == GLFW_PRESS)
    {
        physObj->physics->velocity += glm::vec3(0.0f, 0.0f, -1.0f) * moveStrength * dt;
    }
    state = glfwGetKey(window, GLFW_KEY_K);
    if (state == GLFW_PRESS)
    {
        physObj->physics->velocity += glm::vec3(0.0f, 0.0f, 1.0f) * moveStrength * dt;
    }
    state = glfwGetKey(window, GLFW_KEY_J);
    if (state == GLFW_PRESS)
    {
        physObj->physics->velocity += glm::vec3(-1.0f, 0.0f, 0.0f) * moveStrength * dt;
    }
    state = glfwGetKey(window, GLFW_KEY_L);
    if (state == GLFW_PRESS)
    {
        physObj->physics->velocity += glm::vec3(1.0f, 0.0f, 0.0f) * moveStrength * dt;
    }
    state = glfwGetKey(window, GLFW_KEY_U);
    if (state == GLFW_PRESS)
    {
        physObj->physics->velocity += glm::vec3(0.0f, 1.0f, 0.0f) * moveStrength * dt;
    }
    state = glfwGetKey(window, GLFW_KEY_O);
    if (state == GLFW_PRESS)
    {
        physObj->physics->velocity += glm::vec3(0.0f, -1.0f, 0.0f) * moveStrength * dt;
    }
    state = glfwGetKey(window, GLFW_KEY_B);
    if (state == GLFW_PRESS)
    {
        if (!Blast)
            addForceAtPosition(physObj, glm::vec3(0.0f, 0.0f, 15.0f), physObj->physics->physPosition + glm::vec3(0, -1, 0) * physObj->physics->colliderScale.x * 0.5f);
        Blast = true;
    } //
    else
    {
        Blast = false;
    }
}

float lastG;
double mouseXabsLast, mouseYabsLast;
void Inputs(GLFWwindow *window, float dt)
{
    float mouseSens = 0.00001f;
    float moveSpdWalk = 5 * dt;
    float moveSpdRun = 15 * dt; //
    float moveSpdCur = moveSpdWalk;
    double mouseXabs, mouseYabs;
    double mouseX, mouseY;
    glfwGetCursorPos(window, &mouseXabs, &mouseYabs);
    mouseX = mouseXabs - mouseXabsLast;
    mouseY = mouseYabs - mouseYabsLast;
    mouseXabsLast = mouseXabs;
    mouseYabsLast = mouseYabs;
    camRot_Ptr->y += mouseX * mouseSens * FoV;
    camRot_Ptr->x += mouseY * mouseSens * FoV;
    // FPS cameras rotate around Y then X
    glm::quat camQuat = getQuaternionFromEuler(*camRot_Ptr);
    int state = glfwGetKey(window, GLFW_KEY_LEFT_SHIFT);
    if (state == GLFW_PRESS)
    {
        moveSpdCur = moveSpdRun;
    }
    state = glfwGetKey(window, GLFW_KEY_W);
    if (state == GLFW_PRESS)
    {
        const glm::vec3 forward = glm::vec3(0.0f, 0.0f, -1.0f) * camQuat;
        *camPos_Ptr += forward * moveSpdCur;
    }
    state = glfwGetKey(window, GLFW_KEY_S);
    if (state == GLFW_PRESS)
    {
        const glm::vec3 back = glm::vec3(0.0f, 0.0f, 1.0f) * camQuat;
        *camPos_Ptr += back * moveSpdCur;
    }
    state = glfwGetKey(window, GLFW_KEY_A);
    if (state == GLFW_PRESS)
    {
        const glm::vec3 left = glm::vec3(-1.0f, 0.0f, 0.0f) * camQuat;
        *camPos_Ptr += left * moveSpdCur;
    }
    state = glfwGetKey(window, GLFW_KEY_D);
    if (state == GLFW_PRESS)
    {
        const glm::vec3 right = glm::vec3(1.0f, 0.0f, 0.0f) * camQuat;
        *camPos_Ptr += right * moveSpdCur;
    }
    state = glfwGetKey(window, GLFW_KEY_E);
    if (state == GLFW_PRESS)
    {
        const glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);
        *camPos_Ptr += up * moveSpdCur;
    }
    state = glfwGetKey(window, GLFW_KEY_Q);
    if (state == GLFW_PRESS)
    {
        const glm::vec3 down = glm::vec3(0.0f, -1.0f, 0.0f);
        *camPos_Ptr += down * moveSpdCur;
    }
    state = glfwGetKey(window, GLFW_KEY_EQUAL);
    if (state == GLFW_PRESS)
    {
        FoV -= (3 * FoV) * dt;
    }
    state = glfwGetKey(window, GLFW_KEY_MINUS);
    if (state == GLFW_PRESS)
    {
        FoV += (3 * FoV) * dt;
    }

    if (glfwGetTime() - lastG > 0.1f)
    {
        lastG = glfwGetTime();
        state = glfwGetKey(window, GLFW_KEY_G);
        std::random_device rd;
        std::mt19937 gen(rd()); //
        std::uniform_real_distribution<double> dis(0.4, 1.0);
        if (state == GLFW_PRESS)
        {
            unsigned int newSphere = createObj("mesh/testsphere.glb", allObjects_Ptr);
            (*allObjects_Ptr)[newSphere].physics = new physicsProperties();
            (*allObjects_Ptr)[newSphere].physics->Collidertype = 0;
            physicsObjs.push_back(newSphere); // this needs to be dynamic pointers or something

            float genrand = dis(gen);
            (*allObjects_Ptr)[newSphere].physics->colliderScale = glm::vec3(genrand, genrand, genrand);
            (*allObjects_Ptr)[newSphere].physics->mass = (4 / 3) * 1.0f * glm::pow(genrand, 3);
            (*allObjects_Ptr)[newSphere].scale = glm::vec3((*allObjects_Ptr)[newSphere].physics->colliderScale.x, (*allObjects_Ptr)[newSphere].physics->colliderScale.x, (*allObjects_Ptr)[newSphere].physics->colliderScale.x);
            (*allObjects_Ptr)[newSphere].physics->physPosition = *camPos_Ptr + glm::vec3(0.0f, -1.0f, 0.0f);
            (*allObjects_Ptr)[newSphere].position = (*allObjects_Ptr)[newSphere].physics->physPosition;
            float r = (*allObjects_Ptr)[newSphere].scale.x;
            float m = (*allObjects_Ptr)[newSphere].physics->mass;
            (*allObjects_Ptr)[newSphere].physics->IT = inertiaTensor(
                m * r * r / 2, 0, 0,
                0, m * r * r / 2, 0,
                0, 0, m * r * r / 2);
            glm::quat camQuat = getQuaternionFromEuler(*camRot_Ptr);
            const glm::vec3 forward = glm::vec3(0, 0, -1) * camQuat;
            (*allObjects_Ptr)[newSphere].physics->velocity = forward * 10.0f;
        }
        state = glfwGetKey(window, GLFW_KEY_V);
        if (state == GLFW_PRESS)
        {
            unsigned int newCube = createObj("mesh/testcube.glb", allObjects_Ptr);
            (*allObjects_Ptr)[newCube].physics = new physicsProperties();
            (*allObjects_Ptr)[newCube].physics->Collidertype = 1;
            (*allObjects_Ptr)[newCube].physics->colliderMesh = &allLoadedMeshes[loadMesh("mesh/testcube.glb")].subMeshes[0];

            float genrand = dis(gen);
            (*allObjects_Ptr)[newCube].physics->colliderScale = glm::vec3(genrand, genrand, genrand);
            double volume = glm::pow((*allObjects_Ptr)[newCube].scale.x, 3);
            double mass = 5 * volume;
            (*allObjects_Ptr)[newCube].physics->mass = mass;
            (*allObjects_Ptr)[newCube].scale = glm::vec3((*allObjects_Ptr)[newCube].physics->colliderScale.x, (*allObjects_Ptr)[newCube].physics->colliderScale.x, (*allObjects_Ptr)[newCube].physics->colliderScale.x);
            (*allObjects_Ptr)[newCube].physics->physPosition = *camPos_Ptr + glm::vec3(0.0f, -1.0f, 0.0f);
            (*allObjects_Ptr)[newCube].position = (*allObjects_Ptr)[newCube].physics->physPosition;

            float a = (*allObjects_Ptr)[newCube].scale.x;
            float b = (*allObjects_Ptr)[newCube].scale.x;
            float c = (*allObjects_Ptr)[newCube].scale.x;
            float m = (*allObjects_Ptr)[newCube].physics->mass;
            (*allObjects_Ptr)[newCube].physics->IT = inertiaTensor(
                1 / 12.0f * m * (a * a + b * b), 0, 0,
                0, 1 / 12.0f * m * (a * a + c * c), 0,
                0, 0, 1 / 12.0f * m * (b * b + c * c));

            physicsObjs.push_back(newCube); // this needs to be dynamic pointers or something

            glm::quat camQuat = glm::quat(glm::vec3(0.0f, (*camRot_Ptr).y, 0.0f));
            camQuat = glm::quat(glm::vec3((*camRot_Ptr).x, 0.0f, 0.0f)) * camQuat;
            const glm::vec3 forward = glm::vec3(0.0f, 0.0f, -1.0f) * camQuat;
            (*allObjects_Ptr)[newCube].physics->velocity = forward * 10.0f;
        }
        state = glfwGetKey(window, GLFW_KEY_Z);
        if (state == GLFW_PRESS)
        {
            unsigned int newCube = createObj("mesh/testcube.glb", allObjects_Ptr);
            (*allObjects_Ptr)[newCube].physics = new physicsProperties();
            (*allObjects_Ptr)[newCube].physics->Collidertype = 1;
            (*allObjects_Ptr)[newCube].physics->colliderMesh = &allLoadedMeshes[loadMesh("mesh/testcube.glb")].subMeshes[0];
            physicsObjs.push_back(newCube); // this needs to be dynamic pointers or something

            float genrand = dis(gen);
            // (*allObjects_Ptr)[newCube].physics->colliderScale = glm::vec3(genrand, genrand, genrand);
            (*allObjects_Ptr)[newCube].physics->mass = (4 / 3) * 1.0f * glm::pow(genrand, 3);
            (*allObjects_Ptr)[newCube].scale = glm::vec3((*allObjects_Ptr)[newCube].physics->colliderScale.x, (*allObjects_Ptr)[newCube].physics->colliderScale.x, (*allObjects_Ptr)[newCube].physics->colliderScale.x);
            (*allObjects_Ptr)[newCube].physics->physPosition = *camPos_Ptr + glm::vec3(0.0f, -1.0f, 0.0f);
            (*allObjects_Ptr)[newCube].position = (*allObjects_Ptr)[newCube].physics->physPosition;
            glm::quat camQuat = glm::quat(glm::vec3(0.0f, (*camRot_Ptr).y, 0.0f));
            camQuat = glm::quat(glm::vec3((*camRot_Ptr).x, 0.0f, 0.0f)) * camQuat;
            const glm::vec3 forward = glm::vec3(0.0f, 0.0f, -1.0f) * camQuat;
            // (*allObjects_Ptr)[newCube].physics->velocity = forward * 10.0f;
        }
    }
}

myFullMesh loadGLTF(std::string file)
{
    myFullMesh result;
    result.filename = file;
    cgltf_options options = {cgltf_file_type_invalid};
    cgltf_data *data = NULL;
    cgltf_result loadResult = cgltf_parse_file(&options, file.c_str(), &data);
    cgltf_load_buffers(&options, data, file.c_str());
    if (loadResult == cgltf_result_success)
    {
        auto primCount = &data->meshes->primitives_count;
        for (cgltf_size i = 0; i < *primCount; i++)
        {
            myMesh thisMesh;
            bool hasError = false;
            cgltf_primitive *primitive = &data->meshes->primitives[i];
            // std::cout << "ATTCOUNT: " << primitive->attributes_count << std::endl;
            if (primitive->type == cgltf_primitive_type_triangles)
            {
                thisMesh.submeshIndex = i;
                for (int i = 0; i < primitive->attributes_count; i++)
                {
                    cgltf_attribute attr2 = primitive->attributes[i];
                    if (attr2.data->name != nullptr)
                        std::cout << "attr.data->name:" << attr2.data->name << std::endl;
                }
                for (int i = 0; i < primitive->attributes_count; i++)
                {
                    cgltf_attribute attr = primitive->attributes[i];
                    switch (attr.type)
                    {
                    case cgltf_attribute_type_position:
                    {
                        const cgltf_accessor *position_accessor = attr.data;
                        unsigned int posCount = attr.data->count;
                        float vertPos[posCount * 3];
                        // std::cout << "attIndex: " << data->meshes->primitives->attributes->index << std::endl;
                        // std::cout << "vert_count " << data->meshes->primitives->attributes->data[cgltf_attribute_type_position].count << std::endl;
                        // std::cout << "vertex1 " << newVerts << std::endl;

                        cgltf_accessor_unpack_floats(position_accessor, &vertPos[0], posCount * 3);

                        // for (size_t i = 0; i < posCount; i++)
                        // {
                        //     std::cout << "vert:" << vertPos[i] << std::endl;
                        // }
                        // std::cout << "posCount: " << posCount << std::endl;

                        float lowestX, lowestY, lowestZ, highestX, highestY, highestZ;
                        lowestX = lowestY = lowestZ = std::numeric_limits<float>::max();
                        highestX = highestY = highestZ = std::numeric_limits<float>::min();
                        thisMesh.numVertices = posCount;
                        // size_t mostDistVertIndex = 0;
                        // float mostDistVert = 0.0f;
                        for (size_t i = 0, o = 0; o < posCount; o++, i += 3)
                        {
                            Vertex newVert;
                            newVert.pos[0] = vertPos[i];
                            newVert.pos[1] = vertPos[i + 1];
                            newVert.pos[2] = vertPos[i + 2];
                            thisMesh.verts.push_back(newVert);

                            if (newVert.pos[0] < lowestX)
                                lowestX = newVert.pos[0];
                            if (newVert.pos[1] < lowestY)
                                lowestY = newVert.pos[1];
                            if (newVert.pos[2] < lowestZ)
                                lowestZ = newVert.pos[2];
                            if (newVert.pos[0] > highestX)
                                highestX = newVert.pos[0];
                            if (newVert.pos[1] > highestY)
                                highestY = newVert.pos[1];
                            if (newVert.pos[2] > highestZ)
                                highestZ = newVert.pos[2];
                            thisMesh.bounds.x = highestX - lowestX;
                            thisMesh.bounds.y = highestY - lowestY;
                            thisMesh.bounds.z = highestZ - lowestZ;
                            thisMesh.boundsCenter.x = (highestX + lowestX) / 2;
                            thisMesh.boundsCenter.y = (highestY + lowestY) / 2;
                            thisMesh.boundsCenter.z = (highestZ + lowestZ) / 2;

                            // float vertDist = newVert.pos.x + newVert.pos.y + newVert.pos.z;
                            // if (mostDistVert < vertDist)
                            // {
                            //     mostDistVertIndex = thisMesh.verts.size() - 1;
                            //     mostDistVert = vertDist;
                            // }
                            // std::cout << "vert:" << o << ": " << newVerts[o].pos[0] << " : " << newVerts[o].pos[1] << " : " << newVerts[o].pos[2] << std::endl;
                        }

                        // thisMesh.mostDistantVert = thisMesh.verts[mostDistVertIndex].pos;

                        // indices
                        auto *index_accessor = primitive->indices;
                        thisMesh.numIndices = index_accessor->count;
                        // std::cout << "index_accessor->count: " << thisMesh.numIndices << std::endl;
                        if (index_accessor->component_type == cgltf_component_type_r_16u)
                        {
                            for (int indice = 0; indice < index_accessor->count; indice++)
                            {
                                unsigned short t_indice = cgltf_accessor_read_index(index_accessor, indice);
                                // std::cout << "t_indice: " << t_indice << std::endl;
                                thisMesh.indices.push_back(t_indice);
                            }
                        }
                        break;
                        // case cgltf_attribute_type_texcoord:
                        //     auto *uv_accessor = attr.data;
                        //     break;
                    }
                    case cgltf_attribute_type_normal:
                    {
                        const cgltf_accessor *norm_accessor = attr.data;
                        unsigned int normCount = attr.data->count;
                        float meshNorm[normCount * 3];
                        // std::cout << "attIndex: " << data->meshes->primitives->attributes->index << std::endl;
                        // std::cout << "vert_count " << data->meshes->primitives->attributes->data[cgltf_attribute_type_position].count << std::endl;
                        // std::cout << "vertex1 " << newVerts << std::endl;

                        cgltf_accessor_unpack_floats(norm_accessor, &meshNorm[0], normCount * 3);
                        for (size_t i = 0, o = 0; o < normCount; o++, i = i + 3)
                        {
                            thisMesh.verts[o].norm[0] = meshNorm[i];
                            thisMesh.verts[o].norm[1] = meshNorm[i + 1];
                            thisMesh.verts[o].norm[2] = meshNorm[i + 2];

                            // std::cout << "norm:" << o << ": " << thisMesh.verts[o].norm[0] << " : " << thisMesh.verts[o].norm[1] << " : " << thisMesh.verts[o].norm[2] << std::endl;
                        }
                        break;
                    }
                    case cgltf_attribute_type_texcoord:
                    {
                        const cgltf_accessor *uv_accessor = attr.data;
                        unsigned int uvCount = attr.data->count;
                        float texCoord[uvCount * 2];
                        // std::cout << "attIndex: " << data->meshes->primitives->attributes->index << std::endl;
                        // std::cout << "vert_count " << data->meshes->primitives->attributes->data[cgltf_attribute_type_position].count << std::endl;
                        // std::cout << "vertex1 " << newVerts << std::endl;

                        cgltf_accessor_unpack_floats(uv_accessor, &texCoord[0], uvCount * 2);
                        for (size_t i = 0, o = 0; o < uvCount; o++, i = i + 2)
                        {
                            thisMesh.verts[o].uv[0] = texCoord[i];
                            thisMesh.verts[o].uv[1] = texCoord[i + 1];

                            // std::cout << "vert:" << o << ": " << thisMesh.verts[o].uv[0] << " : " << thisMesh.verts[o].uv[1] << std::endl;
                        }
                        break;
                    }
                    case cgltf_attribute_type_color:
                    {
                        const cgltf_accessor *color_accessor = attr.data;
                        unsigned int colCount = attr.data->count;
                        float color[colCount * 4];

                        cgltf_accessor_unpack_floats(color_accessor, &color[0], colCount * 4);
                        for (size_t i = 0, o = 0; o < colCount; o++, i = i + 4)
                        {
                            thisMesh.verts[o].col[0] = color[i];
                            thisMesh.verts[o].col[1] = color[i + 1];
                            thisMesh.verts[o].col[2] = color[i + 2];
                            thisMesh.verts[o].col[3] = color[i + 3];
                            // std::cout << "vert:" << o << ": " << thisMesh.verts[o].col[0] << " : " << thisMesh.verts[o].col[1] << " : " << thisMesh.verts[o].col[2] << " : " << thisMesh.verts[o].col[3] << std::endl;
                        }
                        break;
                    }
                    }
                }
                bool convexCheck = true;
                for (size_t i = 0; i < thisMesh.verts.size(); i++)
                {
                    // should subtract CoM when its implemented
                    glm::vec3 thisVertDir = glm::normalize(thisMesh.verts[i].pos);
                    float thisVertMag = glm::length(thisMesh.verts[i].pos);
                    for (size_t u = i + 1; u < thisMesh.verts.size(); u++)
                    {
                        // should subtract CoM from thisMesh.verts[u].pos too
                        glm::vec3 thatVertDir = thisMesh.verts[i].pos;
                        if (glm::dot(thisVertDir, thatVertDir) > thisVertMag + 0.0000001f) // xD
                        {
                            // std::cout << "convexCheckDotResult =" << glm::dot(thisVertDir, thatVertDir) << std::endl;
                            convexCheck = false;
                            break;
                        }
                    }
                    if (!convexCheck)
                        break;
                }
                thisMesh.isConvex = convexCheck;
                // std::cout << (convexCheck ? "SubMesh is Convex" : "SubMesh isnt Convex") << std::endl;
            }
            if (!hasError)
            {
                // for (size_t i = 0; i < thisMesh.numVertices; i++)
                // {
                //     std::cout << "&allLoadedMeshes[0].verts: " << *(&thisMesh.verts[5].pos[0] + i) << std::endl;
                // }
                result.subMeshes.push_back(thisMesh);
            }
            else
            {
                std::cout << "subMeshes error" << std::endl;
            }
        }

        std::cout << "loaded model" << std::endl;
        /* TODO make awesome stuff */
        // cgltf_free(data);
    }
    return result;
}

unsigned int loadMesh(std::string file)
{
    for (size_t i = 0; i < allLoadedMeshes.size(); i++)
    {
        if ((allLoadedMeshes)[i].filename == file)
        {
            return i;
        }
    }

    // debug messeges for testing vector push back causing corruption
    // if (physicsObjs.size() > 0)
    //     std::cout << "checking mesh submeshes0 = " << (*allObjects_Ptr)[physicsObjs[0]].mesh->subMeshes.size() << std ::endl;
    // if (physicsObjs.size() > 1)
    //     std::cout << "checking mesh submeshes1 = " << (*allObjects_Ptr)[physicsObjs[1]].mesh->subMeshes.size() << std ::endl;
    // if (physicsObjs.size() > 2)
    //     std::cout << "checking mesh submeshes2 = " << (*allObjects_Ptr)[physicsObjs[2]].mesh->subMeshes.size() << std ::endl;
    // if (physicsObjs.size() > 3)
    //     std::cout << "checking mesh submeshes3 = " << (*allObjects_Ptr)[physicsObjs[3]].mesh->subMeshes.size() << std ::endl;
    allLoadedMeshes.push_back(loadGLTF(file));
    // if (physicsObjs.size() > 0)
    //     std::cout << "checking mesh submeshes0 = " << (*allObjects_Ptr)[physicsObjs[0]].mesh->subMeshes.size() << std ::endl;
    // if (physicsObjs.size() > 1)
    //     std::cout << "checking mesh submeshes1 = " << (*allObjects_Ptr)[physicsObjs[1]].mesh->subMeshes.size() << std ::endl;
    // if (physicsObjs.size() > 2)
    //     std::cout << "checking mesh submeshes2 = " << (*allObjects_Ptr)[physicsObjs[2]].mesh->subMeshes.size() << std ::endl;
    // if (physicsObjs.size() > 3)
    //     std::cout << "checking mesh submeshes3 = " << (*allObjects_Ptr)[physicsObjs[3]].mesh->subMeshes.size() << std ::endl;
    return allLoadedMeshes.size() - 1;
}

unsigned int createObj(std::string filename, std::vector<myObject> *allObjects)
{
    myObject newObj;

    // somehow this line overwrites submeshes of  (*allObjects_Ptr)[physicsObjs[0]].mesh->subMeshes.size()
    auto ass = loadMesh(filename);
    newObj.mesh = &(allLoadedMeshes)[ass];

    newObj.position = glm::vec3(0.0f, 0.0f, 0.0f);
    newObj.rotation = glm::vec3(0.0f, 0.0f, 0.0f);
    newObj.scale = glm::vec3(1.0f, 1.0f, 1.0f);
    newObj.ID = allObjects->size();
    for (size_t i = 0; i < newObj.mesh->subMeshes.size(); i++)
    {
        std::vector<std::string> tex1;
        tex1.push_back("tex/shit.png");
        tex1.push_back("tex/shit2.png");
        DrawCall newDC = createDrawCall("shaders/lighting.vs", "shaders/lighting.fs", &(newObj.mesh->subMeshes[i]), tex1);

        const GLint vpos_location = glGetAttribLocation(newDC.shader.ID, "vPos");
        const GLint vcol_location = glGetAttribLocation(newDC.shader.ID, "vCol");
        const GLint vnorm_location = glGetAttribLocation(newDC.shader.ID, "vNormal");

        GLuint vertex_buffer;
        glGenBuffers(1, &vertex_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
        glBufferData(GL_ARRAY_BUFFER, newDC.mesh->numVertices * sizeof(Vertex), &newDC.mesh->verts[0], GL_STATIC_DRAW);

        //
        glGenVertexArrays(1, &newDC.mesh->vertex_array);
        glBindVertexArray(newDC.mesh->vertex_array);
        glEnableVertexAttribArray(vpos_location);
        glVertexAttribPointer(vpos_location, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void *)(offsetof(Vertex, pos)));
        glEnableVertexAttribArray(vcol_location);
        glVertexAttribPointer(vcol_location, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void *)(offsetof(Vertex, col)));
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void *)(offsetof(Vertex, uv)));
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(vnorm_location, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void *)(offsetof(Vertex, norm)));
        glEnableVertexAttribArray(3);

        newObj.dcs.push_back(newDC);
    }

    allObjects->push_back(newObj);
    std::cout << "Created Obj# " << allObjects->size() - 1 << std::endl;
    return allObjects->size() - 1;
}
glm::quat worldSpaceRotate(glm::quat inQuat, glm::vec3 axis, float degrees)
{
    axis = glm::inverse(inQuat) * axis;
    axis *= degrees;
    inQuat = inQuat * glm::quat(glm::radians(axis));
    return inQuat;
}

float physDeltaTime = 0.02;
double physicsTime;

void convexmeshPlane(myObject *myObj, glm::vec3 norm, glm::vec3 position)
{
    // rotate plane to object space
    position -= myObj->physics->physPosition;
    glm::vec3 localNorm = glm::inverse(myObj->physics->physRotation) * norm;
    position = glm::inverse(myObj->physics->physRotation) * position;

    float lowestDot = FLT_MAX;
    size_t deepestVert = -1;
    std::vector<size_t> verts;
    float planeDist = glm::dot(localNorm, -position);
    bool collided = false;
    for (size_t i = 0; i < myObj->physics->colliderMesh->verts.size(); i++)
    {
        float thisDot = glm::dot(myObj->physics->colliderMesh->verts[i].pos, localNorm);
        if (-thisDot > planeDist)
        {
            verts.push_back(i);
            collided = true;
            if (thisDot < lowestDot)
            {
                lowestDot = thisDot;
            }
        }
    }
    // find how deep into plane it is, teleport out by that much
    if (collided)
    {
        myObj->physics->physPosition += norm * (-lowestDot - planeDist);
        glm::vec3 worldSpaceVertPos;
        for (size_t i = 0; i < verts.size(); i++)
        {
            worldSpaceVertPos += myObj->physics->physPosition + myObj->physics->physRotation * myObj->physics->colliderMesh->verts[verts[i]].pos;
        }
        worldSpaceVertPos /= verts.size();
        glm::vec3 pV = getPointVelocity(myObj, worldSpaceVertPos);
        float normDot = glm::dot(norm, -pV);
        glm::vec3 normVel = norm * normDot;
        glm::vec3 pVFlat = pV - normVel;
        glm::vec3 forceToAdd = (normVel + (normVel * myObj->physics->restitution) + -pVFlat * myObj->physics->friction);

        forceToAdd = glm::inverse(myObj->physics->physRotation) * forceToAdd;
        forceToAdd = myObj->physics->IT.values * forceToAdd;
        forceToAdd = myObj->physics->physRotation * forceToAdd;
        addForceAtPosition(myObj, forceToAdd, worldSpaceVertPos);

        // auto debugsphere = createObj("mesh/testsphere.glb", allObjects_Ptr);
        // (*allObjects_Ptr)[debugsphere].position = worldSpaceVertPos;
        // (*allObjects_Ptr)[debugsphere].scale = glm::vec3(0.3f, 0.3f, 0.3f);
    }
}

void spherePlane(myObject *myObj, glm::vec3 norm, glm::vec3 position)
{
    float radius = myObj->physics->colliderScale.x / 2;

    float posDot = glm::dot(norm, (myObj->physics->physPosition - norm * radius) - position); // check object is inside plane
    float normDot = glm::dot(-norm, myObj->physics->velocity);                                // don't collide if moving away from plane (prevents bounce in wrong direction, but can clip through)
    if (posDot < 0 && normDot > 0)
    {
        // std::cout << "spherePlane: Bounced" << std::endl;
        myObj->physics->physPosition -= norm * posDot;
        glm::vec3 point = myObj->physics->physPosition - norm * radius;
        glm::vec3 pV = getPointVelocity(myObj, point);
        glm::vec3 normVel = norm * normDot;
        glm::vec3 pVFlat = pV + normVel;
        glm::vec3 forceToAdd = (normVel + (normVel * myObj->physics->restitution) + -pVFlat * myObj->physics->friction);
        if (glm::length2(forceToAdd) > 0.00001f)
        {
            addForceAtPosition(myObj, forceToAdd, point);
        }
        // auto debugsphere = createObj("mesh/testsphere.glb", allObjects_Ptr);
        // (*allObjects_Ptr)[debugsphere].scale = 0.1f;
        // (*allObjects_Ptr)[debugsphere].position = myObj->physics->physPosition - norm * radius;
    }
}

void spherePoint(myObject *myObj, glm::vec3 point)
{
    spherePlane(myObj, glm::normalize(myObj->physics->physPosition - point), point);
}

// unused, broken
bool intersectsBounds(myObject *obj1, unsigned int meshIndex1, myObject *obj2, unsigned int meshIndex2)
{
    myMesh &mesh1 = obj1->mesh->subMeshes[meshIndex1];
    myMesh &mesh2 = obj2->mesh->subMeshes[meshIndex2];

    glm::vec3 boundsMin1 = obj1->physics->physPosition + mesh1.boundsCenter - mesh1.bounds / 2.0f; // min bounds of obj1
    glm::vec3 boundsMax1 = obj1->physics->physPosition + mesh1.boundsCenter + mesh1.bounds / 2.0f; // max bounds of obj1

    glm::vec3 boundsMin2 = obj2->physics->physPosition + mesh2.boundsCenter - mesh2.bounds / 2.0f; // min bounds of obj2
    glm::vec3 boundsMax2 = obj2->physics->physPosition + mesh2.boundsCenter + mesh2.bounds / 2.0f; // max bounds of obj2

    // auto debugsphere = createObj("mesh/cube.glb", allObjects_Ptr);
    // (*allObjects_Ptr)[debugsphere].position = obj1->physics->physPosition + mesh1.boundsCenter;
    // (*allObjects_Ptr)[debugsphere].scale = mesh1.bounds / 2.0f;
    // debugsphere = createObj("mesh/cube.glb", allObjects_Ptr);
    // (*allObjects_Ptr)[debugsphere].position = obj2->physics->physPosition + mesh2.boundsCenter;
    // (*allObjects_Ptr)[debugsphere].scale = mesh2.bounds / 2.0f;

    if (check)
    {
        auto debugsphere = createObj("mesh/testsphere.glb", allObjects_Ptr);
        (*allObjects_Ptr)[debugsphere].position = glm::vec3(boundsMin1.x, 0.0f, 0.0f);
        (*allObjects_Ptr)[debugsphere].scale = glm::vec3(0.3f, 0.3f, 0.3f);
        debugsphere = createObj("mesh/testsphere.glb", allObjects_Ptr);
        (*allObjects_Ptr)[debugsphere].position = glm::vec3(boundsMin2.x, 0.0f, 0.0f);
        (*allObjects_Ptr)[debugsphere].scale = glm::vec3(0.3f, 0.3f, 0.3f);
        debugsphere = createObj("mesh/testsphere.glb", allObjects_Ptr);
        (*allObjects_Ptr)[debugsphere].position = glm::vec3(boundsMax2.x, 0.0f, 0.0f);
        (*allObjects_Ptr)[debugsphere].scale = glm::vec3(0.3f, 0.3f, 0.3f);
    }
    // check if any of the axes intersect
    bool xIntersects = (boundsMin1.x >= boundsMin2.x && boundsMin1.x <= boundsMax2.x) || (boundsMax1.x >= boundsMin2.x && boundsMax1.x <= boundsMax2.x);
    bool yIntersects = (boundsMin1.y >= boundsMin2.y && boundsMin1.y <= boundsMax2.y) || (boundsMax1.y >= boundsMin2.y && boundsMax1.y <= boundsMax2.y);
    bool zIntersects = (boundsMin1.z >= boundsMin2.z && boundsMin1.z <= boundsMax2.z) || (boundsMax1.z >= boundsMin2.z && boundsMax1.z <= boundsMax2.z);

    if (check)
    {
        std::cout << xIntersects << std::endl;
    }

    return (xIntersects && yIntersects && zIntersects);
}

// unused, broken, might be faster than sphereInsideBounds if fixed
bool insideBounds2(myObject *myObj, unsigned int meshIndex, glm::vec3 point, float sphereRadius)
{
    myMesh &mesh = myObj->mesh->subMeshes[meshIndex];
    glm::vec3 localPoint = point - myObj->physics->physPosition;

    glm::vec3 boundsMin = (mesh.boundsCenter - mesh.bounds / 2.0f) - sphereRadius;
    glm::vec3 boundsMax = (mesh.boundsCenter + mesh.bounds / 2.0f) + sphereRadius;

    return (boundsMin.x <= localPoint.x && localPoint.x <= boundsMax.x) &&
           (boundsMin.y <= localPoint.y && localPoint.y <= boundsMax.y) &&
           (boundsMin.z <= localPoint.z && localPoint.z <= boundsMax.z);
}

bool sphereInsideBounds(myObject *myObj, unsigned int meshIndex, glm::vec3 point, float sphereRadius)
{
    // designed for use with objects in local space (that have been pre-rotated to local space of myObj)
    // uncomment physPosition and physRotation to use in world space
    // colliderScale is untested.
    myMesh &mesh = myObj->mesh->subMeshes[meshIndex];
    glm::vec3 right = /* myObj->physics->physRotation * */ glm::vec3(1, 0, 0);
    glm::vec3 up = /* myObj->physics->physRotation * */ glm::vec3(0, 1, 0);
    glm::vec3 forward = /* myObj->physics->physRotation * */ glm::vec3(0, 0, 1);

    glm::vec3 localboundsCenter = /* myObj->rotation * */ mesh.boundsCenter;
    glm::vec3 localPoint = point - (/* myObj->physics->physPosition + */ /* myObj->physics->colliderScale * */ localboundsCenter);
    glm::vec3 bounds = (mesh.bounds /* * myObj->physics->colliderScale */ * 0.5f) + sphereRadius;

    // debug bounding box
    // if (check)
    // {
    //     auto debugsphere = createObj("mesh/cube.glb", allObjects_Ptr);

    //     (*allObjects_Ptr)[debugsphere].scale = glm::vec3(bounds.x, bounds.y, bounds.z);
    //     (*allObjects_Ptr)[debugsphere].rotation = myObj->rotation;
    //     (*allObjects_Ptr)[debugsphere].position = myObj->position + localboundsCenter;
    // }

    float upDist = glm::dot(up, localPoint);
    bool insideY = glm::abs(upDist) < bounds.y;
    float rightDist = glm::dot(right, localPoint);
    bool insideX = glm::abs(rightDist) < bounds.x;
    float forwardDist = glm::dot(forward, localPoint);
    bool insideZ = glm::abs(forwardDist) < bounds.z;
    return insideY && insideX && insideZ;
}

// unused, dunno if working
bool cubeInside(myObject *myObj, glm::vec3 point)
{
    point -= myObj->physics->physPosition;
    glm::vec3 right = myObj->physics->physRotation * glm::vec3(1, 0, 0);
    glm::vec3 up = myObj->physics->physRotation * glm::vec3(0, 1, 0);
    glm::vec3 forward = myObj->physics->physRotation * glm::vec3(0, 0, 1);

    float upDist = glm::dot(up, point);
    bool insideX, insideY, insideZ;
    float rightDist = glm::dot(right, point);
    if (glm::abs(rightDist) < myObj->scale.x)
    {
        insideX = true;
    }
    if (glm::abs(upDist) < myObj->scale.y)
    {
        insideY = true;
    }
    float forwardDist = glm::dot(forward, point);
    if (glm::abs(forwardDist) < myObj->scale.z)
    {
        insideZ = true;
    }
    return (insideX + insideY + insideZ) == 3;
}

void cubeTri(myObject *myObj, collisionTriangle tri)
{
}

void sphereSphere(myObject *thisObj, myObject *thatObj)
{
    float dist = glm::distance(thisObj->physics->physPosition, thatObj->physics->physPosition);
    if (dist < (thisObj->physics->colliderScale.x + thatObj->physics->colliderScale.x) * 0.5f)
    {
        float collisionSpeed = glm::distance(thisObj->physics->velocity, thatObj->physics->velocity);
        glm::vec3 collisionNormal = glm::normalize(thisObj->physics->physPosition - thatObj->physics->physPosition);

        glm::vec3 res = ((thisObj->physics->colliderScale.x + thatObj->physics->colliderScale.x) * 0.5f - dist) * collisionNormal;

        float sizeRatio = thisObj->physics->colliderScale.x / (thisObj->physics->colliderScale.x + thatObj->physics->colliderScale.x);

        thisObj->physics->physPosition += res * (1 - sizeRatio);
        thatObj->physics->physPosition -= res * sizeRatio;

        glm::vec3 velDif = thisObj->physics->velocity - thatObj->physics->velocity;

        float velDot = glm::dot(velDif, collisionNormal);

        float massRatio = thisObj->physics->mass / (thisObj->physics->mass + thatObj->physics->mass);

        glm::vec3 reflection = collisionNormal * velDot;
        thisObj->physics->velocity -= reflection * (1 - massRatio);
        thatObj->physics->velocity += reflection * massRatio;

        // TODO: Friction
    }
}

void sphereTri(myObject *myObj, collisionTriangle tri)
{
    glm::vec3 nearestPoint = glm::vec3(0, 0, 0);
    float radius = myObj->physics->colliderScale.x / 2;
    glm::vec3 tripos0 = tri.pos[0];
    glm::vec3 tripos1 = tri.pos[1];
    glm::vec3 tripos2 = tri.pos[2];
    glm::vec3 edge0 = tripos0 - tripos2;
    glm::vec3 edge1 = tripos2 - tripos1;
    glm::vec3 edge2 = tripos1 - tripos0;
    glm::vec3 norm = glm::cross(tripos1 - tripos0, tripos2 - tripos1);
    glm::vec3 physPos = myObj->physics->physPosition;
    if (glm::dot(norm, physPos - tripos0) < 0)
    {
        return;
    }
    if (glm::dot(glm::cross(edge0, norm), physPos - tripos2) > 0)
    {
        float endOfEdge = glm::dot(edge0, edge0);
        float spherePos2 = glm::dot(edge0, physPos - tripos2);
        if (spherePos2 > endOfEdge)
        {
            nearestPoint = tripos0;
            // std::cout << "CollisionCheck: edge1 Point0" << std::endl;
        }
        else if (spherePos2 < 0)
        {
            nearestPoint = tripos2;
            // std::cout << "CollisionCheck: edge1 Point2" << std::endl;
        }
        else
        {
            glm::vec3 edge0Norm = glm::normalize(edge0);
            nearestPoint = tripos2 + edge0Norm * glm::dot(edge0Norm, physPos - tripos2);
            // std::cout << "CollisionCheck: edge1" << std::endl;
        }
        spherePoint(myObj, nearestPoint);
        return;
    }
    else if (glm::dot(glm::cross(edge1, norm), physPos - tripos1) > 0)
    {
        float endOfEdge = glm::dot(edge1, edge1);
        float spherePos = glm::dot(edge1, physPos - tripos1);
        if (spherePos > endOfEdge)
        {
            nearestPoint = tripos2;
            // std::cout << "CollisionCheck: edge1 Point2" << std::endl;
        }
        else if (spherePos < 0)
        {
            nearestPoint = tripos1;
            // std::cout << "CollisionCheck: edge1 Point1" << std::endl;
        }
        else
        {
            glm::vec3 edge1Norm = glm::normalize(edge1);
            nearestPoint = tripos1 + edge1Norm * glm::dot(edge1Norm, physPos - tripos1);
            // std::cout << "CollisionCheck: edge1" << std::endl;
        }
        spherePoint(myObj, nearestPoint);
        return;
    }
    else if (glm::dot(glm::cross(edge2, norm), physPos - tripos0) > 0)
    {
        float endOfEdge = glm::dot(edge2, edge2);
        float spherePos = glm::dot(edge2, physPos - tripos0);
        if (spherePos > endOfEdge)
        {
            nearestPoint = tripos1;
            // std::cout << "CollisionCheck: edge2 Point1" << std::endl;
        }
        else if (spherePos < 0)
        {
            nearestPoint = tripos0;
            // std::cout << "CollisionCheck: edge2 Point0" << std::endl;
        }
        else
        {
            glm::vec3 edge2Norm = glm::normalize(edge2);
            nearestPoint = tripos0 + edge2Norm * glm::dot(edge2Norm, physPos - tripos0);
            // std::cout << "CollisionCheck: edge2" << std::endl;
        }
        spherePoint(myObj, nearestPoint);
        return;
    }

    if (glm::dot(norm, physPos - tripos0) > 0)
    {
        norm = glm::normalize(norm);
        // std::cout << "CollisionCheck: Plane: " << std::endl;
        spherePlane(myObj, norm, tripos0);
    }
}

void convexMeshConvexMesh(myObject *thisMesh, myObject *thatMesh)
{
    // TODO:
}

double checkTimesAdded;
unsigned int numchecks;
void sphereMesh(myObject *thisSphere, myObject *thatMesh)
{
    glm::vec3 thatMeshPos = thatMesh->physics->physPosition;
    glm::quat thatMeshRot = thatMesh->physics->physRotation;
    physicsProperties &thisPhys = *thisSphere->physics;

    // rotate pos and vel to mesh local space so we don't have to rotate every tri
    thisPhys.velocity = glm::inverse(thatMeshRot) * thisPhys.velocity;
    thisPhys.AngVel = glm::inverse(thatMeshRot) * thisPhys.AngVel;
    thisPhys.physPosition = glm::inverse(thatMeshRot) * (thisPhys.physPosition - thatMeshPos);
    // I don't think theres a way to do this for scale, as a sphere is only defined by radius and we need to scale all axis
    // Could be done for spheres if all axis are equally scaled
    // wold require an if statement in the for loop, not sure if worth it

    for (size_t o = 0; o < thatMesh->mesh->subMeshes.size(); o++)
    {
        glm::vec3 &colScale = thatMesh->physics->colliderScale;
        if (!sphereInsideBounds(thatMesh, o, thisPhys.physPosition / colScale, thisPhys.colliderScale.x))
            continue;
        myMesh &theMesh = thatMesh->mesh->subMeshes[o];
        int numtris = theMesh.numIndices / 3;

        // double startTime = glfwGetTime();

        for (size_t i = 0; i < numtris; i++)
        {
            int index = i * 3;
            // commented stuff isn't needed because we rotate pos/vel of thisSphere instead
            // glm::vec3 thatMeshPos = thatMesh->physics->physPosition;
            collisionTriangle tri;
            tri.pos[0] = (/* thatMesh->rotation *  */ (theMesh.verts[theMesh.indices[index]].pos * colScale)) /* + thatMeshPos */;
            tri.pos[1] = (/* thatMesh->rotation *  */ (theMesh.verts[theMesh.indices[index + 1]].pos * colScale)) /* + thatMeshPos */;
            tri.pos[2] = (/* thatMesh->rotation *  */ (theMesh.verts[theMesh.indices[index + 2]].pos * colScale)) /* + thatMeshPos */;
            sphereTri(thisSphere, tri);
        }
        // double endTime = glfwGetTime();
        // if (check)
        // {
        //     checkTimesAdded += endTime - startTime;
        //     numchecks++;
        //     std::cout << checkTimesAdded / numchecks << std::endl;
        // }
    }
    // unrotate after
    thisPhys.velocity = thatMeshRot * thisPhys.velocity;
    thisPhys.AngVel = thatMeshRot * thisPhys.AngVel;
    thisPhys.physPosition = (thatMeshRot * thisPhys.physPosition) + thatMeshPos;
}

void updatePhysics(std::vector<myObject> *allObjects)
{
    physicsTime += physDeltaTime;
    for (size_t i = 0; i < physicsObjs.size(); i++)
    {
        myObject &thisObj = (*allObjects)[physicsObjs[i]];
        if (thisObj.physics->Collidertype < 0)
        {
            continue;
        }
        thisObj.physics->velocity += gravity * physDeltaTime;
        thisObj.physics->velocity = glm::lerp(thisObj.physics->velocity, glm::vec3(0, 0, 0), 1.0f - (float)glm::pow(0.5f, physDeltaTime * thisObj.physics->Drag));
        thisObj.physics->AngVel = glm::lerp(thisObj.physics->AngVel, glm::vec3(0, 0, 0), 1.0f - (float)glm::pow(0.5f, physDeltaTime * thisObj.physics->angDrag));
        thisObj.physics->physPosition += thisObj.physics->velocity * physDeltaTime;
        if (glm::length2(thisObj.physics->AngVel) > 0)
        {
            glm::vec3 rotation_vector = thisObj.physics->AngVel;
            glm::vec3 normalized_rotation_vector = glm::normalize(rotation_vector);
            float rotation_rate = glm::dot(normalized_rotation_vector, (rotation_vector));

            thisObj.physics->physRotation = glm::angleAxis(-rotation_rate * physDeltaTime, normalized_rotation_vector) * thisObj.physics->physRotation;
        }
        if (thisObj.physics->Collidertype == 1)
        {
            convexmeshPlane(&thisObj, glm::vec3(0, 1, 0), glm::vec3(0, 0, 0));
        }
        for (size_t o = 0; o < physicsObjs.size(); o++)
        {
            if (physicsObjs[i] == physicsObjs[o])
            {
                continue;
            }
            switch (thisObj.physics->Collidertype)
            {
            case 0: // sphere
            {
                switch ((*allObjects)[physicsObjs[o]].physics->Collidertype)
                {
                case 0: // sphere
                {
                    sphereSphere(&thisObj, &(*allObjects)[physicsObjs[o]]);
                    break;
                }
                case 1: // convexmesh
                {
                    // convexmeshConvexmesh(&thisObj, &(*allObjects)[physicsObjs[o]]);
                    break;
                }
                }
                break;
            }
            default:
                break;
            }
        }
        for (size_t o = 0; o < staticColliderObjs.size(); o++)
        {
            switch (thisObj.physics->Collidertype)
            {
            case 0: // static mesh
            {
                switch ((*allObjects)[staticColliderObjs[o]].physics->Collidertype)
                {
                case -1: // static mesh
                {
                    sphereMesh(&thisObj, &(*allObjects)[staticColliderObjs[o]]);
                    break;
                }
                }
            }
            default:
                break;
            }
        }
        thisObj.rotation = thisObj.physics->physRotation;
        thisObj.position = thisObj.physics->physPosition;
    }
}

void physicsExtrapolation(std::vector<myObject> *allObjects, float amount)
{
    for (size_t i = 0; i < physicsObjs.size(); i++)
    {
        myObject &thisObj = (*allObjects)[physicsObjs[i]];
        if (glm::length2(thisObj.physics->AngVel) > 0)
        {
            thisObj.position = thisObj.physics->physPosition + (thisObj.physics->velocity * amount);

            glm::vec3 rotation_vector = thisObj.physics->AngVel;
            glm::vec3 normalized_rotation_vector = glm::normalize(rotation_vector);
            float rotation_rate = glm::dot(normalized_rotation_vector, (rotation_vector));

            thisObj.rotation = glm::angleAxis((-rotation_rate * amount), normalized_rotation_vector) * thisObj.physics->physRotation;
        }
    }
}

int main(void)
{
    // TODO: fix values getting corrupted when vector length is increased
    allLoadedMeshes.reserve(1000);

    bool use_GLFW_DOUBLEBUFFER = true; // true needed for vsyncoff @ fullscreen

    glfwSetErrorCallback(error_callback);
    if (!glfwInit())
        exit(EXIT_FAILURE);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_DECORATED, true); // !borderless
    glfwWindowHint(GLFW_FLOATING, false); // alwaysontop
    // glfwWindowHint(GLFW_MAXIMIZED, GLFW_TRUE);
    glfwWindowHint(GLFW_TRANSPARENT_FRAMEBUFFER, false); // transparent
    glfwWindowHint(GLFW_MOUSE_PASSTHROUGH, false);
    glfwWindowHint(GLFW_POSITION_X, 500);
    glfwWindowHint(GLFW_POSITION_Y, 500);
    glfwWindowHint(GLFW_SAMPLES, 0);
    glfwWindowHint(GLFW_WIN32_KEYBOARD_MENU, true); // allow alt then spacebar to open menu
    glfwWindowHint(GLFW_RESIZABLE, true);
    if (!use_GLFW_DOUBLEBUFFER)
        glfwWindowHint(GLFW_DOUBLEBUFFER, GL_FALSE);

    GLFWwindow *window = glfwCreateWindow(640, 480, "OpenGL Triangle", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    if (glfwRawMouseMotionSupported())
        glfwSetInputMode(window, GLFW_RAW_MOUSE_MOTION, GLFW_TRUE);

    glfwSetKeyCallback(window, key_callback);

    glfwMakeContextCurrent(window);
    gladLoadGL();
    glfwSwapInterval(2);
    glEnable(GL_DEPTH_TEST);

    std::vector<myObject> allObjects;
    allObjects_Ptr = &allObjects;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(-100.0, 100.0);

    // 1000 balls
    // for (size_t i = 0; i < 1000; i++)
    // {
    //     auto physobj = createObj("mesh/testsphere.glb", &allObjects);
    //     allObjects[physobj].position = glm::vec3(0, 0, 5);
    //     allObjects[physobj].physics = new physicsProperties();
    //     allObjects[physobj].physics->physPosition = glm::vec3((int)(i / 100) * 2, 150 * (i / (float)100 / 10), (i * 2) % 200);
    //     // allObjects[physobj].physics->velocity = glm::vec3(0, 0, 1);
    //     physicsObjs.push_back(physobj); // this needs to be dynamic pointers or something
    // }

    // just one ball
    auto physobj1 = createObj("mesh/testsphere.glb", &allObjects);
    allObjects[physobj1].physics = new physicsProperties();
    allObjects[physobj1].physics->Collidertype = 0;
    allObjects[physobj1].physics->AngVel = glm::radians(glm::vec3(0, 0, 0));
    allObjects[physobj1].physics->physPosition = glm::vec3(0, 3, 0);
    physicsObjs.push_back(physobj1); // this needs to be dynamic pointers or something

    unsigned int numBoundTex = 0;

    glm::vec3 camPos = glm::vec3(0.0f, 0.0f, 3.0f), camRot, camScale;
    camPos_Ptr = &camPos;
    camRot_Ptr = &camRot;
    glm::mat4 view;
    glm::mat4 projection;

    double curTime = glfwGetTime();

    int numships = 0;

    glfwGetCursorPos(window, &mouseXabsLast, &mouseYabsLast);

    while (!glfwWindowShouldClose(window))
    {
        double now = glfwGetTime();
        deltaTime = now - curTime;
        if (deltaTime < FPS_Limit_FrameLen)
        {
            double expectedTime = curTime + FPS_Limit_FrameLen;
            now = expectedTime;
            deltaTime = now - curTime;
            float curTimeDif = now - expectedTime;
            // time is set to the time of the frame length set by the frame rate limit
            // we wait at the end of the frame so that we don't waste time in which we could be creating the frame
        }
        curTime = now;
        double curtimeminusphys = now - physDeltaTime;
        Inputs_MoveObject(window, deltaTime, &camPos, &camRot, &allObjects[physobj1]);
        while (physicsTime < curtimeminusphys)
        {
            updatePhysics(allObjects_Ptr);
        }
        physicsExtrapolation(&allObjects, curTime - physicsTime);

        std::stringstream ss;
        ss << "SaccEngine"
           << " [" << 1 / deltaTime << " FPS]";

        glfwSetWindowTitle(window, ss.str().c_str());

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        double xpos, ypos;
        mouseRelative(window, &xpos, &ypos);

        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        const float ratio = width / (float)height;

        glViewport(0, 0, width, height);
        glClear(GL_COLOR_BUFFER_BIT);

        mat4x4 m, p, mvp;

        glfwPollEvents();
        Inputs(window, deltaTime);
        if (numships < 100)
        {
            numships++;
            int newShip = createObj("mesh/shipmulti.glb", &allObjects);
            allObjects[newShip].position = glm::vec3(dis(gen), dis(gen), dis(gen));
            allObjects[newShip].physics = new physicsProperties();
            allObjects[newShip].physics->physPosition = allObjects[newShip].position;
            allObjects[newShip].physics->Collidertype = -1;
            staticColliderObjs.push_back(newShip);
        }

        for (size_t u = 0; u < allObjects.size(); u++)
        {
            for (size_t i = 0; i < allObjects[u].dcs.size(); i++)
            {
                allObjects[u].dcs[i].shader.use();
                float div = (1 - i * .15f);
                vec4 myCol = {1.0f * div, 1.0f * div, 1.0f * div, 1.0f * div};
                // myCol *= greenValue;
                allObjects[u].dcs[i].shader.setVec4("myColor", myCol);

                // mat4x4_identity(m);
                // mat4x4_rotate_Z(m, m, (float)glfwGetTime());
                // mat4x4_ortho(p, -ratio, ratio, -1.f, 1.f, 1.f, -1.f);
                // mat4x4_mul(mvp, p, m);

                allObjects[u].dcs[i].shader.setMat4x4("MVP", mvp);

                allObjects[u].dcs[i].shader.setInt("ourTexture", 0);
                allObjects[u].dcs[i].shader.setInt("ourTexture2", 1);
                // allObjects[u].dcs[i].shader.setInt("ourTexture3", 2);

                for (size_t o = 0; o < 16; o++)
                {
                    glActiveTexture(GL_TEXTURE0 + o);
                    if (o < allObjects[u].dcs[i].numTextures)
                    {
                        glBindTexture(GL_TEXTURE_2D, allObjects[u].dcs[i].textures[o].texID);
                    }
                    else if (o < numBoundTex)
                    {
                        glBindTexture(GL_TEXTURE_2D, 0);
                    }
                    else
                    {
                        numBoundTex = o;
                    }
                }

                glm::mat4 model = glm::mat4(1.0f);
                model = glm::translate(model, allObjects[u].position);
                glm::mat4 modelRot = glm::mat4(1.0f);
                glm::mat4 rotationMatrix = glm::toMat4(allObjects[u].rotation);

                modelRot = rotationMatrix;
                model *= modelRot;
                model = glm::scale(model, allObjects[u].scale);

                glm::mat4 view = glm::mat4(1.0f);
                projection = glm::perspective(glm::radians(FoV), ratio, 0.1f, 10000.0f);
                // view matrix is inverse camPos
                view = glm::translate(view, -camPos);
                // Rotation order is also inverse
                projection = glm::rotate(projection, camRot.x, glm::vec3(1.0, 0.0, 0.0));
                projection = glm::rotate(projection, camRot.y, glm::vec3(0.0, 1.0, 0.0));

                int modelLoc = glGetUniformLocation(allObjects[u].dcs[i].shader.ID, "model");
                glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));
                int viewLoc = glGetUniformLocation(allObjects[u].dcs[i].shader.ID, "view");
                glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
                int projectionLoc = glGetUniformLocation(allObjects[u].dcs[i].shader.ID, "projection");
                glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
                int rotationLoc = glGetUniformLocation(allObjects[u].dcs[i].shader.ID, "rotation");
                glUniformMatrix4fv(rotationLoc, 1, GL_FALSE, glm::value_ptr(modelRot));
                // int camPosLoc = glGetUniformLocation(allObjects[u].dcs[i].shader.ID, "cameraPos");
                // glUniform3fv(camPosLoc, 1, glm::value_ptr(camPos));

                glBindVertexArray(allObjects[u].dcs[i].mesh->vertex_array);

                glDrawElements(GL_TRIANGLES, allObjects[u].dcs[i].mesh->numIndices, GL_UNSIGNED_SHORT, &allObjects[u].dcs[i].mesh->indices[0]);
            }
        }
        if (use_GLFW_DOUBLEBUFFER)
            glfwSwapBuffers(window);
        else
            glFlush();
        // using Sleep for this does NOT work properly.
        while (glfwGetTime() < curTime)
        {
        }
    }

    glfwDestroyWindow(window);

    glfwTerminate();
    exit(EXIT_SUCCESS);
}