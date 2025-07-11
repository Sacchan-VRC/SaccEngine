#ifndef MAIN_H
#define MAIN_H

struct physicsProperties;
struct myObject;
unsigned int createObj(std::string filename, std::vector<myObject> *allObjects);
struct myFullMesh;
std::vector<myFullMesh> allLoadedMeshes;
unsigned int loadMesh(std::string file);
glm::vec3 getPointVelocity(const myObject *myObj, glm::vec3 pos);

#endif