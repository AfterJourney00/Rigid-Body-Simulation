#pragma once

#include <utility>
#include <vector>
#include <omp.h>
#include <glm/glm.hpp>

class Contact {
public:
    Contact () {};
    RigidBody *a;
    RigidBody *b;

    glm::vec3 particle_position;
    glm::vec3 face_normal;

    glm::vec3 edge1;
    glm::vec3 edge2;

    bool is_face_vertex;
    bool is_valid;


};