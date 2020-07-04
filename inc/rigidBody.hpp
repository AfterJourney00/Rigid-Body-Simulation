#pragma once

#include <utility>
#include <vector>
#include <omp.h>
#include <glm/glm.hpp>

class Particle {
public:
    Particle(glm::vec3 s_pos, glm::vec3 c_pos, float m){       //construct the particle with the distance and mass
        this->self_pos = s_pos;
        this->center_pos = c_pos;
        this->mass = m;
    }

    //compute the dis from the particle to the center
    float get_dis(){
        glm::vec3 coor = this->self_pos - this->center_pos;
        float dis = glm::length(coor);
        return dis;
    }

    //retrieve the position in body-space of the particle
    glm::vec3 get_pos(){
        return self_pos;
    }

    glm::vec3 get_force(){
        return force;
    }

    //compute the total force of this single particle
    void computeForce(std::vector<glm::vec3> f){
        glm::vec3 Force_tmp = glm::vec3 (0.0f, 0.0f, 0.0f);

        #pragma omp parallel for reduction(+:Force_tmp)
        for(auto i : f){
            Force_tmp += i;
        }
        this->force = Force_tmp;
    }

private:
    glm::vec3 self_pos{};
    glm::vec3 center_pos{};
    float mass;
    glm::vec3 force{};
};


class RigidBody {
public:
    RigidBody(glm::vec3 Xt, glm::mat4 Rt, double m, std::vector<Particle> particles) {
        this->Transformation = Xt;
        this->Rotation = Rt;
        this->mass = m;
        this->body_struct = std::move(particles);      //put the particles into the rigid-body
    }
    
    void UpdateForce(std::vector<glm::vec3> f) {
        glm::vec3 Force_tmp = glm::vec3 (0.0f, 0.0f, 0.0f);

        #pragma omp parallel for reduction(+:Force_tmp)
        for(int i = 0; i < f.size(); i++){
            Force_tmp += f[i];           //sum the total Linear Force
        }
        this->Ft = Force_tmp;       //update the linear force

        Force_tmp = glm::vec3 (0.0f, 0.0f, 0.0f);
        #pragma omp parallel for reduction(+:Force_tmp)
        for(auto & i : body_struct){
            Force_tmp += glm::cross(i.get_pos(), (i.get_force()));        //sum the total angular Force
        }
        this->Tt = Force_tmp;        //update the angular force
    }

    glm::vec3 get_transformation () {
        return this->Transformation;
    }

    glm::mat3 get_rotation () {
        return this->Rotation;
    }

private:
    //components
    std::vector<Particle> body_struct;

    // const quantities
    double mass;
    glm::mat3 Ibody;

    // sate variables
    glm::vec3 Transformation;
    glm::mat3 Rotation;

    glm::vec3 Pt;
    glm::vec3 Lt;

    // Derived quantities
    glm::mat3 Iinv;
    glm::vec3 Vt;
    glm::vec3 Wt;

    // computed quantities
    glm::vec3 Ft;
    glm::vec3 Tt;
};
