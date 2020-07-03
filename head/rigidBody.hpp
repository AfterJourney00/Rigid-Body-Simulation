#pragma once

#include <utility>
#include <vector>
#include <omp.h>

class Particle {
public:
    Particle(Eigen::Vector3f s_pos, Eigen::Vector3f c_pos, float m){       //construct the particle with the distance and mass
        this->self_pos = s_pos;
        this->center_pos = c_pos;
        this->mass = m;
    }

    //compute the dis from the particle to the center
    float get_dis(){                                                    
        Eigen::Vector3f coor = this->self_pos - this->center_pos;
        float dis = coor.norm();
        return dis;
    }

    //retrieve the position in body-space of the particle
    Eigen::Vector3f get_pos(){                            
        return self_pos;
    }

    Eigen::Vector3f get_force(){                            
        return force;
    }

    //compute the total force of this single particle
    void computeForce(std::vector<Eigen::Vector3f> f){
        Eigen::Vector3f Force_tmp = Eigen::Vector3f (0.0f, 0.0f, 0.0f);

        #pragma omp parallel for reduction(+:Force_tmp)
        for(int i=0; i<f.size(); i++){
            Force_tmp += f[i];
        }
        this->force = Force_tmp;
    }

private:
    Eigen::Vector3f self_pos;
    Eigen::Vector3f center_pos;
    float mass;
    Eigen::Vector3f force;
};
class RigidBody {
public:
    RigidBody(Eigen::Vector3f Xt, Eigen::Matrix3f Rt, double m, std::vector<Particle> particles) {
        this->Transformation = std::move(Xt);
        this->Rotation = std::move(Rt);
        this->mass = m;
        this->body_struct = particles;      //put the particles into the rigid-body
    }
    
    void UpdateForce(std::vector<Eigen::Vector3f> f) {
        Eigen::Vector3f Force_tmp = Eigen::Vector3f (0.0f, 0.0f, 0.0f);

        #pragma omp parallel for reduction(+:Force_tmp)
        for(int i = 0; i < f.size(); i++){
            Force_tmp += f[i];           //sum the total Linear Force
        }
        this->Ft = Force_tmp;       //update the linear force

        Force_tmp = Eigen::Vector3f (0.0f, 0.0f, 0.0f);
        #pragma omp parallel for reduction(+:Force_tmp)
        for(int i = 0; i < body_struct.size(); i++){
            Force_tmp += body_struct[i].get_pos().cross(body_struct[i].get_force());        //sum the total angular Force
        }
        this->Tt = Force_tmp;        //update the angular force
    }

private:
    //components
    std::vector<Particle> body_struct;

    // const quantities
    double mass;
    Eigen::Matrix3f Ibody;

    // sate variables
    Eigen::Vector3f Transformation;
    Eigen::Matrix3f Rotation;
    Eigen::Vector3f Pt;
    Eigen::Vector3f Lt;

    // Derived quantities
    Eigen::Matrix3f Iinv;
    Eigen::Vector3f Vt;
    Eigen::Vector3f Wt;

    // computed quantities
    Eigen::Vector3f Ft;
    Eigen::Vector3f Tt;
};
