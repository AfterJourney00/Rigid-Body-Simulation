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
    RigidBody(glm::vec3 Xt,  double m) {
        this->Transformation = Xt;
        this->mass = m;
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

    glm::vec3 get_rotation_dir () {
        return glm::normalize(this->Wt);
    }

	float get_rotation_angle() {
		return this->angle;
	}

    glm::vec3 get_vt() {
        return Vt;
    }
    glm::vec3 get_wt() {
        return Wt;
    }

	glm::vec3 get_Pt() {
		return this->Pt;
	}

	glm::vec3 get_Lt() {
		return this->Lt;
	}

	double get_Ibody() {
		return this->Ibody;
	}

    float get_mass() {
        return mass;
    }

    void sum_Pt(glm::vec3 new_p) {
        this->Pt += new_p;
    }

    void sum_Lt(glm::vec3 new_l) {
        // 角动量变化值 = I−1(t) * τ impulse
        this->Lt += Iinv * new_l;
    }

	void setForce(glm::vec3 f) {
		this->Ft = f;
	}

	void setIbody(double mass, float r) {
		this->Ibody = mass * powf(r,2) / 6.0f;		//正方体转动惯量的公式
	}

	void UpdateStates(glm::vec3 vt, glm::vec3 wt) {		//物体状态改变更新函数
		this->Vt = vt;					//更新 linear velocity
		this->Wt = wt;					//更新 angular velocity

		this->Transformation = vt * time_interval;		//更新质心位移（位置改变）			time_interval还未定义
		this->angle = wt * time_interval;		//更新物体旋转角度				time_interval还未定义
	}

    glm::mat4 to_world() {
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, this->get_transformation());
        model = glm::rotate(model, glm::radians(this->get_rotation_angle()), this->get_rotation_dir());
        return model;
    }


	// һ��rigidbody��Ӧһ����ײ��
	std::vector<RigidBody> possible_collision;
private:
    //components
    std::vector<Particle> body_struct;

    // const quantities
    double mass;
    double Ibody;

    // sate variables
    glm::vec3 Transformation;
    glm::mat3 Rotation;

    glm::vec3 Pt;
    glm::vec3 Lt;

    // Derived quantities
    glm::mat3 Iinv;
    glm::vec3 Vt;
    glm::vec3 Wt;
	float angle;

    // computed quantities
    glm::vec3 Ft;
    glm::vec3 Tt;
};

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
