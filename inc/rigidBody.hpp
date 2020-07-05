#pragma once

#include <utility>
#include <vector>
#include <omp.h>
#include <glm/glm.hpp>

float time_interval = ((float) (1))/30;

class RigidBody {
public:
    RigidBody(glm::vec3 Xt,  double m) {
        this->Transformation = Xt;
        this->mass = m;
        this->Pt = glm::vec3(0.0f);
        this->Wt = glm::vec3(0.0f, 0.0f, 1.0f);
    }

    RigidBody(glm::vec3 Xt, double m, glm::vec3 rotate, float an) {
        this->Transformation = Xt;
        this->mass = m;
        this->Pt = glm::vec3(0.0f);
        this->Wt = rotate;
        this->angle = an;
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
        //std::cout<<new_p.z<<std::endl;
        //std::cout<<Pt.z<<std::endl;
        this->Pt += new_p;
        //std::cout<<Pt.z<<std::endl;
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

		this->Transformation += vt * time_interval;		//更新质心位移（位置改变）			time_interval还未定义
		this->angle += glm::length(wt) * time_interval;		//更新物体旋转角度				time_interval还未定义
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
