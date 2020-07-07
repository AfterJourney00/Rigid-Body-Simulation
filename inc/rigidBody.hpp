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
        this->Lt = glm::vec3(0.0f);
        this->angle = 0;
    }

    RigidBody(glm::vec3 Xt, double m, glm::vec3 rotate, float an) {
        this->Transformation = Xt;
        this->mass = m;
        this->Pt = glm::vec3(0.0f);
        this->Lt = rotate;
        this->angle = an;
    }

    glm::vec3 get_transformation () {
        return this->Transformation;
    }

    glm::vec3 get_rotation_dir () {
        if (this->Wt == glm::vec3(0)){
            return glm::vec3(0);
        }
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

    void reset_Pt() {
        this->Pt = glm::vec3(0);
    }

    void reset_Lt() {
        this->Lt = glm::vec3(0);
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
        this->Iinv = 1 / Ibody;
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
        if (glm::dot(this->get_rotation_dir(), glm::vec3(1)) != 0) {
            model = glm::rotate(model, glm::radians(this->get_rotation_angle()), this->get_rotation_dir());
        }
        return model;
    }


	// һ��rigidbody��Ӧһ����ײ��
	std::vector<RigidBody*> possible_collision;
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
    float Iinv;
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

class Line {
public:
	Line(glm::vec3 o, glm::vec3 d) {
		dir = glm::normalize(d);
		ori = o;
	};
	Line() {};
	glm::vec3 dir;
	glm::vec3 ori;
};

class Segment {
public:
	Segment(glm::vec3 s, glm::vec3 e) {
		start = s;
		end = e;
	};
	Segment(){};
	glm::vec3 start;
	glm::vec3 end;
};
