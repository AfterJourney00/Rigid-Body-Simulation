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
        this->last_dir = glm::vec3(0,1,0);
        this->Rotation_state = glm::mat4(1);
    }

    RigidBody(glm::vec3 Xt, double m, glm::vec3 rotate, float an) {
        this->Transformation = Xt;
        this->mass = m;
        this->Pt = glm::vec3(0.0f);
        this->angle = an;
        this->last_dir = rotate;
        this->Rotation_state = glm::rotate(glm::mat4(1), glm::radians(an), rotate);
    }

    glm::vec3 get_transformation () {
        return this->Transformation;
    }

    glm::vec3 get_rotation_dir () {
        if (this->Wt == glm::vec3(0)) {
            return glm::normalize(this->last_dir);
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

    void setIbody(float m, float r) {
        this->Ibody = m * powf(r,2) / 6.0f;		//正方体转动惯量的公式
        this->Iinv = 1 / Ibody;
    }

    void UpdateStates(glm::vec3 vt, glm::vec3 wt) {		//物体状态改变更新函数
        this->Vt = vt;					//更新 linear velocity
        this->Wt = wt;					//更新 angular velocity

        // record last dir
        if (wt != glm::vec3(0)) {
            this->last_dir = glm::normalize(wt);
        }

        this->Transformation += vt * time_interval;		//更新质心位移（位置改变）			time_interval还未定义
        this->angle = -glm::length(wt);		//更新物体旋转角度
        // 				time_interval还未定义
        //std::cout<<"angle: "<<angle<<std::endl;
        this->update_rotation();
    }

    glm::mat4 to_world() {
        glm::mat4 model = glm::mat4(1);
        model = glm::translate(model, this->get_transformation());
        model = glm::rotate(model, glm::radians(this->get_rotation_angle()), this->get_rotation_dir());
        model *= Rotation_state;
        // 更新rotate state
        // 不能写在这里，会导致不更新，未知原因
        //        this->Rotation_state = glm::rotate(this->Rotation_state, glm::radians(this->get_rotation_angle()), this->get_rotation_dir());
        if (model == Rotation_state) {
            std::cout<<"ERROR"<<std::endl;
        }
        return model;
    }

    void update_rotation() {
        this->Rotation_state = glm::rotate(this->Rotation_state, glm::radians(this->get_rotation_angle()), this->get_rotation_dir());
    }


    // һ��rigidbody��Ӧһ����ײ��
    std::vector<RigidBody*> possible_collision;
    int id;
private:

    // const quantities
    double mass;
    float Ibody;

    // sate variables
    glm::vec3 Transformation;
    glm::mat4 Rotation_state;


    glm::vec3 Pt;
    glm::vec3 Lt;

    // Derived quantities
    float Iinv;
    glm::vec3 Vt;
    glm::vec3 Wt;
    float angle;
    glm::vec3 last_dir;

    // computed quantities
    glm::vec3 Ft;
    glm::vec3 Tt;
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
