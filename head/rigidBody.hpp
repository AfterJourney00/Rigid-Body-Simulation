#pragma once

#include <utility>
#include <vector>

class RigidBody {
public:
    RigidBody(Eigen::Vector3f Xt, Eigen::Matrix3f Rt, double m) {
        this->Transformation = std::move(Xt);
        this->Rotation = std::move(Rt);
        this->mass = m;
    }

private:
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
