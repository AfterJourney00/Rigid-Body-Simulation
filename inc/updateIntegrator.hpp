#pragma once

#include <utility>
#include <vector>
#include <omp.h>
#include <glm/glm.hpp>


const float GRAVITY = 10.0f;
const float MASS = 10.0f;
const float elasticity = 0.4;
const float slow = 0.06;

class updateIntegrator {
public:
    updateIntegrator (std::vector<RigidBody>* pos) {
        cubes = pos;
    };

    void update_cube_positions() {
        // todo: 计算所有的物体碰撞后改变的速度和角速度
        // todo: 按照输入的时间间隔移动cubes位置（修改RigidBody->transformation以及旋转

        // std::cout<<"new round"<<std::endl;
        // step1:
        for (int i = 0; i < (*cubes).size(); i++) {
            for (int j = i + 1; j < (*cubes).size(); j++) {
                float dis = glm::length((*cubes)[i].get_transformation() - (*cubes)[j].get_transformation());
                if (dis < sqrtf(3)) {
                    (*cubes)[i].possible_collision.push_back(&(*cubes)[j]);
                    (*cubes)[j].possible_collision.push_back(&(*cubes)[i]);
                }
            }
        }
        std::vector<Contact> contacts;
        // step2: 计算碰撞点
        for (auto &cube : (*cubes)) {
            while (!cube.possible_collision.empty()) {
                RigidBody *tem = cube.possible_collision[cube.possible_collision.size() - 1];
                Contact tem_contact = check_collision(*tem, cube);
                if (tem_contact.is_valid) {
                    contacts.push_back(tem_contact);
                } else if (glm::length(tem->get_transformation() - cube.get_transformation()) < 1.2f) {
                    Contact ball_volume;
                    ball_volume.a = tem;
                    ball_volume.b = &cube;
                    ball_volume.is_face_vertex = true;
                    ball_volume.is_valid = true;
                    ball_volume.face_normal = glm::normalize((tem->get_transformation() - cube.get_transformation()));
                    ball_volume.particle_position = 0.5f * (tem->get_transformation() - cube.get_transformation());

                    contacts.push_back(ball_volume);
                }
                int i = 0;
                for (i = 0; i < tem->possible_collision.size(); i++) {
                    if (tem->possible_collision[i]->get_transformation() == cube.get_transformation()) {
                        break;
                    }
                }
                tem->possible_collision.erase(tem->possible_collision.begin() + i);
                cube.possible_collision.pop_back();
            }
        }

        // step3: 计算动量变化
        for (auto &con : contacts) {
            process_collision(con);
        }

        // step4： 根据动量移动物体的transformation和Wt
        for (auto &cube: (*cubes)) {
            // 重力和地板
            //std::cout<<"before inside: "<<std::endl;
            //print_vec3(cubes[0].get_Pt());
            process_gravity_floor(cube);
            // 判定是否需要rest
            process_rest(cube);
            //std::cout<<"after inside1: "<<std::endl;
            //print_vec3(cubes[0].get_Pt());
            move_bodies(cube);
            //std::cout<<"after inside2: "<<std::endl;
            //print_vec3(cubes[0].get_Pt());
        }
    }

    void move_bodies(RigidBody &body) {
        // 根据动量角动量移动物体
        glm::vec3 curr_v = body.get_Pt() / MASS;
        //计算并更新此时刻物体的 linear velocity
        glm::vec3 curr_w;
        if (body.get_Lt() != glm::vec3(0)) {
            curr_w = body.get_Lt() / (float) body.get_Ibody();
        } else {
            curr_w = glm::vec3(0);
        }
        //计算并更新此时刻物体的 angular velocity
        body.UpdateStates(curr_v, curr_w);
    }

    void remove_noise(glm::vec3 &tao_impulse) {
        float accuracy = 0.01f;
        tao_impulse.x = tao_impulse.x < accuracy && tao_impulse.x > -accuracy ? 0 : tao_impulse.x;
        tao_impulse.y = tao_impulse.y < accuracy && tao_impulse.y > -accuracy ? 0 : tao_impulse.y;
        tao_impulse.z = tao_impulse.z < accuracy && tao_impulse.z > -accuracy ? 0 : tao_impulse.z;
    }

    void process_rest(RigidBody &body) {
        //print_vec3(body.get_Pt());
        //std::cout<< glm::length(body.get_Pt())<<std::endl;
        if (glm::length(body.get_Pt()) <= 0.1f) {
            body.reset_Pt();
        } else if (glm::length(body.get_Lt()) <= 0.01f) {
            //std::cout<<"rest Lt"<<std::endl;
            body.reset_Lt();
        } else {
            // 空气阻力等衰减
            body.sum_Pt(-0.005f*body.get_Pt());

            body.sum_Lt(-0.01f*body.get_Lt());
        }
    }

    static bool compare_height(glm::vec3 a,glm::vec3 b) {
        return a.y < b.y;
    }


    void process_gravity_floor(RigidBody &body) {
        if (body.get_transformation().y < sqrt(2)/2) {
            //std::cout<<"bounce"<<std::endl;
            // 物体接触了地面, 去除重力影响， 给重力一半的支持力
            // 初始化八个顶点body space
            std::vector<glm::vec3> points;
            glm::mat4 model = body.to_world();
            points.reserve(8);
            glm::vec3 tem = glm::vec4(0.5, 0.5, 0.5, 1) * model;
            points.push_back(tem);
            tem = glm::vec4(0.5, 0.5, -0.5, 1) * model;
            points.push_back(tem);
            tem = glm::vec4(0.5, -0.5, 0.5, 1) * model;
            points.push_back(tem);
            tem = glm::vec4(0.5, -0.5, -0.5, 1) * model;
            points.push_back(tem);
            tem = glm::vec4(-0.5, 0.5, 0.5, 1) * model;
            points.push_back(tem);
            tem = glm::vec4(-0.5, 0.5, -0.5, 1) * model;
            points.push_back(tem);
            tem = glm::vec4(-0.5, -0.5, 0.5, 1) * model;
            points.push_back(tem);
            tem = glm::vec4(-0.5, -0.5, -0.5, 1) * model;
            points.push_back(tem);

            double min_y = 10;
            std::vector<int> repeat_index;
            std::vector<glm::vec3> bottom;
            repeat_index.push_back(-1);
            for (int i = 0;i < points.size();i++) {
                bottom.push_back(points[i]);
                if (points[i].y < min_y) {
                    repeat_index.clear();
                    repeat_index.push_back(i);
                    min_y = points[i].y;
                } else if (points[i].y - min_y < 0.0001) {
                    repeat_index.push_back(i);
                }
            }

            // 稳定力
            std::sort(bottom.begin(), bottom.end(), compare_height);
            glm::vec3 bottom_center(0);
            for (int j = 0; j < 4; ++j) {
                bottom_center += bottom[j];
            }
            bottom_center /= 4.0f;
            glm::vec3 bottom_normal = glm::cross((bottom[0] - bottom[1]), (bottom[0] - bottom[2]));
            if (glm::dot((body.get_transformation() - bottom_center), bottom_normal) > 0) {
                bottom_normal *= -1.0f;
            }


            if (min_y <= 0.1) {
                // 一条边和地板相撞
                glm::vec3 hit_point;
                switch(repeat_index.size()) {
                    // 一个角撞击地面
                    case 1:
                        hit_point = points[repeat_index[0]];
                        break;
                        // 一条边撞击地面
                    case 2:
                        hit_point = 0.5f * (points[repeat_index[0]] + points[repeat_index[1]]);
                        break;
                        // 一个面撞击地面
                    case 4:
                        hit_point = 0.25f * (points[repeat_index[0]] + points[repeat_index[1]] + points[repeat_index[2]] + points[repeat_index[3]]);
                        break;
                }


                // 更新动量
                float v_rel = glm::dot(glm::vec3(0, -1.0f, 0.0f), body.get_vt());
                glm::vec3 delt_v = glm::vec3(0, 1.0f, 0.0f) * v_rel - body.get_vt();
                delt_v *= elasticity;
                glm::vec3 J = delt_v * body.get_mass();
                if (body.get_Pt().y < 0) {
                    J += glm::vec3(0, -body.get_Pt().y, 0);
                    //std::cout<<"J: "<<std::endl;
                    //print_vec3(J);
                    body.sum_Pt(J);
                }


                J *= 0.1f;

                delt_v = glm::vec3(0, -1.0f, 0.0f);
                delt_v *= GRAVITY * time_interval * slow;
                glm::vec3 J_gravity = delt_v * body.get_mass();
                J += J_gravity * 0.1f;

                glm::vec3 tao_steady(0);

                tao_steady = -glm::cross(glm::normalize(bottom_normal), glm::vec3(0,1,0)) * MASS;

                // 更新角动量
                glm::vec3 re_xt = hit_point - body.get_transformation();
                //remove_noise(re_xt);
                glm::vec3 tao_impulse = glm::cross(J, re_xt) * 0.5f;
                //std::cout<<"re_xt"<<std::endl;
                //print_vec3(re_xt);
                //std::cout<<"J"<<std::endl;
                //print_vec3(J);

                remove_noise(tao_impulse);

                //std::cout<<glm::dot(glm::normalize(bottom_normal), glm::vec3(0,1,0))<<std::endl;


                if (glm::length(body.get_Pt()) < 20.0f || glm::length(body.get_Lt()) < 10.0f) {
                    if (glm::length(body.get_Pt()) < 6.0f && glm::dot(glm::normalize(bottom_normal), glm::vec3(0,1,0)) < 0.95) {
                        // 慢速归位
                        if (glm::length(body.get_Pt()) < 2.0f && glm::dot(glm::normalize(bottom_normal), glm::vec3(0,1,0)) < 0.95) {
                            body.sum_Lt(body.get_Lt() * -0.9f);
                            body.sum_Lt(1.0f * tao_steady);
                            //std::cout<<"tao after"<<std::endl;
                            //print_vec3(body.get_Lt());
                        } else {
                            body.sum_Lt(body.get_Lt() * -0.9f);
                            //std::cout<<"tao_steady"<<std::endl;
                            //print_vec3(0.4f * tao_steady);
                            body.sum_Lt(0.4f * tao_steady);
                            //std::cout<<"tao after"<<std::endl;
                            //print_vec3(body.get_Lt());
                        }
                    } else {
                        // 减速过程
                        body.sum_Lt(body.get_Lt() * -0.8f);
                        //std::cout<<"tao_steady"<<std::endl;
                        //print_vec3(0.5f * tao_steady + 0.5f * tao_impulse);

                        // 刹车过程
                        body.sum_Lt(0.4f * tao_steady);
                    }

                } else {
                    tao_impulse += tao_steady;
                    //std::cout<<"tao_impulse"<<std::endl;
                    //print_vec3(tao_impulse);
                    body.sum_Lt(tao_impulse * 0.4f);
                }
                //body.sum_Lt(body.get_Lt() * -0.9f);
                //body.sum_Lt( tao_steady);

                return ;
            }


        }
        if (body.get_transformation().y >= sqrt(2)/2) {
            // 物体没有接触地面
            // std::cout<<"gravity"<<std::endl;
            glm::vec3 delt_v = glm::vec3(0, -1.0f, 0.0f);
            delt_v *= GRAVITY * time_interval * slow;
            glm::vec3 J = delt_v * body.get_mass();
            //std::cout<<"before: "<<std::endl;
            //print_vec3(body.get_Pt());
            body.sum_Pt(J);
            //std::cout<<"after: "<<std::endl;
            //print_vec3(body.get_Pt());
        }

    }

    glm::vec4 face_function(glm::vec3 p, glm::vec3 n) {
        return glm::vec4(n, glm::dot(p, n));
    }

    Contact check_collision(RigidBody& a, RigidBody& b) {
        Contact f_v;
        f_v.a = &a;
        f_v.b = &b;
        // 面的法向量
        std::vector<glm::vec3> face_normal;
        // 六个面分别的四个点
        std::vector<std::vector<glm::vec3>> face_point;
        std::vector<glm::vec3> tem;

        tem.emplace_back(0.5f, 0.5f, -0.5f);
        tem.emplace_back(0.5f, -0.5f, -0.5f);
        tem.emplace_back(0.5f, -0.5f, 0.5f);
        tem.emplace_back(0.5f, 0.5f, 0.5f);
        face_point.push_back(tem);
        tem.clear();

        face_normal.emplace_back(1.0f, 0.0f, 0.0f);

        tem.emplace_back(-0.5f, 0.5f, -0.5f);
        tem.emplace_back(-0.5f, -0.5f, -0.5f);
        tem.emplace_back(-0.5f, -0.5f, 0.5f);
        tem.emplace_back(-0.5f, 0.5f, 0.5f);

        face_point.push_back(tem);
        tem.clear();

        face_normal.emplace_back(-1.0f, 0.0f, 0.0f);

        tem.emplace_back(0.5f, 0.5f, -0.5f);
        tem.emplace_back(-0.5f, 0.5f, -0.5f);
        tem.emplace_back(-0.5f, 0.5f, 0.5f);
        tem.emplace_back(0.5f, 0.5f, 0.5f);

        face_point.push_back(tem);
        tem.clear();
        face_normal.emplace_back(0.0f, 1.0f, 0.0f);

        tem.emplace_back(0.5f, -0.5f, -0.5f);
        tem.emplace_back(-0.5f, -0.5f, -0.5f);
        tem.emplace_back(-0.5f, -0.5f, 0.5f);
        tem.emplace_back(0.5f, -0.5f, 0.5f);

        face_point.push_back(tem);
        tem.clear();
        face_normal.emplace_back(0.0f, -1.0f, 0.0f);

        tem.emplace_back(0.5f, 0.5f, 0.5f);
        tem.emplace_back(0.5f, -0.5f, 0.5f);
        tem.emplace_back(-0.5f, -0.5f, 0.5f);
        tem.emplace_back(-0.5f, 0.5f, 0.5f);

        face_point.push_back(tem);
        tem.clear();
        face_normal.emplace_back(0.0f, 0.0f, 1.0f);


        tem.emplace_back(0.5f, 0.5f, -0.5f);
        tem.emplace_back(0.5f, -0.5f, -0.5f);
        tem.emplace_back(-0.5f, -0.5f, -0.5f);
        tem.emplace_back(-0.5f, 0.5f, -0.5f);

        face_point.push_back(tem);
        tem.clear();

        face_normal.emplace_back(0.0f, 0.0f, -1.0f);

        // face information including the face vertices(0-3) and normal(4)
        std::vector<std::vector<glm::vec3>> face_p_n_a;
        std::vector<std::vector<glm::vec3>> face_p_n_b;
        std::vector<glm::vec3> temp;
        // 初始化两个vector
        for (int j = 0; j < 5; j++) {
            temp.emplace_back(0.0f, 0.0f, 0.0f);
        }
        for (int i = 0; i < 6; i++) {
            face_p_n_a.push_back(temp);
            face_p_n_b.push_back(temp);
        }

        // 面上的点和法向量转换成世界坐标系
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 4; j++) {
                face_p_n_a[i][j] = a.to_world() * glm::vec4(face_point[i][j], 1.0f);
                face_p_n_b[i][j] = b.to_world() * glm::vec4(face_point[i][j], 1.0f);
            }
            face_p_n_a[i][4] = a.to_world() * glm::vec4(face_normal[i], 1.0f);
            face_p_n_b[i][4] = a.to_world() * glm::vec4(face_normal[i], 1.0f);
        }
        glm::vec4 func_1f, func_2f;
        // line[face_index][0-2n] for n lines each face
        Line line;
        std::vector<Segment> line_possible;
        std::vector<std::vector<int>> face_line_number_a;
        std::vector<std::vector<int>> face_line_number_b;
        std::vector<int> temp_int;
        temp_int.clear();
        face_line_number_a.reserve(6);
        face_line_number_b.reserve(6);
        for (int i = 0; i < 6; i++) {
            face_line_number_a.push_back(temp_int);
        }
        for (int i = 0; i < 6; i++) {
            face_line_number_b.push_back(temp_int);
        }
        for (int i = 0; i < 6; i++) {
            for (int j = i; j < 6; j++) {
                // 传入face_function面上的一个点和一个法向量
                // 计算出面的表达式ax+by+cz=d 结果存在一个vec4里面
                func_1f = face_function(face_p_n_a[i][0], face_p_n_a[i][4]);	//计算面a表达式
                func_2f = face_function(face_p_n_b[j][0], face_p_n_b[j][4]);	//计算面b表达式
                // result[0] 一个点 result[1] 方向   直线表达形式(a,b,c) + t * (d, e, f)
                line = calculate_line(func_1f, func_2f, face_p_n_a[i][4], face_p_n_b[j][4]);
                Segment line_segment;
                int state = judge_line_possibility(face_p_n_a[i], face_p_n_b[j], line, line_segment);
                if (state == -1) {
                    f_v.is_valid = false;
                    return f_v;
                }
                if (state) {
                    line_possible.push_back(line_segment);
                    face_line_number_a[i].push_back(line_possible.size() - 1);
                    face_line_number_b[j].push_back(line_possible.size() - 1);
                }
            }
        }
        // the num of face a, b needed to be targeted
        int num_of_face_a = 0;
        int num_of_face_b = 0;
        // the face index targeted
        std::vector<int> target_face_a;
        std::vector<int> target_face_b;
        for (int i = 0; i < 6; i++) {
            if (!face_line_number_a[i].empty()) {
                target_face_a.push_back(i);
                num_of_face_a += 1;
            }
        }

        for (int i = 0; i < 6; i++) {
            if (!face_line_number_b[i].empty()) {
                target_face_b.push_back(i);
                num_of_face_b += 1;
            }
        }

        f_v.is_valid = false;
        if (num_of_face_a == 1) {
            // 单面相交情况
            Segment current_lines = line_possible[target_face_a[0]];
            f_v.is_valid = true;
            f_v.is_face_vertex = true;
            f_v.face_normal = face_normal[target_face_a[0]];
            if (line_possible.size() == 3) {
                // 点面相交
                // 输入三个线段， 返回一个顶点
                f_v.particle_position = calculate_face_vertex(line_possible[face_line_number_a[target_face_a[0]][0]],line_possible[face_line_number_a[target_face_a[0]][1]], line_possible[face_line_number_a[target_face_a[0]][2]]);
            }
        } else if (num_of_face_b == 1) {
            Segment current_lines = line_possible[target_face_b[0]];
            f_v.is_valid = true;
            f_v.is_face_vertex = true;
            f_v.face_normal = face_normal[target_face_b[0]];
            if (line_possible.size() == 3) {
                // 点面相交
                // 输入三个线段， 返回一个顶点
                f_v.particle_position = calculate_face_vertex(line_possible[face_line_number_b[target_face_b[0]][0]],line_possible[face_line_number_b[target_face_b[0]][1]], line_possible[face_line_number_b[target_face_b[0]][2]]);
            }
        } else if (num_of_face_b == 2) {
            if (num_of_face_a == 2) {
                f_v.is_valid = true;
                f_v.is_face_vertex = false;
                // 线线相交
                // 输入四条线段 返回两条棱
                if (face_line_number_a[target_face_a[0]].size() != 2 || face_line_number_a[target_face_a[1]].size() != 2) {
                    f_v.is_valid = false;
                    return f_v;
                }
                int i_1 = face_line_number_a[target_face_a[0]][0], i_2 = face_line_number_a[target_face_a[0]][1], i_3 = face_line_number_a[target_face_a[1]][0], i_4 = face_line_number_a[target_face_a[1]][1];
                std::vector<glm::vec3> edges = calculate_edge(line_possible[i_1],line_possible[i_2],line_possible[i_3],line_possible[i_4]);
                f_v.edge1 = edges[0];
                f_v.edge1 = edges[1];
            } else {
                f_v.is_valid = false;
                printf("error\n");
                return f_v;
            }
        }

        return f_v;
    }

    std::vector<glm::vec3> calculate_edge(Segment line1, Segment line2, Segment line3, Segment line4) {
        // 两个三角形中共六个点，其中在一个棱上有两组重复点
        // 两组重复点形成一条edge  一队single 点形成一条edge
        glm::vec3 s1 = line1.start, s2;
        glm::vec3 e1 = line1.end, e2;

        if (s1 == line2.start) {
            e2 = line2.end;
        } else if (s1 == line2.end) {
            e2 = line2.start;
        } else if (s1 == line3.start) {
            e2 = line3.end;
        } else if (s1 == line3.end) {
            e2 = line3.start;
        } else if (s1 == line4.start) {
            e2 = line4.end;
        } else if (s1 == line4.end) {
            e2 = line4.start;
        }

        if (e1 == line2.start) {
            s2 = line2.end;
        } else if (e1 == line2.end) {
            s2 = line2.start;
        } else if (e1 == line3.start) {
            s2 = line3.end;
        } else if (e1 == line3.end) {
            s2 = line3.start;
        }  else if (e1 == line4.start) {
            s2 = line4.end;
        } else if (e1 == line4.end) {
            s2 = line4.start;
        }

        std::vector<glm::vec3> result;
        result.push_back(s1 - s2);
        result.push_back(e1 - e2);

        return result;
    }

    glm::vec3 calculate_face_vertex(Segment line1, Segment line2, Segment line3) {
        return ((line1.end + line1.start) + (line2.end + line2.start) +(line3.end + line3.start)) / 6.0f;
    }

    int judge_line_possibility(std::vector<glm::vec3> face_a, std::vector<glm::vec3> face_b, Line line, Segment line_segment) {
        std::vector<Line> face_a_lines;
        std::vector<Line> face_b_lines;
        std::vector<Segment> face_a_segs;
        std::vector<Segment> face_b_segs;
        for (int i = 1; i < 4; i++) {			//将面a和面b的四条边各自存放到face_a_lines和face_b_lines中		以及以线段形式存放在face_a_segs和face_b_segs
            Line curr_a_line(face_a[i - 1], face_a[i] - face_a[i - 1]);		//Line这个类的dir变量或许需要取normalize
            Line curr_b_line(face_b[i - 1], face_b[i] - face_b[i - 1]);
            Segment curr_a_seg(face_a[i - 1], face_a[i]);
            Segment curr_b_seg(face_b[i - 1], face_b[i]);

            face_a_lines.push_back(curr_a_line);
            face_b_lines.push_back(curr_b_line);
            face_a_segs.push_back(curr_a_seg);
            face_b_segs.push_back(curr_b_seg);
        }

        Line curr_a_line(face_a[0], face_a[3] - face_a[0]);
        Line curr_b_line(face_b[0], face_b[3] - face_b[0]);
        Segment curr_a_seg(face_a[0], face_a[3]);
        Segment curr_b_seg(face_b[0], face_b[3]);

        face_a_lines.push_back(curr_a_line);
        face_b_lines.push_back(curr_b_line);
        face_a_segs.push_back(curr_a_seg);
        face_b_segs.push_back(curr_b_seg);


        //std::vector<glm::vec3> intersection_vec_a;
        //std::vector<glm::vec3> intersection_vec_b;
        std::vector<float> intersection_vec_a;
        std::vector<float> intersection_vec_b;
        for (int i = 0; i < face_a_lines.size(); i++) {
            glm::vec3 intersection_a = solution_lines_vertex(face_a_lines[i], line);		//计算交线与面a第i条边的交点坐标
            glm::vec3 intersection_b = solution_lines_vertex(face_b_lines[i], line);		//计算交线与面b第i条边的交点坐标
            if (intersection_a == glm::vec3(-1,-1,-1) || intersection_b == glm::vec3(-1,-1,-1)) {
                return -1;
            }

            glm::vec3 check_line1_a = face_a_segs[i].start - intersection_a;			//得到从交点到面a第i条边第一个顶点的向量
            glm::vec3 check_line2_a = face_a_segs[i].end - intersection_a;			//得到从交点到面a第i条边第二个顶点的向量

            glm::vec3 check_line1_b = face_b_segs[i].start - intersection_b;			//得到从交点到面b第i条边第一个顶点的向量
            glm::vec3 check_line2_b = face_b_segs[i].end - intersection_b;			//得到从交点到面b第i条边第二个顶点的向量

            if (!Is_Same_dir(check_line1_a, check_line2_a) && !Is_Same_dir(check_line1_b, check_line2_b)) {	//两个向量反向，说明交点intersection在线段两端中间 是一个有效的交点
                float t1 = (intersection_a - line.ori).x / line.dir.x;			//计算交点对应line方程的参数t1
                float t2 = (intersection_b - line.ori).x / line.dir.x;			//计算交点对应line方程的参数t2
                //float t2 = (intersection_a - line.ori).y() / line.dir.y();
                //float t3 = (intersection_a - line.ori).z() / line.dir.z();
                intersection_vec_a.push_back(t1);
                intersection_vec_b.push_back(t2);
                if (t1 == NAN || t2 == NAN) return -1;						//遇到nan, 返回-1
            }
        }
        if (intersection_vec_a.empty() || intersection_vec_b.empty()) return false;
        else {
            return SelectSegment(intersection_vec_a, intersection_vec_b, line, line_segment);
        }
    }

    bool SelectSegment(std::vector<float> v1, std::vector<float> v2, Line line, Segment& line_seg) {		//如果存在 就一定是有两个交点？（有待确认）
        float t_min = fmax(fmin(v1[0], v1[1]), fmin(v2[0], v2[1]));		//取交集[t_min, t_max]中的t_min
        float t_max = fmin(fmax(v1[0], v1[1]), fmax(v2[0], v2[1]));		//取交集[t_min, t_max]中的t_max

        if (t_min > t_max) return false;			//如果t_min > t_max, 说明disjoint， 不设置line segment
        else {								//如果t_min <= t_max, 说明交集存在， 设置line segment
            line_seg.start = line.ori + t_min * line.dir;
            line_seg.end = line.ori + t_max * line.dir;
            return true;
        }
    }

    bool Is_Same_dir(glm::vec3 vector1, glm::vec3 vector2) {
        float dx = vector1.x / vector2.x;
        //float dy = vector1.y() / vector2.y();
        //float dz = vector1.z() / vector2.z();
        if (dx > 0) return true;		//比值大于0，同向
        else return false;				//比值小于0，反向
    }

    glm::vec3 solution_lines_vertex(Line line1, Line line2) {
        float x;
        x = solution_of_functions(glm::vec3(line1.dir.x, -line2.dir.x, (line2.ori.x - line1.ori.x)), glm::vec3(line1.dir.y, -line2.dir.y, (line2.ori.y - line1.ori.y)));
        if (x == -1) {
            return { -1,-1,-1};
        }
        return (line1.dir * x + line1.ori);
    }

    float solution_of_functions(glm::vec3 func1, glm::vec3 func2) {
        if (func2.y != 0) {
            if ((func1.x - func1.y / func2.y * func2.x) != 0) {
                return (func1.z - func1.y / func2.y * func2.z) / (func1.x - func1.y / func2.y * func2.x);
            }
            else {
                printf("error1!\n");
                return -1;
            }
        }
        else {
            if (func2.x != 0) {
                return func2.z / func2.x;
            }
            else {
                printf("error2!\n");
                return -1;
            }
        }
    }

    Line calculate_line(glm::vec4 func_1f, glm::vec4 func_2f, glm::vec3 n1, glm::vec3 n2) {

        float x, y, z;
        if (func_1f.y * func_2f.z != func_2f.y * func_1f.z) {
            if (func_2f.z != 0.0f) {
                x = 0.0f;
                y = (func_1f.w - (func_1f.z / func_2f.z) * func_2f.w) / (func_1f.y - (func_1f.z / func_2f.z) * func_2f.y);
                z = (func_2f.w - y * func_2f.y) / func_2f.z;
            }
            else {
                x = 0.0f;
                y = func_2f.w / func_2f.y;
                z = (func_1f.w - y * func_1f.y) / func_1f.z;
            }
        }
        else {
            if (func_2f.y != 0.0f) {
                z = 0.0f;
                x = (func_1f.w - (func_1f.y / func_2f.y) * func_2f.w) / (func_1f.x - (func_1f.y / func_2f.y) * func_2f.x);
                y = (func_2f.w - func_2f.x * x) / func_2f.y;
            }
            else {
                z = 0.0f;
                x = func_2f.w / func_2f.x;
                y = (func_1f.w - func_1f.x * x) / func_1f.y;
            }
        }

        return {glm::vec3(x, y, z), glm::cross(n1, n2)};
    }

    void process_collision(Contact con) {
        if (con.is_valid) {
            glm::vec3 normal;
            if (con.is_face_vertex) {
                normal = glm::normalize(con.face_normal);
            } else {
                normal = glm::normalize(glm::cross(con.edge1, con.edge2));
            }
            glm::vec3 pat = con.a->get_vt() + con.a->get_wt() * (con.particle_position - con.a->get_transformation());
            glm::vec3 pbt= con.b->get_vt() + con.b->get_wt() * (con.particle_position - con.b->get_transformation());
            float v_rel;
            v_rel = glm::dot(normal, (pat - pbt))* elasticity;
            if (v_rel >= 0) {
                // 物体正在远离不需要处理碰撞了
                return ;
            }
            glm::vec3 Ja, Jb;
            glm::vec3 det_va = v_rel * -normal - con.a->get_vt();
            glm::vec3 det_vb = v_rel * normal - con.b->get_vt();
            Ja = det_va * con.a->get_mass();
            Jb = det_vb * con.a->get_mass();
            glm::vec3 tao_a_impulse = glm::cross((con.particle_position - con.a->get_transformation()), Ja);
            glm::vec3 tao_b_impulse = glm::cross((con.particle_position - con.b->get_transformation()), Jb);
            // 更新body a和b的动量和角动量

            remove_noise(tao_a_impulse);
            remove_noise(tao_b_impulse);
            remove_noise(Ja);
            remove_noise(Jb);

            std::vector<glm::vec3> points_a;
            glm::mat4 model_a = con.a->to_world();
            points_a.reserve(8);
            glm::vec3 tem = glm::vec4(0.5, 0.5, 0.5, 1) * model_a;
            points_a.push_back(tem);
            tem = glm::vec4(0.5, 0.5, -0.5, 1) * model_a;
            points_a.push_back(tem);
            tem = glm::vec4(0.5, -0.5, 0.5, 1) * model_a;
            points_a.push_back(tem);
            tem = glm::vec4(0.5, -0.5, -0.5, 1) * model_a;
            points_a.push_back(tem);
            tem = glm::vec4(-0.5, 0.5, 0.5, 1) * model_a;
            points_a.push_back(tem);
            tem = glm::vec4(-0.5, 0.5, -0.5, 1) * model_a;
            points_a.push_back(tem);
            tem = glm::vec4(-0.5, -0.5, 0.5, 1) * model_a;
            points_a.push_back(tem);
            tem = glm::vec4(-0.5, -0.5, -0.5, 1) * model_a;
            points_a.push_back(tem);
            //cube b
            std::vector<glm::vec3> points_b;
            glm::mat4 model_b = con.b->to_world();
            points_b.reserve(8);
            tem = glm::vec4(0.5, 0.5, 0.5, 1) * model_b;
            points_b.push_back(tem);
            tem = glm::vec4(0.5, 0.5, -0.5, 1) * model_b;
            points_b.push_back(tem);
            tem = glm::vec4(0.5, -0.5, 0.5, 1) * model_b;
            points_b.push_back(tem);
            tem = glm::vec4(0.5, -0.5, -0.5, 1) * model_b;
            points_b.push_back(tem);
            tem = glm::vec4(-0.5, 0.5, 0.5, 1) * model_b;
            points_b.push_back(tem);
            tem = glm::vec4(-0.5, 0.5, -0.5, 1) * model_b;
            points_b.push_back(tem);
            tem = glm::vec4(-0.5, -0.5, 0.5, 1) * model_b;
            points_b.push_back(tem);
            tem = glm::vec4(-0.5, -0.5, -0.5, 1) * model_b;
            points_b.push_back(tem);

            std::vector<glm::vec3> bottom_a, bottom_b;
            for (int i = 0; i < points_a.size(); i++) {
                bottom_a.push_back(points_a[i]);
                bottom_b.push_back(points_b[i]);
            }

            // 稳定力
            std::sort(bottom_a.begin(), bottom_a.end(), compare_height);
            std::sort(bottom_b.begin(), bottom_b.end(), compare_height);

            glm::vec3 n_a, n_b;
            n_a = glm::cross((bottom_a[0] - bottom_a[1]), (bottom_a[1] - bottom_a[2]));
            n_b = glm::cross((bottom_b[0] - bottom_b[1]), (bottom_b[1] - bottom_b[2]));

            n_a = glm::normalize(n_a);
            n_b = glm::normalize(n_b);

            float x = glm::length(n_a - n_b);
            if (glm::length(n_a + n_b) < x) {
                x = glm::length(n_a + n_b);
            }

            /*

            glm::vec3 bottom_center_a(0);
            glm::vec3 bottom_center_b(0);
            for (int j = 0; j < 4; ++j) {
                bottom_center_a += bottom_a[j];
            }

            for (int j = 0; j < 4; ++j) {
                bottom_center_b += bottom_b[j];
            }
            bottom_center_a /= 4.0f;
            bottom_center_b /= 4.0f;
            */
            float height = glm::dot(con.a->get_transformation() - con.b->get_transformation(), glm::vec3(0,1,0));
            float dis = glm::length(con.a->get_transformation() - con.b->get_transformation());
            height = height < 0 ? - height : height;
            if (x < 0.1f && glm::length(con.a->get_Pt()) < 5.0f && glm::length(con.b->get_Pt()) < 5.0f && (height <= 1.0f && dis < 1.2f)) {
                glm::vec3 delt_v = glm::vec3(0, 1.0f, 0.0f);
                delt_v *= GRAVITY * time_interval * slow;
                glm::vec3 J = delt_v * MASS;
                con.a->reset_Pt();
                con.b->reset_Pt();
                //con.a->sum_Pt(8.0f*J);
                //con.b->sum_Pt(8.0f*J);
                tao_a_impulse = tao_a_impulse.y * glm::vec3(0, 0.1f, 0);
                tao_b_impulse = tao_b_impulse.y * glm::vec3(0, 0.1f, 0);

                //tao_a_impulse = glm::vec3(0.0f, 0.0f, 0.0f);
                //tao_b_impulse = glm::vec3(0.0f, 0.0f, 0.0f);
            } else {
                con.a->sum_Pt(Ja);
                con.b->sum_Pt(Jb);
            }


            con.a->sum_Lt(tao_a_impulse * 0.1f);
            con.b->sum_Lt(tao_b_impulse* 0.1f);
        }
    }

    void check_calculate_line() {

        std::cout<<"unit test calculate_line"<<std::endl;
        glm::vec3 a(0.5,0.5,0.5);
        glm::vec3 a_n(0,0,1);
        glm::vec3 b(0.5,0.25,0.5);
        glm::vec3 b_n(0,-1,1);

        glm::vec4 func_1f = face_function(a, a_n);
        glm::vec4 a_t(0.25,0.5,0.5, 0);
        // 检查点是否在面上
        if (glm::dot(func_1f, a_t) != func_1f.w){
            std::cout<<"face_function error!"<<std::endl;
        }

        glm::vec4 func_2f = face_function(b,b_n);
        glm::vec4 b_t(0.3,0.45,0.7, 0);
        // 检查点是否在面上
        if (glm::dot(func_2f, b_t) != func_2f.w){
            std::cout<<"face_function error!"<<std::endl;
        }
        Line result = calculate_line(func_1f, func_2f, a_n, b_n);
        glm::vec3 test_p(1.5,0.25,0.5);
        if (glm::normalize(test_p - result.ori) != glm::normalize(result.dir)
            || glm::normalize(b - result.ori) != glm::normalize(result.dir)) {
            std::cout<<"calculate_line error!"<<std::endl;
        }

        Segment seg1(glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, -2.0f, 1.0f));
        Segment seg2(glm::vec3(0.0f, 0.0f, 2.0f), glm::vec3(0.0f, -2.0f, 1.0f));
        Segment seg3(glm::vec3(-2.0f, 0.0f, 1.0f), glm::vec3(0.0f, 0.0f, 2.0f));
        Segment seg4(glm::vec3(-2.0f, 0.0f, 1.0f), glm::vec3(0.0f, 0.0f, 0.0f));
        std::vector<glm::vec3> vec;
        vec = calculate_edge(seg1, seg2, seg3, seg4);

        std::vector<glm::vec3> a_points;
        a_points.emplace_back(0.5,0.5,0.5);
        a_points.emplace_back(0.5,-0.5,0.5);
        a_points.emplace_back(-0.5,0.5,0.5);
        a_points.emplace_back(-0.5,-0.5,0.5);
        a_points.emplace_back(a_n);

        std::vector<glm::vec3> b_points;
        b_points.emplace_back(0.5,0.25,0.5);
        b_points.emplace_back(-0.5,0.25,0.5);
        b_points.emplace_back(0.5,0.45,0.7);
        b_points.emplace_back(-0.5,0.45,0.7);
        b_points.emplace_back(b_n);

        Segment line_segment;

        // 线过两个面
        if(!judge_line_possibility(a_points, b_points, result, line_segment)
           && (line_segment.start != glm::vec3(0.5,0.25,0.5) && line_segment.start != glm::vec3(-0.5,0.25,0.5))
           && (line_segment.end != glm::vec3(0.5,0.25,0.5) && line_segment.end != glm::vec3(-0.5,0.25,0.5))) {
            std::cout<<"judge error 1"<<std::endl;
        }


        a_points.clear();
        a_points.emplace_back(0.5,0,0.5);
        a_points.emplace_back(0.5,-0.5,0.5);
        a_points.emplace_back(-0.5,0,0.5);
        a_points.emplace_back(-0.5,-0.5,0.5);
        a_points.emplace_back(a_n);

        b_points.clear();
        b_points.emplace_back(0.5,0.25,0.5);
        b_points.emplace_back(-0.5,0.25,0.5);
        b_points.emplace_back(0.5,0.45,0.7);
        b_points.emplace_back(-0.5,0.45,0.7);
        b_points.emplace_back(b_n);

        // 线不过一个面
        if(judge_line_possibility(a_points, b_points, result, line_segment)) {
            std::cout<<"judge error 2"<<std::endl;
        }

        a_points.clear();
        a_points.emplace_back(0.5,0,0.5);
        a_points.emplace_back(0.5,-0.5,0.5);
        a_points.emplace_back(-0.5,0,0.5);
        a_points.emplace_back(-0.5,-0.5,0.5);
        a_points.emplace_back(a_n);

        b_points.clear();
        b_points.emplace_back(0.5,0.35,0.6);
        b_points.emplace_back(-0.5,0.35,0.6);
        b_points.emplace_back(0.5,0.45,0.7);
        b_points.emplace_back(-0.5,0.45,0.7);
        b_points.emplace_back(b_n);

        // 线不过两个面
        if(judge_line_possibility(a_points, b_points, result, line_segment)) {
            std::cout<<"judge error 3"<<std::endl;
        }
    }

private:
    std::vector<RigidBody> *cubes;
};