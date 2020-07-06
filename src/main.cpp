
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include "glm/gtc/matrix_transform.hpp"
#include <iostream>
#include <vector>
#include <inc/camera.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <inc/rigidBody.hpp>

#include "../inc/my_texture.h"
#include "../inc/shader_m.h"
#include "tiny_obj_loader.h"
#include <math.h>



/*-----------------------------------------------------------------------*/
//Here are some mouse and keyboard function. You can change that.
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow* window);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);

const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
float degreeX = (360 * lastX / 400);
float degreeY = (360 * lastY / 300);
bool firstMouse = true;
float deltaTime = 0.0f;	// time between current frame and last frame
float lastFrame = 0.0f;
float OX = 0;//should be update to a new coordinate
float OY = 0;
float OZ = 0;
float currentFrame;

#define MASS 10
#define GRAVITY 10
const float slow = 0.01;

Camera *camera;

void print_vec3(glm::vec3 cube) {
    std::cout<<cube.x<<std::endl;
    std::cout<<cube.y<<std::endl;
    std::cout<<cube.z<<std::endl;
}

void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera->ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera->ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera->ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera->ProcessKeyboard(RIGHT, deltaTime);
}
// 处理鼠标移动
void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera->ProcessMouseMovement(xoffset, yoffset);
}


void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void initPMV(Shader my_shader, Shader light_shader, const glm::vec3 pointLightPositions[]) {
    camera = new Camera (glm::vec3(0.0f, 4.0f, 3.0f));

    glm::mat4 projection = glm::perspective(glm::radians(45.0f), 800.0f / 600.0f, 0.1f, 100.0f);
    my_shader.use();
    my_shader.setMat4("projection", projection);
    my_shader.setVec3("viewPos", camera->Position);

    my_shader.setVec3("dirLight.direction", glm::vec3(1.01f, 1.01f, 1.01f));
    my_shader.setVec3("dirLight.ambient", glm::vec3(0.01f, 0.01f, 0.02f));
    my_shader.setVec3("dirLight.diffuse", glm::vec3(1.0f, 1.0f, 1.0f));
    my_shader.setVec3("dirLight.specular", glm::vec3(1.0f, 1.0f, 1.0f));
    // point light 1
    my_shader.setVec3("pointLights[0].position", pointLightPositions[0]);
    my_shader.setVec3("pointLights[0].ambient", 0.05f, 0.05f, 0.05f);
    my_shader.setVec3("pointLights[0].diffuse", 0.8f, 0.8f, 0.8f);
    my_shader.setVec3("pointLights[0].specular", 1.0f, 1.0f, 1.0f);
    my_shader.setFloat("pointLights[0].constant", 1.0f);
    my_shader.setFloat("pointLights[0].linear", 0.09);
    my_shader.setFloat("pointLights[0].quadratic", 0.032);
    // point light 2
    my_shader.setVec3("pointLights[1].position", pointLightPositions[1]);
    my_shader.setVec3("pointLights[1].ambient", 0.05f, 0.05f, 0.05f);
    my_shader.setVec3("pointLights[1].diffuse", 0.8f, 0.8f, 0.8f);
    my_shader.setVec3("pointLights[1].specular", 1.0f, 1.0f, 1.0f);
    my_shader.setFloat("pointLights[1].constant", 1.0f);
    my_shader.setFloat("pointLights[1].linear", 0.09);
    my_shader.setFloat("pointLights[1].quadratic", 0.032);
    // point light 3
    my_shader.setVec3("pointLights[2].position", pointLightPositions[2]);
    my_shader.setVec3("pointLights[2].ambient", 0.05f, 0.05f, 0.05f);
    my_shader.setVec3("pointLights[2].diffuse", 0.6f, 0.1f, 0.8f);
    my_shader.setVec3("pointLights[2].specular", 1.0f, 1.0f, 1.0f);
    my_shader.setFloat("pointLights[2].constant", 1.0f);
    my_shader.setFloat("pointLights[2].linear", 0.09);
    my_shader.setFloat("pointLights[2].quadratic", 0.032);
    // point light 4
    my_shader.setVec3("pointLights[3].position", pointLightPositions[3]);
    my_shader.setVec3("pointLights[3].ambient", 0.05f, 0.05f, 0.05f);
    my_shader.setVec3("pointLights[3].diffuse", 0.1f, 1.1f, 0.8f);
    my_shader.setVec3("pointLights[3].specular", 1.0f, 1.0f, 1.0f);
    my_shader.setFloat("pointLights[3].constant", 1.0f);
    my_shader.setFloat("pointLights[3].linear", 0.09);
    my_shader.setFloat("pointLights[3].quadratic", 0.032);

    my_shader.setInt("texture1", 0);
    my_shader.setInt("texture2", 1);

    light_shader.use();
    glm::mat4 view = glm::mat4(1.0f);
    light_shader.setMat4("view", view);
    light_shader.setMat4("projection", projection);
    light_shader.setVec3("lightColor", glm::vec3(1.0, 1.0, 1.0));
}

void changePMV(Shader shader, Shader light_shader) {
    glm::mat4 view;
    view = camera->GetViewMatrix();
    shader.use();
    // 更新shader中观察坐标系
    shader.setMat4("view", view);
    // shader更新视角位置
    shader.setVec3("viewPos", camera->Position);
    light_shader.use();
    light_shader.setMat4("view", view);
}
/*-----------------------------------------------------------------------*/

// calculate tangent space
glm::vec3 calculate_tangent(glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, glm::vec2 uv1, glm::vec2 uv2, glm::vec2 uv3) {
    glm::vec3 edge1 = p2 - p1;
    glm::vec3 edge2 = p3 - p1;
    glm::vec2 deltaUV1 = uv2 - uv1;
    glm::vec2 deltaUV2 = uv3 - uv1;
    glm::vec3 proTangent;

    GLfloat proF = 1.0f / (deltaUV1.x * deltaUV2.y - deltaUV2.x * deltaUV1.y);

    proTangent.x = proF * (deltaUV2.y * edge1.x - deltaUV1.y * edge2.x);
    proTangent.y = proF * (deltaUV2.y * edge1.y - deltaUV1.y * edge2.y);
    proTangent.z = proF * (deltaUV2.y * edge1.z - deltaUV1.y * edge2.z);
    proTangent = glm::normalize(proTangent);
    return proTangent;
}



void get_vec3(std::vector<float> list, std::vector<glm::vec3> &vec)
{
    int n = list.size() / 3;
    for (int i = 0; i < n; i++)
    {
        vec.emplace_back(list[i * 3], list[i * 3 + 1], list[i * 3 + 2]);
    }
}
void get_vec2(std::vector<float> list, std::vector<glm::vec2>& vec)
{
    int n = list.size() / 2;
    for (int i = 0; i < n; i++)
    {
        vec.emplace_back(list[i * 2], list[i * 2 + 1]);
    }
}

void drawLights(Shader light_shader, glm::vec3 pointLightPositions[], unsigned int VAO) {
    light_shader.use();
    // 点光源
    for (int i = 0;i < 4;i++) {
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, pointLightPositions[i]);
        model = glm::scale(model, glm::vec3(0.2f));
        light_shader.setMat4("model", model);
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 36);
    }
    // 模拟线光源
    glm::mat4 model = glm::mat4(1.0f);
    model = glm::translate(model,  glm::vec3(-20.01f, -20.01f, -20.01f));
    model = glm::scale(model, glm::vec3(3.0f));
    light_shader.setMat4("model", model);
    glBindVertexArray(VAO);
    glDrawArrays(GL_TRIANGLES, 0, 36);
}

void drawCubes(Shader shader, const std::vector<RigidBody>& cubes, unsigned int VAO) {
    shader.use();

    for (auto cube : cubes) {
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, cube.get_transformation());
        model = glm::rotate(model, glm::radians(cube.get_rotation_angle()), cube.get_rotation_dir());
        shader.setMat4("model", model);
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 36);
    }
}



// 生成rigid body的函数
// 输入为初始位置
RigidBody create_body(glm::vec3 init_pos, glm::vec3 dir = glm::vec3(0,1,0), float angle = 0) {
    RigidBody new_body(init_pos, MASS, dir, angle);		//初始时init_pos为x(t)，质量设置为10
    new_body.setForce(glm::vec3(0.0f, 0.0f, -GRAVITY * MASS));		//初始时刻受到重力
    new_body.setIbody(MASS, 1.0f);
    return new_body;
}


// 返回值修改a，b中的成员变量
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
        v_rel = glm::dot(normal, (pat - pbt));
        glm::vec3 Ja, Jb;
        glm::vec3 det_va = v_rel * normal - con.a->get_vt();
        glm::vec3 det_vb = -v_rel * normal - con.b->get_vt();
        Ja = det_va * con.a->get_mass();
        Jb = det_vb * con.a->get_mass();
        glm::vec3 tao_a_impulse = glm::cross((con.particle_position - con.a->get_transformation()), Ja);
        glm::vec3 tao_b_impulse = glm::cross((con.particle_position - con.b->get_transformation()), Jb);
        // 更新body a和b的动量和角动量
        con.a->sum_Pt(Ja);
        con.a->sum_Lt(tao_a_impulse);
        con.b->sum_Pt(Jb);
        con.b->sum_Lt(tao_b_impulse);
    }
}
glm::vec4 face_function(glm::vec3 p, glm::vec3 n) {
    return glm::vec4(n, glm::dot(p, n));
}

glm::vec3 calculate_edge() {
    return glm::vec3(0.0f, 0.0f, 0.0f);
}

glm::vec3 calculate_face_vertex(std::vector<glm::vec3> line1, std::vector<glm::vec3> line2, std::vector<glm::vec3> line3) {
    glm::vec3 a, b, c, result;
    a = line1[0] + (line1[0] - line2[0]) / (line2[1] - line1[1]) * line1[1];
    b = line2[0] + (line2[0] - line3[0]) / (line3[1] - line2[1]) * line2[1];
    c = line1[0] + (line1[0] - line3[0]) / (line3[1] - line1[1]) * line1[1];
    result = (a + b + c);
    result /= 3;
    return result;
}

std::vector<glm::vec3> calculate_line(glm::vec4 func_1f, glm::vec4 func_2f, glm::vec3 n1, glm::vec3 n2) {
    std::vector<glm::vec3> result_line;
    result_line.emplace_back(0.0f);
    result_line.emplace_back(0.0f);
    float x, y, z;
    result_line[1] = glm::cross(n1, n2);
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
    result_line[0] = glm::vec3(x, y, z);
    return result_line;
}

bool judge_line_possibility(std::vector<glm::vec3> face_a, std::vector<glm::vec3> face_b, std::vector<glm::vec3> line) {
    bool line_face_a = false;
    glm::vec3 k;
    for (int i = 0; i < 4; i++) {
        if (i == 3) {
            k = face_a[0] - face_a[3];
        }
        else {
            k = face_a[i + 1] - face_a[i];
        }
        if ((face_a[i].x - line[0].x) / (line[1].x - k.x) >= 0 && (face_a[i].x - line[0].x) / (line[1].x - k.x) <= 1) {
            for (int j = 0; j < 4; j++) {
                if (j == 3) {
                    k = face_b[0] - face_b[3];
                }
                else {
                    k = face_b[j + 1] - face_b[j];
                }
                if ((face_b[j].x - line[0].x) / (line[1].x - k.x) >= 0 && (face_b[j].x - line[0].x) / (line[j].x - k.x) <= 1) {
                    line_face_a = true;
                }
            }
        }
    }
    return line_face_a;
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
        temp.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
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
    std::vector<glm::vec3> line;
    std::vector<std::vector<glm::vec3>> line_possible;
    std::vector<std::vector<int>> face_line_number;
    std::vector<int> temp_int;
    temp_int.clear();
    face_line_number.reserve(6);
    for (int i = 0; i < 6; i++) {
        face_line_number.push_back(temp_int);
    }
    for (int i = 0; i < 6; i++) {
        for (int j = i; j < 6; j++) {
            // 传入face_function面上的一个点和一个法向量
            // 计算出面的表达式ax+by+cz=d 结果存在一个vec4里面
            func_1f = face_function(face_p_n_a[i][0], face_p_n_a[i][4]);
            func_2f = face_function(face_p_n_b[j][0], face_p_n_b[j][4]);
            // result[0] 一个点 result[1] 方向   直线表达形式(a,b,c) + t * (d, e, f)
            line = calculate_line(func_1f, func_2f, face_p_n_a[i][4], face_p_n_b[j][4]);
            if (judge_line_possibility(face_p_n_a[i], face_p_n_b[j], line)) {
                line_possible.push_back(line);
                face_line_number[i].push_back(line_possible.size());
                face_line_number[j].push_back(line_possible.size());
            }
        }
    }
    int num_of_face = 0;
    std::vector<int> target_face;
    for (int i = 0; i < 6; i++) {
        if (face_line_number[i].size() > 1) {
            target_face.push_back(i);
            num_of_face += 1;
        }
    }
    if (num_of_face == 1) {
        f_v.is_valid = 1;
        f_v.is_face_vertex = 1;
        f_v.face_normal = face_normal[target_face[0]];
        f_v.particle_position = calculate_face_vertex(line_possible[face_line_number[target_face[0]][0]], line_possible[face_line_number[target_face[0]][1]], line_possible[face_line_number[target_face[0]][2]]);
    }
    else if (num_of_face == 2) {
        f_v.is_valid = 1;
        f_v.is_face_vertex = 0;
        f_v.particle_position = calculate_edge();
    }
    else {
        f_v.is_valid = 0;
    }
    return f_v;
}



void process_gravity_floor(RigidBody &body) {
    if (body.get_transformation().y < 0.8) {
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
        repeat_index.push_back(-1);
        for (int i = 0;i < points.size();i++) {
            if (points[i].y < min_y) {
                repeat_index.clear();
                repeat_index.push_back(i);
                min_y = points[i].y;
            } else if (points[i].y - min_y < 0.01) {
                repeat_index.push_back(i);
            }
        }
        if (min_y <= 0.01) {
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

            if (body.get_Pt().y > 2) {
                // 初次碰撞瞬间失去一部分动量
                body.sum_Pt(glm::vec3(0,-(body.get_Pt().y - 2.5),0));
            }


            glm::vec3 delt_v = glm::vec3(0, 1.0f, 0.0f);
            // 反弹为0.5的重力
            delt_v *= GRAVITY * time_interval * 0.5f;
            glm::vec3 J = delt_v * body.get_mass();
            body.sum_Pt(J);

            // 更新角动量
            glm::vec3 tao_impulse = glm::cross((hit_point - body.get_transformation()), J);
            tao_impulse.x = tao_impulse.x < 0.0001 ? 0 : tao_impulse.x;
            tao_impulse.y = tao_impulse.y < 0.0001 ? 0 : tao_impulse.y;
            tao_impulse.z = tao_impulse.z < 0.0001 ? 0 : tao_impulse.z;

            body.sum_Lt(tao_impulse * slow * 10.0f);
            return ;
        }


    }
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

void move_bodies(RigidBody &body) {
    // 根据动量角动量移动物体
    glm::vec3 curr_v = body.get_Pt() / (float) MASS;				//计算并更新此时刻物体的 linear velocity
    glm::vec3 curr_w = body.get_Lt() / (float) body.get_Ibody();	//计算并更新此时刻物体的 angular velocity
    body.UpdateStates(curr_v, curr_w);
}


void update_cube_positions(std::vector<RigidBody> &cubes) {
    // todo: 计算所有的物体碰撞后改变的速度和角速度
    // todo: 按照输入的时间间隔移动cubes位置（修改RigidBody->transformation以及旋转

    //std::cout<<"new round"<<std::endl;
    // step1:
    for (int i = 0; i < cubes.size(); i++) {
        for (int j = i + 1; j < cubes.size(); j++) {
            float dis = glm::length(cubes[i].get_transformation() - cubes[j].get_transformation());
            if (dis < sqrtf(3)) {
                cubes[i].possible_collision.push_back(&cubes[j]);
                cubes[j].possible_collision.push_back(&cubes[i]);
            }
        }
    }
    std::vector<Contact> contacts;
    // step2: 计算碰撞点
    for (auto &cube : cubes) {
        while (!cube.possible_collision.empty()) {
            RigidBody *tem = cube.possible_collision[cube.possible_collision.size() - 1];
            Contact tem_contact = check_collision(*tem, cube);
            if (tem_contact.is_valid) {
                contacts.push_back(tem_contact);
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
    for (auto &cube: cubes) {
        // 重力和地板
        //std::cout<<"before inside: "<<std::endl;
        //print_vec3(cubes[0].get_Pt());
        process_gravity_floor(cube);
        //std::cout<<"after inside1: "<<std::endl;
        //print_vec3(cubes[0].get_Pt());
        move_bodies(cube);
        //std::cout<<"after inside2: "<<std::endl;
        //print_vec3(cubes[0].get_Pt());
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
    std::vector<glm::vec3> result = calculate_line(func_1f, func_2f, a_n, b_n);
    glm::vec3 test_p(1.5,0.25,0.5);
    if (glm::normalize(test_p - result[0]) != glm::normalize(result[1])
    || glm::normalize(b - result[0]) != glm::normalize(result[1])) {
        std::cout<<"calculate_line error!"<<std::endl;
    }
}



int main()
{
    check_calculate_line();


    std::string root_dir = "/Users/TT/Desktop/CS171/RIgif-Body-Simulation";
    int len = root_dir.length();
    std::string model_dir = root_dir + "/model";

    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this statement to fix compilation on OS X
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    glEnable(GL_DEPTH_TEST);
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Here you need to fill construct function of class Shader. And you need to understand other funtions in Shader.//
    // Then, write code in shader_m.vs, shader_m.fs and shader_m.gs to finish the tasks.                             //
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Shader my_shader(
            (root_dir + "/src/shader_m.vs").c_str(),
            (root_dir + "/src/shader_m.fs").c_str()
    );
    //A shader for light visiable source
    Shader lampShader((root_dir + "/src/lamp.vs").c_str(), (root_dir + "/src/lamp.fs").c_str());


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Render a box to show nice normal mapping.                                                                   //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    float raw_vertices_cube_0[] = {

            // positions          // normals           // texture coords

            -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  0.0f,  0.0f,
            0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  1.0f,  0.0f,
            0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  1.0f,  1.0f,
            0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  1.0f,  1.0f,
            -0.5f,  0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  0.0f,  1.0f,
            -0.5f, -0.5f, -0.5f,  0.0f,  0.0f, -1.0f,  0.0f,  0.0f,

            -0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f,
            0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  1.0f,  0.0f,
            0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  1.0f,  1.0f,
            0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  1.0f,  1.0f,
            -0.5f,  0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  0.0f,  1.0f,
            -0.5f, -0.5f,  0.5f,  0.0f,  0.0f,  1.0f,  0.0f,  0.0f,

            -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,  1.0f,  0.0f,
            -0.5f,  0.5f, -0.5f, -1.0f,  0.0f,  0.0f,  1.0f,  1.0f,
            -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
            -0.5f, -0.5f, -0.5f, -1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
            -0.5f, -0.5f,  0.5f, -1.0f,  0.0f,  0.0f,  0.0f,  0.0f,
            -0.5f,  0.5f,  0.5f, -1.0f,  0.0f,  0.0f,  1.0f,  0.0f,

            0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,  1.0f,  0.0f,
            0.5f,  0.5f, -0.5f,  1.0f,  0.0f,  0.0f,  1.0f,  1.0f,
            0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
            0.5f, -0.5f, -0.5f,  1.0f,  0.0f,  0.0f,  0.0f,  1.0f,
            0.5f, -0.5f,  0.5f,  1.0f,  0.0f,  0.0f,  0.0f,  0.0f,
            0.5f,  0.5f,  0.5f,  1.0f,  0.0f,  0.0f,  1.0f,  0.0f,

            -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,  0.0f,  1.0f,
            0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,  1.0f,  1.0f,
            0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,  1.0f,  0.0f,
            0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,  1.0f,  0.0f,
            -0.5f, -0.5f,  0.5f,  0.0f, -1.0f,  0.0f,  0.0f,  0.0f,
            -0.5f, -0.5f, -0.5f,  0.0f, -1.0f,  0.0f,  0.0f,  1.0f,

            -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,  0.0f,  1.0f,
            0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,  1.0f,  1.0f,
            0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,  1.0f,  0.0f,
            0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,  1.0f,  0.0f,
            -0.5f,  0.5f,  0.5f,  0.0f,  1.0f,  0.0f,  0.0f,  0.0f,
            -0.5f,  0.5f, -0.5f,  0.0f,  1.0f,  0.0f,  0.0f,  1.0f

    };

    float vertices_cube_0[36 * 11];
    for (int l = 0; l < 36; l += 3) {
        std::vector<glm::vec3> points;
        std::vector<glm::vec2> uvs;
        std::vector<float> cod;
        std::vector<float> uv;
        for (int j = 0; j < 3; ++j) {
            for (int i = 0; i < 3; ++i) {
                cod.push_back(raw_vertices_cube_0[(l + j) * 8 + i]);
                // std::cout<<"push point: "<<raw_vertices_cube_0[(l + j) * 8 + i];
            }
            // std::cout<<std::endl;
            for (int i = 0; i < 2; ++i) {
                uv.push_back(raw_vertices_cube_0[(l + j) * 8 + i + 6]);
                // std::cout<<"push uvs: "<<raw_vertices_cube_0[(l + j) * 8 + i + 6];
            }
            // std::cout<<std::endl;
        }
        get_vec3(cod, points);
        get_vec2(uv, uvs);
        glm::vec3 proTangent = calculate_tangent(points[0], points[1], points[2], uvs[0], uvs[1], uvs[2]);
        for (int i = 0; i < 8; ++i) {
            vertices_cube_0[l * 11 + i] = raw_vertices_cube_0[l * 8 + i];
            vertices_cube_0[(l + 1) * 11 + i] = raw_vertices_cube_0[(l + 1) * 8 + i];
            vertices_cube_0[(l + 2) * 11 + i] = raw_vertices_cube_0[(l + 2) * 8 + i];
        }
        for (int k = 0; k < 3; ++k) {
            vertices_cube_0[(l + k) * 11 + 8] = proTangent.x;
            vertices_cube_0[(l + k) * 11 + 9] = proTangent.y;
            vertices_cube_0[(l + k) * 11 + 10] = proTangent.z;
        }
    }

    unsigned int VBO, VAO;
    // the first argument is id
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).
    glBindVertexArray(VAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices_cube_0), vertices_cube_0, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)nullptr);
    glEnableVertexAttribArray(0);
    // normals
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    // texture coords
    glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
    // tangent
    glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)(8 * sizeof(float)));
    glEnableVertexAttribArray(3);

    unsigned int lightVAO;
    glGenVertexArrays(1, &lightVAO);

    glBindVertexArray(lightVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)nullptr);
    glEnableVertexAttribArray(0);

    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 11 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(1);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // You need to fill this function which is defined in my_texture.h. The parameter is the path of your image.   //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    unsigned int texture1 = loadTexture((model_dir + "/brickwall.jpg").c_str());
    unsigned int texture2 = loadTexture((model_dir + "/normal_map.jpg").c_str());

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Here we defined pointlights in shader and passed some parameter for you. You can take this as an example.   //
    // Or you can change it if you like.                                                                           //
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    glm::vec3 pointLightPositions[] = {
            glm::vec3(5.7f,  5.2f,  2.0f),
            glm::vec3(0.0f, 2.3f, 0.0f),
            glm::vec3(-4.0f,  2.0f, -12.0f),
            glm::vec3(0.0f,  0.0f, -3.0f)
    };

    // todo: 这行以上的代码都不用看了， 是设置cube的模型数据的, 下面的是最重要的!!!!


    // 用来存cube的位置，每次循环位置都会改变。
    std::vector<RigidBody> CubePositions;

    // 以下为使用方法
    CubePositions.push_back(create_body(glm::vec3(0.0f, 3.5f, 0.0f)));
    CubePositions.push_back(create_body(glm::vec3(0.0f, 6.0f, 0.0f),glm::vec3(1,1,0), 45));

    initPMV(my_shader, lampShader, pointLightPositions);

    while (!glfwWindowShouldClose(window))
    {
        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        changePMV(my_shader, lampShader);

        //std::cout<<"before: "<<std::endl;
        //print_vec3(CubePositions[0].get_transformation());
        update_cube_positions(CubePositions);
        //std::cout<<"after: "<<std::endl;
        //print_vec3(CubePositions[0].get_transformation());


        //Update Camera Matrix
        glFlush();
        glEnable(GL_MULTISAMPLE);
        glEnable(GL_LIGHTING);
        glEnable(GL_COLOR_MATERIAL);
        glLightModeli(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //  Render an object using texture and normal map.                                                             //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 现在有的带normal map 的图形在这里画的
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture1);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, texture2);

        // 画地板
        my_shader.use();
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::scale(model, glm::vec3(20.0f));
        model = glm::translate(model, glm::vec3(0.0f, -0.5f, 0.0f));
        my_shader.setMat4("model", model);
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 36);


        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //  Render the object in .obj file. You need to set materials and wrap texture for objects.                    //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        /*
        my_shader.use();
        model = glm::mat4(1.0f);
        my_shader.setMat4("model", model);
        glBindVertexArray(obj_VAO_l[0]);
        glDrawArrays(GL_TRIANGLES, 0, shapes[0].mesh.indices.size());

        my_shader.use();
        model = glm::mat4(1.0f);
        my_shader.setMat4("model", model);
        glBindVertexArray(obj_VAO_l[1]);
        glDrawArrays(GL_TRIANGLES, 0, shapes[1].mesh.indices.size());
         */



        /////////////////////////////////////////////////////////////////////

        /////////////////////////////end/////////////////////////////////////
        // Visualize point lights and dir lights
        glActiveTexture(GL_TEXTURE0);

        // 画出的光源都是没带normal map的砖头墙
        drawLights(lampShader, pointLightPositions, VAO);
        drawCubes(my_shader, CubePositions, VAO);

        // todo: 从天上生成新的cube落下
        // CubePositions.push_back(create_body(glm::vec3(2.0f, 3.5f, 4.5f)));

        glfwSwapBuffers(window);
        glfwPollEvents();
        glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();

    glDeleteVertexArrays(1, &VAO);
    glDeleteVertexArrays(1, &lightVAO);
    glDeleteBuffers(1, &VBO);

    delete camera;
    return 0;
}

