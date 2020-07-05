
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

Camera *camera;

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
    camera = new Camera (glm::vec3(0.0f, 0.0f, 3.0f));

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
        // todo: add rotation
        shader.setMat4("model", model);
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 36);
    }
}



// 生成rigid body的函数
// 输入为初始位置
RigidBody create_body(glm::vec3 init_pos) {
    // todo: 创建一个rigidbody并且生成其中的点云等初始化数据
    return RigidBody(glm::vec3(2,2,2), glm::mat4(1.0f), 150, std::vector<Particle> {});
}

void update_cube_positions(std::vector<RigidBody> cubes, double time_interval) {
    // todo: 计算所有的物体碰撞后改变的速度和角速度
    // todo: 按照输入的时间间隔移动cubes位置（修改RigidBody->transformation以及旋转
}



int main()
{
    std::string root_dir = "C:/Users/38182/Desktop/cg learning OpenGL/project/Rigid-Body-Simulation";
    int len = root_dir.length();
    // cmakeLists.txt所在文件目录绝对位置
    std::string model_dir = root_dir + "/model";
    // std::cout<<root_dir<<std::endl;
    // todo: 如果这里print出来的model_dir或root_dir值不对，无法修改，可以改成绝对路径写死。

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
            glm::vec3(2.3f, -3.3f, -4.0f),
            glm::vec3(-4.0f,  2.0f, -12.0f),
            glm::vec3(0.0f,  0.0f, -3.0f)
    };

    // todo: 这行以上的代码都不用看了， 是设置cube的模型数据的, 下面的是最重要的!!!!


    // 用来存cube的位置，每次循环位置都会改变。
    std::vector<RigidBody> CubePositions;

    // 以下为使用方法
    CubePositions.push_back(create_body(glm::vec3(2.0f, 3.5f, 4.5f)));

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
        update_cube_positions(CubePositions, ((double) (1))/30);
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

        my_shader.use();
        glm::mat4 model = glm::mat4(1.0f);
        model = glm::translate(model, glm::vec3(-1.5f, -1.5f, -1.5f));
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
        CubePositions.push_back(create_body(glm::vec3(2.0f, 3.5f, 4.5f)));

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

