#include "Triangle.hpp"
#include "rasterizer.hpp"
#include <eigen3/Eigen/Eigen>
#include <iostream>
#include <opencv2/opencv.hpp>

constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
    float angle = MY_PI * rotation_angle / 180.0;
    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.
    model << cos(angle), -sin(angle), 0, 0,
             sin(angle),  cos(angle), 0, 0,
                         0,             0, 1, 0,
                         0,             0, 0, 1;
    return model;
}

Eigen::Matrix4f get_rotation(Vector3f axis, float angle) {
    Eigen::Matrix4f rotation = Eigen::Matrix4f::Identity();

    // 1. 角度：degree → radian
    float rad = angle * M_PI / 180.0f;

    // 2. 轴向量单位化（非常重要）
    axis.normalize();

    float x = axis.x();
    float y = axis.y();
    float z = axis.z();

    // 3. Rodrigues 公式
    Eigen::Matrix3f R;
    R << cos(rad) + x*x*(1 - cos(rad)),      x*y*(1 - cos(rad)) - z*sin(rad),  x*z*(1 - cos(rad)) + y*sin(rad),
         y*x*(1 - cos(rad)) + z*sin(rad),    cos(rad) + y*y*(1 - cos(rad)),    y*z*(1 - cos(rad)) - x*sin(rad),
         z*x*(1 - cos(rad)) - y*sin(rad),    z*y*(1 - cos(rad)) + x*sin(rad),  cos(rad) + z*z*(1 - cos(rad));

    // 4. 放进 4×4 齐次矩阵
    rotation.block<3,3>(0,0) = R;

    return rotation;
}

Eigen::Matrix4f get_model_matrix(Vector3f axis, float rotation_angle)
{
    return get_rotation(axis, rotation_angle);
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    // Students will implement this function

    Eigen::Matrix4f projection = Eigen::Matrix4f::Zero();
    float fov_rad = eye_fov * M_PI / 180.0f;

    float tan_half_fov = tan(fov_rad / 2.0f);
    projection(0,0) = 1.0 / (tan_half_fov * aspect_ratio);
    projection(1,1) = 1.0 / tan_half_fov;
    projection(2,2) = (zNear + zFar) / (zNear - zFar);
    projection(2,3) = -2 * zNear * zFar / (zNear - zFar);
    projection(3,2) = 1.0; 
    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.

    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";
    Eigen::Vector3f axis(0,0,1); // axis-z by default
    
    
if (argc >= 3) {
    command_line = true;
    angle = std::stof(argv[2]);

    int arg_idx = 3;

    // 如果还有参数，看看是不是 axis
    if (argc > arg_idx) {
        char c = argv[arg_idx][0];
        if (c == 'x' || c == 'X') {
            axis = Eigen::Vector3f(1, 0, 0);
            arg_idx++;
        }
        else if (c == 'y' || c == 'Y') {
            axis = Eigen::Vector3f(0, 1, 0);
            arg_idx++;
        }
        else if (c == 'z' || c == 'Z') {
            axis = Eigen::Vector3f(0, 0, 1);
            arg_idx++;
        }
        // 否则：不是轴，就当它是文件名
    }

    // 剩下的参数当作文件名
    if (argc > arg_idx) {
        filename = std::string(argv[arg_idx]);
    }
}


    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0, 0, 5};

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};

    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    if (command_line) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(axis, angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));
        
        r.draw(pos_id, ind_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);

        cv::imwrite(filename, image);

        return 0;
    }

    while (key != 27) {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(axis,angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';

        if (key == 'a') {
            angle += 10;
        }
        else if (key == 'd') {
            angle -= 10;
        }
    }

    return 0;
}
