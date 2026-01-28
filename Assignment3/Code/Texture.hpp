//
// Created by LEI XU on 4/27/19.
//

#ifndef RASTERIZER_TEXTURE_H
#define RASTERIZER_TEXTURE_H
#include "global.hpp"
#include <algorithm>
#include <cmath>
#include <eigen3/Eigen/Eigen>
#include <opencv2/opencv.hpp>
class Texture{
private:
    cv::Mat image_data;

public:
    Texture(const std::string& name)
    {
        image_data = cv::imread(name);
        cv::cvtColor(image_data, image_data, cv::COLOR_RGB2BGR);
        width = image_data.cols;
        height = image_data.rows;
    }

    int width, height;

    Eigen::Vector3f getColor(float u, float v)
    {
        auto u_img = u * width;
        auto v_img = (1 - v) * height;
        auto color = image_data.at<cv::Vec3b>(v_img, u_img);
        return Eigen::Vector3f(color[0], color[1], color[2]);
    }

    Eigen::Vector3f getColorBilinear(float u, float v)
    {
        u = std::clamp(u, 0.0f, 1.0f);
        v = std::clamp(v, 0.0f, 1.0f);

        float u_img = u * width;
        float v_img = (1 - v) * height;

        int x0 = static_cast<int>(std::floor(u_img));
        int y0 = static_cast<int>(std::floor(v_img));
        int x1 = std::min(x0 + 1, width - 1);
        int y1 = std::min(y0 + 1, height - 1);

        float sx = u_img - x0;
        float sy = v_img - y0;

        const cv::Vec3b c00 = image_data.at<cv::Vec3b>(y0, x0);
        const cv::Vec3b c10 = image_data.at<cv::Vec3b>(y0, x1);
        const cv::Vec3b c01 = image_data.at<cv::Vec3b>(y1, x0);
        const cv::Vec3b c11 = image_data.at<cv::Vec3b>(y1, x1);

        const float w00 = (1.0f - sx) * (1.0f - sy);
        const float w10 = sx * (1.0f - sy);
        const float w01 = (1.0f - sx) * sy;
        const float w11 = sx * sy;

        return Eigen::Vector3f(
            w00 * c00[0] + w10 * c10[0] + w01 * c01[0] + w11 * c11[0],
            w00 * c00[1] + w10 * c10[1] + w01 * c01[1] + w11 * c11[1],
            w00 * c00[2] + w10 * c10[2] + w01 * c01[2] + w11 * c11[2]
        );
    }

};
#endif //RASTERIZER_TEXTURE_H
