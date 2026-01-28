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
    inline float lerp(float x, float v0, float v1)
    {
        return v0 + x * (v1 - v0);
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

        // s, t exactly as in the slide
        float s = u_img - x0;
        float t = v_img - y0;

        const cv::Vec3b u00 = image_data.at<cv::Vec3b>(y0, x0);
        const cv::Vec3b u10 = image_data.at<cv::Vec3b>(y0, x1);
        const cv::Vec3b u01 = image_data.at<cv::Vec3b>(y1, x0);
        const cv::Vec3b u11 = image_data.at<cv::Vec3b>(y1, x1);

        Eigen::Vector3f result;
        for (int c = 0; c < 3; c++)
        {
            float u0 = lerp(s, u00[c], u10[c]);
            float u1 = lerp(s, u01[c], u11[c]);
            result[c] = lerp(t, u0, u1);
        }

        return result;
    }


};
#endif //RASTERIZER_TEXTURE_H
