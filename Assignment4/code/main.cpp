#include <chrono>
#include <iostream>
#include <opencv2/opencv.hpp>

std::vector<cv::Point2f> control_points;

void blend_pixel(cv::Mat &window, int x, int y, float intensity)
{
    if (x < 0 || x >= window.cols || y < 0 || y >= window.rows)
    {
        return;
    }

    auto &pixel = window.at<cv::Vec3b>(y, x);
    const float blended = std::min(255.0f, pixel[1] + 255.0f * intensity);
    pixel[1] = static_cast<uchar>(blended);
}

void mouse_handler(int event, int x, int y, int flags, void *userdata) 
{
    if (event == cv::EVENT_LBUTTONDOWN && control_points.size() < 4) 
    {
        std::cout << "Left button of the mouse is clicked - position (" << x << ", "
        << y << ")" << '\n';
        control_points.emplace_back(x, y);
    }     
}

void naive_bezier(const std::vector<cv::Point2f> &points, cv::Mat &window) 
{
    auto &p_0 = points[0];
    auto &p_1 = points[1];
    auto &p_2 = points[2];
    auto &p_3 = points[3];

    for (double t = 0.0; t <= 1.0; t += 0.001) 
    {
        auto point = std::pow(1 - t, 3) * p_0 + 3 * t * std::pow(1 - t, 2) * p_1 +
                 3 * std::pow(t, 2) * (1 - t) * p_2 + std::pow(t, 3) * p_3;

        window.at<cv::Vec3b>(point.y, point.x)[2] = 255;
    }
}

cv::Point2f recursive_bezier(const std::vector<cv::Point2f> &control_points, float t) 
{
    if (control_points.size() == 1)
    {
        return control_points[0];
    }

    std::vector<cv::Point2f> next_control_points;
    next_control_points.reserve(control_points.size() - 1);

    for (size_t i = 0; i + 1 < control_points.size(); ++i)
    {
        next_control_points.emplace_back(
            (1.0f - t) * control_points[i] + t * control_points[i + 1]);
    }

    return recursive_bezier(next_control_points, t);

}

void bezier(const std::vector<cv::Point2f> &control_points, cv::Mat &window) 
{
    for (float t = 0.0f; t <= 1.0f; t += 0.001f)
    {
        const auto point = recursive_bezier(control_points, t);
        const int x0 = static_cast<int>(std::floor(point.x));
        const int y0 = static_cast<int>(std::floor(point.y));

        for (int dy = 0; dy <= 1; ++dy)
        {
            for (int dx = 0; dx <= 1; ++dx)
            {
                const int px = x0 + dx;
                const int py = y0 + dy;
                const float center_x = px + 0.5f;
                const float center_y = py + 0.5f;
                const float weight_x = std::max(0.0f, 1.0f - std::abs(point.x - center_x));
                const float weight_y = std::max(0.0f, 1.0f - std::abs(point.y - center_y));

                blend_pixel(window, px, py, weight_x * weight_y);
            }
        }
    }

}

int main() 
{
    cv::Mat window = cv::Mat(700, 700, CV_8UC3, cv::Scalar(0));
    cv::cvtColor(window, window, cv::COLOR_BGR2RGB);
    cv::namedWindow("Bezier Curve", cv::WINDOW_AUTOSIZE);

    cv::setMouseCallback("Bezier Curve", mouse_handler, nullptr);

    int key = -1;
    while (key != 27) 
    {
        for (auto &point : control_points) 
        {
            cv::circle(window, point, 3, {255, 255, 255}, 3);
        }

        if (control_points.size() == 4) 
        {
            naive_bezier(control_points, window);
            bezier(control_points, window);

            cv::imshow("Bezier Curve", window);
            cv::imwrite("my_bezier_curve.png", window);
            key = cv::waitKey(0);

            return 0;
        }

        cv::imshow("Bezier Curve", window);
        key = cv::waitKey(20);
    }

return 0;
}
