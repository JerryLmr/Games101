#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include<iostream>

int main(){

    // Basic Example of cpp
    std::cout << "Example of cpp \n";
    float a = 1.0, b = 2.0;
    std::cout << a << std::endl;
    std::cout << a/b << std::endl;
    std::cout << std::sqrt(b) << std::endl;
    std::cout << std::acos(-1) << std::endl;
    std::cout << std::sin(30.0/180.0*acos(-1)) << std::endl;

    // Example of vector
    std::cout << "Example of vector \n";
    // vector definition
    Eigen::Vector3f v(1.0f,2.0f,3.0f);
    Eigen::Vector3f w(1.0f,0.0f,0.0f);
    // vector output
    std::cout << "Example of output \n";
    std::cout << v << std::endl;
    // vector add
    std::cout << "Example of add \n";
    std::cout << v + w << std::endl;
    // vector scalar multiply
    std::cout << "Example of scalar multiply \n";
    std::cout << v * 3.0f << std::endl;
    std::cout << 2.0f * v << std::endl;

    // Example of matrix
    std::cout << "Example of matrix \n";
    // matrix definition
    Eigen::Matrix3f i,j;
    i << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    j << 2.0, 3.0, 1.0, 4.0, 6.0, 5.0, 9.0, 7.0, 8.0;
    // matrix output
    std::cout << "Example of output \n";
    std::cout << i << std::endl;
    // matrix add i + j
    // matrix scalar multiply i * 2.0
    // matrix multiply i * j
    // matrix multiply vector i * v


    /* 
    * PA 0
    */

    // TO DO: Define point P
    Eigen::Vector3d P(2,1,1);

    // TO DO: Define rotation matrix M

    // 旋转角度45度转弧度
    double theta = 45.0 * M_PI / 180;
    Eigen::Matrix3d R;
    R << cos(theta), -sin(theta), 0,
         sin(theta),  cos(theta), 0,
         0,            0,              1;
    Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
    T(0, 2) = 1;  // x方向平移1
    T(1, 2) = 2;  // y方向平移2

    // TO DO: M * P
    Eigen::Matrix3d M = T * R;
    Eigen::Vector3d P_transformed = M * P;

    std::cout << "变换后的点坐标 (齐次坐标):\n" << P_transformed << std::endl;
    std::cout << "变换后的点坐标 (普通二维坐标): (" 
              << P_transformed(0) / P_transformed(2) << ", "
              << P_transformed(1) / P_transformed(2) << ")\n";
    return 0;
}