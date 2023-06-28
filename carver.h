// Computer Graphics Project Assignment 6 Digital Photography Seam Carving
// carver.h
#pragma once
#ifndef CARVER_H

#include <Eigen/Dense>
#include <string>
#include <vector>

enum Direction {
    VERTICAL,
    HORTIZONTAL
};

class Carver {
public:
    Carver();
    Carver(std::string fileName);
    ~Carver();

    void load(std::string fileName);
    void save(std::string fileName);

    void resize(int width, int height);

    Eigen::Vector3d getPixel(int x, int y) const;
    void setPixel(int x, int y, Eigen::Vector3d color);
    
    double getEnergy(int x, int y) const;
    void setEnergy(int x, int y, double strength);

private:
    void calculateEnergy();
    double* calculateFlow(const Direction direction);
    std::vector<int> calculateSeam(const Direction direction);
    void removeSeam(const Direction direction);
    void addSeam(const Direction direction);
    double * convolve(Eigen::Matrix3d kernel);

    double *energy;
    Eigen::Vector3d *image;

    int currentWidth;
    int currentHeight;
    int originalWidth;
    int originalHeight;
    int depth;
    int spectrum;
};

#endif
