// Computer Graphics Project Assignment 6 Digital Photography Seam Carving
// carver.cpp
#define cimg_display 0
#include "CImg.h"
#include "carver.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

int count = 0;

// Create Scharr Convolution operators for horizontal and vertical operation.
const Eigen::Matrix3d scharrHorizontal{
    {47, 162, 47},
    {0, 0, 0},
    {-47, -162, -47}};
const Eigen::Matrix3d scharrVertical{
    {47, 0, -47},
    {162, 0, -162},
    {47, 0, -47}};

// Create Sobel Convolution operators for horizontal and vertical operation.
const Eigen::Matrix3d sobelHorizontal{
    {3, 10, 3},
    {0, 0, 0},
    {-3, -10, -3}};
const Eigen::Matrix3d sobelVertical{
    {3, 0, -3},
    {10, 0, -10},
    {3, 0, -3}};

// Default Constructor.
Carver::Carver() : image(nullptr), energy(nullptr), currentWidth(0), currentHeight(0), originalWidth(0), originalHeight(0), depth(0), spectrum(0) {}

// Parameterized Constructor used to load a file.
Carver::Carver(std::string fileName) : image(nullptr), energy(nullptr), currentWidth(0), currentHeight(0), originalWidth(0), originalHeight(0), depth(0), spectrum(0)
{
    load(fileName);
}

// Create a Destructor for clearing allocated memory in the Image and Energy Buffers.
Carver::~Carver()
{
    if (image != nullptr)
    {
        delete[] image;
        image = nullptr;
    }
    if (energy != nullptr)
    {
        delete[] energy;
        energy = nullptr;
    }
}

// The Public Member Function "void Carver::load(std::string filename)" loads the image specified in String Variable that is "string fileName" and calculates the
// energy.
void Carver::load(std::string fileName)
{
    // Use the Cimg Library to load the image.
    cimg_library::CImg<double> input(fileName.c_str());

    // Record the size, depth, and spectrum information from the image.
    originalWidth = currentWidth = input.width();
    originalHeight = currentHeight = input.height();
    depth = input.depth();
    spectrum = input.spectrum();

    // Convert the image from RGB to LAB to make working with it simpler.
    cimg_library::CImg<double> lab = input.RGBtoLab();

    // Generate the buffers for image and energy.
    image = new Eigen::Vector3d[currentWidth * currentHeight];
    energy = new double[currentWidth * currentHeight];
    for (int i = 0; i < currentWidth; i++)
    {
        for (int j = 0; j < currentHeight; j++)
        {
            // Set the Pixel Values from the image.
            setPixel(i, j, Eigen::Vector3d(lab(i, j, 0), lab(i, j, 1), lab(i, j, 2)));

            // Set the energy at the Pixel to default.
            setEnergy(i, j, 255.0);
        }
    }

    // Calculate the Energy Buffer in order to do Seam Carving.
    calculateEnergy();
}


// The Public Member Function "void Carver::save(std::string fileName)" saves the image that is generated.
void Carver::save(std::string fileName)
{
    // Prepare a CImg object to save.
    cimg_library::CImg<double> output(currentWidth, currentHeight, depth, spectrum, 0);

    // Create a for() Loop to loop through and record the Pixels in the desired area of the image to the output.
    double maxEnergy = 0;
    for (int iIndex = 0; iIndex < currentWidth; iIndex++)
    {
        for (int jIndex = 0; jIndex < currentHeight; jIndex++)
        {
            maxEnergy = std::max(getEnergy(iIndex, jIndex), maxEnergy);
            Eigen::Vector3d pixel = getPixel(iIndex, jIndex);
            output(iIndex, jIndex, 0) = pixel[0];
            output(iIndex, jIndex, 1) = pixel[1];
            output(iIndex, jIndex, 2) = pixel[2];
        }
    }

    // Convert the output from LAB to RGB for saving.
    cimg_library::CImg<double> rgb = output.LabtoRGB();
    if (strstr(fileName.c_str(), "png"))
        rgb.save_png(fileName.c_str());
    else if (strstr(fileName.c_str(), "jpg"))
        rgb.save_jpeg(fileName.c_str());

    // In the same manner as converting the output from LAB to RGB for saving, do the same process to the Energy Map for viewing.
    cimg_library::CImg<double> energyOutput(currentWidth, currentHeight, depth, spectrum, 0);

    for (int iIndex = 0; iIndex < currentWidth; iIndex++)
    {
        for (int jIndex = 0; jIndex < currentHeight; jIndex++)
        {
            double grey = ((getEnergy(iIndex, jIndex) / maxEnergy), 1.0/3.0);
            energyOutput(iIndex, jIndex, 0) = grey;
            energyOutput(iIndex, jIndex, 1) = grey;
            energyOutput(iIndex, jIndex, 2) = grey;
        }
    }

    // Always save the Energy Map as a Portable Network Graphics PNG Image for lossless viewing.
    energyOutput.save_png("Energy.png");
}

// The Public Member Function "void Carver::resize(int width, int height)" performs the resize.
void Carver::resize(int width, int height)
{
    // Check if the image becomes larger and generate a New Buffer set for that size.
    if (width > originalWidth || height > originalHeight)
    {
        int nWidth = std::max(width, originalWidth);
        int nHeight = std::max(height, originalHeight);
        if (image != nullptr)
        {
            // If the image is generated, then copy the image into the New Buffer.
            Eigen::Vector3d *newImage = new Eigen::Vector3d[nWidth * nHeight];
            double *newEnergy = new double[nWidth * nHeight];

            for (int i = 0; i < nWidth; ++i)
            {
                for (int j = 0; j < nHeight; ++j)
                {
                    int oIndex = i * originalHeight + j;
                    int nIndex = i * nHeight + j;

                    if (i >= originalWidth || j >= originalHeight)
                    {
                        newImage[nIndex] = Eigen::Vector3d::Zero();
                        newEnergy[nIndex] = 255.0;
                        continue;
                    }
                    newImage[nIndex] = image[oIndex];
                    newEnergy[nIndex] = energy[oIndex];
                }
            }

            delete[] image;
            delete[] energy;

            image = newImage;
            energy = newEnergy;

            originalHeight = nHeight;
            originalWidth = nWidth;
        }
        else
        {
            // If an image is not generated, then generate a blank image in order to prevent errors.
            image = new Eigen::Vector3d[nWidth * nHeight];
            energy = new double[nWidth * nHeight];

            for (int iIndex = 0; iIndex < nWidth; ++iIndex)
            {
                for (int jIndex = 0; jIndex < nHeight; ++jIndex)
                {
                    int nIndex = iIndex * nHeight + jIndex;
                    image[nIndex] = Eigen::Vector3d::Zero();
                    energy[nIndex] = 255.0;
                }
            }
        }
    }

    // Calculate the Difference in Dimensions.
    int differenceInDimensionX = abs(width - currentWidth);
    int differenceInDimensionY = abs(height - currentHeight);

    // Create a while() Loop in order to loop until getting to the desired size.
    while (differenceInDimensionX != 0 || differenceInDimensionY != 0)
    {
        std::cout << differenceInDimensionX << " " << differenceInDimensionY << std::endl;

        // Determine the direction.
        Direction direction;
        int sign;
        if (differenceInDimensionX > differenceInDimensionY)
        {
            // If the width needs to be modified, then a Vertical Seam is required.
            direction = VERTICAL;
            differenceInDimensionX--;
            sign = width > currentWidth ? 1 : -1;
        }
        else
        {

	    // If the height needs to be modified, then a Horizontal Seam is required.
            direction = HORTIZONTAL;
            differenceInDimensionY--;
            sign = height > currentHeight ? 1 : -1;
        }

        // If the sign is negative, then remove the Seam.
        if (sign < 0)
        {
            removeSeam(direction);
            continue;
        }

        // If the sign is not negative, then add the Seam.
        addSeam(direction);
    }
}

// The Public Member Function "Eigen::Vector3d Carver::getPixel(int x, int y) const" that is an Accessor Function, or in other words the Getter Function, is defined
// for retrieving the Pixel information.
Eigen::Vector3d Carver::getPixel(int x, int y) const
{
    return image[x * originalHeight + y];
}

// The Public Member Function "void Carver::setPixel(int x, int y, Eigen::Vector3d color)" that is a Mutator Function, or in other words the Setter
// Function, is defined for applying a New Color Value to the Pixel.
void Carver::setPixel(int x, int y, Eigen::Vector3d color)
{
    image[x * originalHeight + y] = color;
}

// The Public Member Function "double Carver::getEnergy(int x, int y) const" that is an Accessor Function, or in other words the Getter Function, is defined for
// retrieving the energy information.
double Carver::getEnergy(int x, int y) const
{
    return energy[x * originalHeight + y];
}

// The Public Member Function "void Carver::setEnergy(int x, int y, double strength)" that is a Mutator Function, or in other words the Setter Function, is defined
// for applying a new energy value to the pixel.
void Carver::setEnergy(int x, int y, double strength)
{
    energy[x * originalHeight + y] = strength;
}

// The Private Member Function "void Carver::calculateEnergy()" calculates the initial energy value for the image.
void Carver::calculateEnergy()
{
    // Convert to Smoothed Normalization of the Pixel Data.
    for (int iIndex = 0; iIndex < currentWidth; iIndex++)
    {
        for (int jIndex = 0; jIndex < currentHeight; jIndex++)
        {
            Eigen::Vector3d pixel = getPixel(iIndex, jIndex);
            double energyAtPixel = pixel.norm();
            int countOfPixels = 0;
            if (iIndex > 0)
            {
                if (jIndex > 0)
                {
                    energyAtPixel += (getPixel((iIndex - 1), (jIndex - 1)) + pixel).norm();
                    countOfPixels++;
                }
                energyAtPixel += (getPixel((iIndex - 1), jIndex) + pixel).norm();
                countOfPixels++;
                if (jIndex < currentHeight - 1)
                {
                    energyAtPixel += (getPixel((iIndex - 1), jIndex + 1) + pixel).norm();
                    countOfPixels++;
                }
            }
            if (iIndex < currentWidth - 1)
            {
                if (jIndex > 0)
                {
                    energyAtPixel += (getPixel((iIndex + 1), (jIndex - 1)) + pixel).norm();
                    countOfPixels++;
                }
                energyAtPixel += (getPixel((iIndex + 1), jIndex) + pixel).norm();
                countOfPixels++;
                if (jIndex < currentHeight - 1)
                {
                    energyAtPixel += (getPixel((iIndex + 1), (jIndex + 1)) + pixel).norm();
                    countOfPixels++;
                }
            }
            if (jIndex > 0)
            {
                energyAtPixel += (getPixel(iIndex, (jIndex - 1)) + pixel).norm();
                countOfPixels++;
            }
            if (jIndex < currentHeight - 1)
            {
                energyAtPixel += (getPixel(iIndex, (jIndex + 1)) + pixel).norm();
                countOfPixels++;
            }
            setEnergy(iIndex, jIndex, energyAtPixel / countOfPixels);
        }
    }

    // Apply the Sobel Operator.
    /**
     * [-3  0  3] [ 3  10  3]
     * [-10 0 10] [ 0   0  0]
     * [-3  0  3] [-3 -10 -3]
     */
    // Implement the Scharr Gradient Operator.
    double *horizontal = convolve(scharrHorizontal);
    double *vertical = convolve(scharrVertical);

    // Create a for() Loop to loop through and apply the halving function.
    for (int iIndex = 0; iIndex < currentWidth; ++iIndex)
    {
        for (int jIndex = 0; jIndex < currentHeight; ++jIndex)
        {
            int index = iIndex * currentHeight + jIndex;

            // Get the Absolute Value from the Horizontal Convolution Kernel and Vertical Convolution Kernel.
            double hIsForHorizontalConvolutionKernel = abs(horizontal[index]);
            double vIsForVerticalConvolutionKernel = abs(vertical[index]);

            // Compute the Average of the values for the Horizontal Convolution Kernel and Vertical Convolution Kernel.
            energy[index] = (hIsForHorizontalConvolutionKernel + vIsForVerticalConvolutionKernel) * 0.5;
        }
    }

    // Clean up the allocated data.
    delete[] horizontal;
    delete[] vertical;
}

// The Private Member Function "double *Carver::calculateFlow(const Direction direction)" calculates the energy's flow to help determine the lowest energy path.
double *Carver::calculateFlow(const Direction direction)
{
    // Create the Flow Buffer that will hold the information.
    double *flowBuffer = new double[originalWidth * originalHeight];

    for (int iIndex = 0; iIndex < (originalWidth * originalHeight); ++iIndex)
    {
        flowBuffer[iIndex] = std::numeric_limits<int>::max();
    }

    // Calculate the flow differently in the case of looking for a Vertical Seam or a Horizontal Seam.
    if (direction == VERTICAL)
    {
        // Set the initial energy value of the bottom edge.
        for (int iIndex = 0; iIndex < currentWidth; ++iIndex)
        {
            int index = iIndex * originalHeight + (currentHeight - 1);
            flowBuffer[index] = energy[index];
        }

        // Create a for() Loop in order to loop from the bottom to the top edge.
        for (int xValueColumn = 0; xValueColumn < currentWidth; ++xValueColumn)
        {
            for (int y = currentHeight - 2; y >= 0; --y)
            {
                int index = xValueColumn * originalHeight + y;
                int indexA = std::max(0, xValueColumn - 1) * originalHeight + y + 1;
                int indexB = xValueColumn * originalHeight + y + 1;
                int indexC = std::min((int)currentWidth - 1, xValueColumn + 1) * originalHeight + y + 1;

                // Grab the three energy values associated with this Pixel.
                double aEnergyValue = flowBuffer[indexA];
                double bEnergyValue = flowBuffer[indexB];
                double cEnergyValue = flowBuffer[indexC];

                // Add the minimum energy to indicate the path to follow.
                flowBuffer[index] = energy[index] + std::min(aEnergyValue, std::min(bEnergyValue, cEnergyValue));
            }
        }
    }
    else
    {
        // Set the initial energy value of the right edge.
        for (int iIndex = 0; iIndex < currentHeight; ++iIndex)
        {
            int index = (currentWidth - 1) * originalHeight + iIndex;
            flowBuffer[index] = energy[index];
        }

        // Create for() Loop to loop from the right edge to the left edge.
        int end = currentHeight - 1;
        for (int xValueColumn = currentWidth - 2; xValueColumn >= 0; --xValueColumn)
        {
            for (int yValueRow = 0; yValueRow < currentHeight; ++yValueRow)
            {
                int index = xValueColumn * originalHeight + yValueRow;
                int indexA = (xValueColumn + 1) * originalHeight + std::max(yValueRow - 1, 0);
                int indexB = (xValueColumn + 1) * originalHeight + yValueRow;
                int indexC = (xValueColumn + 1) * originalHeight + std::min(yValueRow + 1, end);

                // Grab the three energy values associated with this Pixel.
                double aEnergyValue = flowBuffer[indexA];
                double bEnergyValue = flowBuffer[indexB];
                double cEnergyValue = flowBuffer[indexC];

                // Add the minimum energy to indicate the path to follow.
                flowBuffer[index] = energy[index] + std::min(aEnergyValue, std::min(bEnergyValue, cEnergyValue));
            }
        }
    }

    // Return the Flow Buffer.
    return flowBuffer;
}

// The Private Member Function "std::vector<int> Carver::calculateSeam(const Direction direction)" uses the Flow Buffer to calculate an appropriate seam.
std::vector<int> Carver::calculateSeam(const Direction direction)
{
    // Generate Flow Buffer.
    double *flowBuffer = calculateFlow(direction);

    // Start a vector to hold the seam
    std::vector<int> path;

    // If we're generating a vertical seam move top to bottom
    if (direction == VERTICAL)
    {
        path.reserve(currentHeight);

        // Get the Minimum Value along the top edge.
        int start = 0;
        double minimumValue = std::numeric_limits<double>::max();
        int cHeight = currentHeight - 1;
        for (int iIndex = 0; iIndex < currentWidth; ++iIndex)
        {
            double energyFlowAtPixel = flowBuffer[iIndex * originalHeight];
            if (energyFlowAtPixel < minimumValue)
            {
                minimumValue = energyFlowAtPixel;
                start = iIndex;
            }
        }

        // Start from the top edge and move towards the bottom edge by following the Minimum Flow Value.
        path.push_back(start);
        int jIndex = 1;
        int iIndex = start;
        int nextLevel = 0;
        while (path.size() < currentHeight)
        {
            minimumValue = std::numeric_limits<double>::max();
            nextLevel = iIndex;
            for (int kIndex = -1; kIndex < 2; ++kIndex)
            {
                // Find the Smallest Value in the Next Level.
                if ((kIndex + iIndex) < 0 || (kIndex + iIndex) >= currentWidth)
                {
                    // If it is off the edge, then skip it.
                    continue;
                }

                // If the energy is the smallest for the next value, then choose that energy that is the smallest for the Next Value.
                double energyFlowAtPixel = flowBuffer[(kIndex + iIndex) * originalHeight + jIndex];
                if (energyFlowAtPixel < minimumValue)
                {
                    minimumValue = energyFlowAtPixel;
                    nextLevel = (kIndex + iIndex);
                }
            }

            // Add the Seam Value and then continue to the Next Pixel.
            path.push_back(nextLevel);
            jIndex++;
            iIndex = nextLevel;
        }
    }
    else
    {
        path.reserve(currentWidth);

        // Get the Minimum Value along the left edge.
        int start = 0;
        double minimumValue = std::numeric_limits<double>::max();
        for (int jIndex = 0; jIndex < currentHeight; ++jIndex)
        {
            double energyFlowAtPixel = flowBuffer[jIndex];
            if (energyFlowAtPixel < minimumValue)
            {
                minimumValue = energyFlowAtPixel;
                start = jIndex;
            }
        }

        // Begin from the left edge. Then, move towards the right edge following the Ninimum Flow Value.
        path.push_back(start);
        int jIndex = start;
        int iIndex = 1;
        int nextLevel = 0;
        while (path.size() < currentWidth)
        {
            minimumValue = std::numeric_limits<double>::max();
            nextLevel = jIndex;
            for (int kIndex = -1; kIndex < 2; ++kIndex)
            {
                // Find the Smallest Value in the Next Level.
                if ((kIndex + jIndex) < 0 || (kIndex + jIndex) >= currentHeight)
                {
                    // If it is off the edge, then skip it.
                    continue;
                }

                // If the energy is the smallest for the Next Value, then choose that energy that is the smallest for the Next Value.
                double energyFlowAtPixel = flowBuffer[iIndex * originalHeight + (kIndex + jIndex)];
                if (energyFlowAtPixel < minimumValue)
                {
                    minimumValue = energyFlowAtPixel;
                    nextLevel = (kIndex + jIndex);
                }
            }

            // Add the Seam Value. Then, continue to the Next Pixel.
            path.push_back(nextLevel);
            iIndex++;
            jIndex = nextLevel;
        }
    }

    // Clear the memory of the Flow Buffer.
    delete[] flowBuffer;
    flowBuffer = nullptr;

    // Return the Seam.
    return path;
}

// The Private Member Function "void Carver::removeSeam(const Direction direction)" removes the calculated Seam from the Image and Energy Buffer.
void Carver::removeSeam(const Direction direction)
{
    // Calculate the Seam.
    std::vector<int> seam = calculateSeam(direction);

    // Create an if/else Statement that defines the condition for removing a Vertical Seam.
    if (direction == VERTICAL)
    {
        // Create a for() Loop to loop through the Size of the Seam.
        for (int yValueRow = 0; yValueRow < seam.size(); ++yValueRow)
        {
            // Get the x value that is represented by the Integer Variable "int xValueColumn" for each seam's y value that is represented by the Integer
            // Variable "int yValueRow".
            int xValueColumn = seam[yValueRow];

            // Create a for() Loop to Loop from the Seam Pixel to the right edge shifting everything left by 1.
            for (int iIndex = xValueColumn; iIndex < currentWidth - 1; ++iIndex)
            {
                setPixel(iIndex, yValueRow, getPixel(iIndex + 1, yValueRow));
                setEnergy(iIndex, yValueRow, getEnergy(iIndex + 1, yValueRow));
            }
        }

        // Update the Width.
        currentWidth--;
    }
    else
    {
        // Create a for() Loop to loop through the Size of the Seam.
        for (int xValueColumn = 0; xValueColumn < seam.size(); ++xValueColumn)
        {
            // Get the y value that is represented by the Integer Variable "yValueRow" for each seam's x value that is represented by the Integer
	    // Variable "int xValueColumn".
            int yValueRow = seam[xValueColumn];

            // Create a for() Loop to loop from the Seam Pixel to the bottom edge shifting everything up by 1.
            for (int jIndex = yValueRow; jIndex < currentHeight - 1; ++jIndex)
            {
                setPixel(xValueColumn, jIndex, getPixel(xValueColumn, jIndex + 1));
                setEnergy(xValueColumn, jIndex, getEnergy(xValueColumn, jIndex + 1));
            }
        }

        // Update our Height.
        currentHeight--;
    }
}

// The Private Member Function "void Carver::addSeam(const Directoin direction)" adds the Calculated Seam from the Image and Energy Buffer.
void Carver::addSeam(const Direction direction)
{
    // Calculate the Seam.
    std::vector<int> seam = calculateSeam(direction);

    // Create an if/else statement that defines the condition for adding a Vertical Seam.
    if (direction == VERTICAL)
    {
        // Modify the Width.
        currentWidth++;

        // Create a for() Loop to loop through the Size of the Seam.
        for (int yValueRow = 0; yValueRow < seam.size(); ++yValueRow)
        {
            // Get the x value that is represented by the Integer Variable "int xValueColumn" for each seam's y value that is represented by the Integer
	    // Variable "int yValueRow".
            int xValueColumn = seam[yValueRow];

            // Create a for() Loop to loop from the right edge to the Seam Pixel shifting everything right by 1.
            for (int iIndex = currentWidth - 1; iIndex > xValueColumn; --iIndex)
            {
                setPixel(iIndex, yValueRow, getPixel(std::max(0, iIndex - 1), yValueRow));
                setEnergy(iIndex, yValueRow, getEnergy(std::max(0, iIndex - 1), yValueRow));
            }

            // Determine the Pixels to the left of the Seam and to the right of the Seam.
            int lLeftX = std::max(xValueColumn - 1, 0);
            int rRightX = std::min(xValueColumn + 1, currentWidth - 1);

            // Retrieve the image and energy values.
            Eigen::Vector3d left = getPixel(lLeftX, yValueRow);
            Eigen::Vector3d right = getPixel(rRightX, yValueRow);
            double leftEnergy = getEnergy(lLeftX, yValueRow);
            double rightEnergy = getEnergy(rRightX, yValueRow);

            // Set the Pixel to an Average of the Two Pixels around it.
            setPixel(xValueColumn, yValueRow, (left + right) / 2);

            // Set the Three Pixels to the combined energy to make the Seam unlikely to be picked multiple times in succession.
            setEnergy(lLeftX, yValueRow, (leftEnergy + rightEnergy));
            setEnergy(xValueColumn, yValueRow, (leftEnergy + rightEnergy));
            setEnergy(rRightX, yValueRow, (leftEnergy + rightEnergy));
        }
    }
    else
    {
        // Modify the Height.
        currentHeight++;

        // Create a for() Loop to loop through the Size of the Seam.
        for (int xValueColumn = 0; xValueColumn < seam.size(); ++xValueColumn)
        {
            // Get the y value that is representd by the Integer Variable "int yValueRow" for each seam's x value that is representedd by the Integer
	    // Variable "int xValueColumn".
            int yValueRow = seam[xValueColumn];

            // Create a for() Loop to loop from the bottom edge to the Seam Pixel shifting everything down by 1.
            for (int jIndex = currentHeight - 1; jIndex > yValueRow; --jIndex)
            {
                setPixel(xValueColumn, jIndex, getPixel(xValueColumn, std::max(0, jIndex - 1)));
                setEnergy(xValueColumn, jIndex, getEnergy(xValueColumn, std::max(0, jIndex - 1)));
            }

            // Determine the Pixels above the Seam and below the Seam.
            int lLeftY = std::max(yValueRow - 1, 0);
            int rRightY = std::min(yValueRow + 1, currentHeight - 1);

            // Retrieve the image and energy values.
            Eigen::Vector3d left = getPixel(xValueColumn, lLeftY);
            Eigen::Vector3d right = getPixel(xValueColumn, rRightY);
            double leftEnergy = getEnergy(xValueColumn, lLeftY);
            double rightEnergy = getEnergy(xValueColumn, rRightY);

            // Set the Pixel to an average of the Two pixels around it.
            setPixel(xValueColumn, yValueRow, (left + right) / 2);

            // Set the Three Pixels to the combined energy to make the Seam unlikly to be picked multiple times in succession.
            setEnergy(xValueColumn, lLeftY, (leftEnergy + rightEnergy));
            setEnergy(xValueColumn, yValueRow, (leftEnergy + rightEnergy));
            setEnergy(xValueColumn, rRightY, (leftEnergy + rightEnergy));
        }
    }
}

// The Private Member Function "double *Carver::convolve(Eigen::Matrix3d kernel)" performs a convolution on the Energy Buffer using the kernel provided.
double *Carver::convolve(Eigen::Matrix3d kernel)
{
    // Generate the Convolution Buffer.
    double *convolvedBuffer = new double[currentWidth * currentHeight];

    double sum;

    // Create a for() Loop to loop through the Width of the Image.
    for (int xValueColumn = 0; xValueColumn < currentWidth; ++xValueColumn)
    {

        // Create a for() Loop to loop through the Height of the Image.
        for (int yValueRow = 0; yValueRow < currentHeight; ++yValueRow)
        {
            // Next, generate the Pixel Index. Then, start the Convolution.
            int index = xValueColumn * currentHeight + yValueRow;
            sum = 0;

            // Create a for() Loop to loop through the 3x3 kernel, in other words loop through the kernel with 3 by 3 dimensions, and apply to the Pixels
            // surrounding the current one.
            for (int iIndex = -1; iIndex < 2; ++iIndex)
            {
                for (int jIndex = -1; jIndex < 2; ++jIndex)
                {
                    // If it is off the edge, then do the calculation on the Current Pixel.
                    if (xValueColumn + iIndex >= currentWidth || xValueColumn + iIndex < 0 || yValueRow + jIndex >= currentHeight || yValueRow + jIndex < 0)
                    {
                        sum += kernel.coeff(iIndex + 1, jIndex + 1) * getEnergy(xValueColumn, yValueRow);
                        continue;
                    }

                    // Add the calculation to the the sum.
                    sum += kernel.coeff(iIndex + 1, jIndex + 1) * getEnergy(xValueColumn + iIndex, yValueRow + jIndex);
                }
            }

            // Average the result of the Convolution.
            convolvedBuffer[index] = sum / 9;
        }
    }

    // Return the generated values.
    return convolvedBuffer;
}
