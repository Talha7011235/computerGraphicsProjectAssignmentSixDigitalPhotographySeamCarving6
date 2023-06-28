Talha Hussain, myDigitalPhotographySeamCarvingProjectSix

Computer Graphics Project Assignment 6 Digital Photography Seam Carving

Websites:
To complete my Computer Graphics Project Assignment 6 Digital Photography and Seam Carving, I received helpful inspiration from the following website:
https://towardsdatascience.com/computer-vision-convolution-basics-2d0ae3b79346 - I used this website article "Computer Vision: Convolution Basics" by Harsh Yadav
to get inspiration in applying Convolution to my "Computer Graphics Project Assignment 6 Digital Photography Seam Carving".
https://en.wikipedia.org/wiki/Sobel_operator#Alternative_operators - I used this website to get inspiration and a better understanding of the Scharr Gradient
Operator in order to come up with my own idea of implementing the Scharr Gradient Operator to achieve one of the Extra Credit options of "Computer Graphics Project
Assignment 6 Digital Photography Seam Carving". 
https://www.todocanada.ca/25-summer-activities-adventures-in-edmontons-river-valley/ - I used a random .jpg image of a River from this website and converted the .jpg
image of the River to a Portable Network Graphics .png image to show that loading and saving of .png Files worked with my Seam Carving program.

Usage:
make seamcarving - Runs the compilation of the primary built target that is the "seamcarving" application.
make clean - cleans the intermediate files and compiled target

After running the command "make seamcarving" to build the target that is the "seamcarving" application, the project can be run using the following command rule to
test any single image of one's choosing along with changing the dimensions of that single image in terms of resizing the width of the image and the height of the
image: 
./seamcarving input_image output_image output_width output_height

Run the following command to test some image commands that have been included and built in the "Makefile" for my Computer Graphics Project Assignment 6 Digital
Photography Seam Carving to test the various adding and removing of seams:
make test

Summary:
For this assignment, I started with the energy function and built a basic convolution system using Sobel operators in order to make adding the Sharr operators easier
to implement when I got to the extra credit section. This energy function starts by doing a simple 9 pixel smoothing to determine the energy. I then use that
generated value to perform a horizontal and vertical convolution with the 9x9 operators. The results of that convolution is then averaged into the final energy
of the pixel. Once the energy is calculated, I am then able to start the resize process. I begin with checking to see if a dimension needs to be expanded. If a
dimension that is either the width or the height needs to be expanded, then I generate a New Image Buffer with the largest bounds to properly hold the data. After
that, a while loop is used to ensure that if the width and height do not match the desired resize value, then the carving would continue. Inside the carving loop, I
check to see which dimension is the most distant in order to determine if I need to do a vertical carving or a horizontal carving. When the direction of the
carving is determined, I start by generating a Flow Buffer and that is a cumulative energy value from the right or bottom edge to its opposite. This final value
tells me where to start when generating a seam with the smallest energy value from start to end. With that, I can create a seam vector by following that path back
to the source. With the seam generated, I am then able to loop through the [x,y] coordinates of the seam to find the pixels that need to be added or removed. In
removing pixels in the seam, my C++ code for my "Computer Graphics Project Assignment 6 Digital Photography Seam Carving" shifts the values of the pixel one space
from the desired edge. In addition, the energy is shifted and this results in keeping the Energy Buffer valid for the image. Furthermore, this ensures that I do not
need to recalculate the energy in every frame. Instead, I just need to modify it in the same way that I modify the image. Adding along the seam starts out the same as
removing pixels in the seam just in the opposite direction that means in other words, moving the pixels away from the seam and leaving a duplicated value. I then
set the new pixels to be an average of the two pixels on its sides. In addition, the energy is modified, however, the energy is modified in a different way
because the seam's pixels and the pixels on either side of the seam are changed to be the sum of the energies of the two side pixels. This makes it unlikely for the
seam or the pixels on either side to be reused for a new seam immediately and that stops the image from being constantly stretched using the same pixel information
and creating a smear. After these functions were complete, my C++ code for my "Computer Graphics Project Assignment 6 Digital Photography Seam Carving" then operates
to continue until the final resolution was achieved and then saving the resulting image. So far I have found the result to be quite pleasing so long as the images are
not stretched or shrunk too far which will still cause artifacts to appear. Finally, I substituted the Scharr Gradient Operator values in place of the Sobel operator
and created a simple HTML Web Page to showcase two examples using two of the images that applies the concepts of Seam Carving
for my "Computer Graphics Project Assignment 6 Digital Photography Seam Carving".

