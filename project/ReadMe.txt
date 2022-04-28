function format:

ipf opf #roi x1 y1 Sx1 Sy1 func p1 p2 pX x2 y2 Sx2 Sy2 func p1 p2 pX

ipf - input image file
opf - output image file
func - function
#roi - number of rois [1, 3]
x - col of first pixel in roi (left to right)
y - row of first pixel in roi (top to bottom)
sx -  total number of pixels in x axis
sy - total number of pixels in y axis


available functions:
add
binarize
scale
dualthres
reg2dsmooth
sep2dsmooth
colorbright
colorvisual