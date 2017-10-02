# WPG_bspline_-_DE

Developped for Matlab 2012a and more recent

WPG based on deboor bspline and a deformation estimator to take into account soft sole deformation detailed in https://www.researchgate.net/publication/318412580_Optimized_Humanoid_Walking_with_Soft_Soles


In main.m is the WPG built-in as a programming object in matlab.

Uncomment the two last lines to enable the Deformation estimator.

matlab compiler option for the deformation estimator:

mex '-IC:\Users\Giovanni\Documents\eigen-eigen-1306d75b4a21\eigen-eigen-1306d75b4a21' Simulation_3D\splineBasis.cpp

mex '-IC:\Users\Giovanni\Documents\eigen-eigen-1306d75b4a21\eigen-eigen-1306d75b4a21' GaussFtotZMP.cpp Gausstot.cpp

In \fsqp : 
mex mexfsqp.c cfsqp.c qld.c


WARNING : no interpolation of foot in the air are made in this program to generate the foot trajectory in single support phases when the foot in not in contact with the floor.

In test deboor bspline can be found a Mathematica method to generate the polynomials of deboor bspline of any order.