README

Required files are Tubeplot.py , first_der.py , and seconder_der.py.
DriveTubePlot.py is a test function which provides inputs to tubeplot.

This code requires very few inputs in order to make a nice tube plot.
The MATLAB version of this code is much nicer looking than the matplotlib python version. 
first_der and second_der files are required as they are the finite element 1st and 2nd derivatives used in creating the shape. 

**Inputs for filament:**
- r(N,3) = Euclidian location of N number of points which define a backbone of a filament in 3 Euclidian dimensions
- a(N,1) = radius of each node. This assumes each node is a cylinder. 
- Omega3(N,1) = $\Omega_3$ tangential component of a strain vector (twist / unit length). Using zeros everywhere assumes that your shapes is bent into its conformation

**Inputs for plotting (Python):**
- opacity = zero to 1 representation of opacity
- TubeColor = [r,g,b] color of the tube
- StripeColor = [r,g,b] color of the stripe painted on the tube
- RadPoints = number of points which represent the circle at radius = a. More points here makes the tube appear to have a more circular cross section
