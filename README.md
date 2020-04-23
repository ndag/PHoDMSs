# Spatiotemporal Persistent Homology for Dynamic Data
The code in this repository works to fill out the computational pipeline of taking in a discrete dynamic metric space (DMS), generating a spatiotemporal Betti-0 function for it, and then compute erosion distance between these Betti-0 functions. This is based on theoretical work by Woojin Kim and Facundo Memoli of The Ohio State University Department of Mathematics. The details can be found in the paper: "Spatiotemporal Persistent Homology for Dynamic Metric
Spaces", published in Discrete & Computational Geometry, 2020, link: https://link.springer.com/article/10.1007/s00454-019-00168-w . 
The code is authored by Nate Clause (The Ohio State University Department of Mathematics) in collaboration with Woojin Kim.

## Spatiotemporal Betti-0 function generation
All of the codes are currently set up to be run from a command line, for example, the current approach to run "betti_generator.py" would be to enter command line or the terminal in a python IDE and write something such as:

```
python betti_generator.py [dmsfile] [bettifile] 40 0 50 5
```

All of these codes have portions at the top of main where the user can change the approach to how these functions are called/executed to another method as desired. Later, code will also be posted that exemplifies a larger-scale version with parallelization for generating Betti-0 functions as well as computing erosion distance between a large number of Betti-0 functions.

"boids_simulation.py" is code that uses a standard flocking boids model to generate DMSs, and is added to give the user an example of how to generate DMS, as well as to have easily-made DMS available to test with the other codes.
This code has two sets of initial parameters (again, this can be changed at the top of main):

```
python boids_simulation.py num_points separation_force alignment_force cohesion_force dmsfile
```

or:

```
python boids_simulation.py num_points separation_force separation_radius alignment_force alignment_radius cohesion_force cohesion_radius dmsfile
```

The forces and radii are parameters that alter the motion of the boids based on the rules based on the rules of the model. For a general explanation of the flocking boids model, see https://en.wikipedia.org/wiki/Boids. Then dmsfile is the filename that the DMS data is saved out to.

"betti_generator.py" is code that takes in a dmsfile and generates a Betti-0 function for that DMS. It currently has one set of initial parameters:

```
python betti_generator.py [dmsfile] [bettifile] num_points start_threshold end_threshold spacing
```

"dmsfile" and "bettifile" are the names of the files containing the input DMS and the output Betti-0 function, respectively. num_points is the number of points in the DMS (will soon be unecessary to input manually). The threshold and spacing inputs determine the discretization values to be used for the Vietoris-Rips filtration when generating the Betti-0 function. For example, with input:

```
python betti_generator.py dmsfile bettifile num_points 0 50 5
```

The VR distance thresholds for finding Betti-0 numbers would be taken as 0, 5, 10, ..., 50
The code is currently built so that the DMS file it reads in is a dynamic point cloud and then it computes the time series of distance matrices using Euclidean distance in the function "get_dist". This "get_dist" function is the one to change should the DMS be based on other metrics. In the future, additional functionality will be added so the user could input either a DMS as a dynamic point cloud or as a time series of distance matrices already. Of note is that the "boids_simulation.py" output is a dynamic point cloud, as the current "betti_generator.py" code is setup to handle.

