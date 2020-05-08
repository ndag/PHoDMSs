# Spatiotemporal Persistent Homology for Dynamic Data

A **dynamic metric space (DMS)** is a time series of distance functions over a fixed underlying set. Instance of DMSs include collective behaviors of animals (a flock of birds or a school of fish), and social networks in the human sphere.

The code in this repository works to 
1. Compute the **spatiotemporal persistent Betti-0 function** of **DMSs**, and in turn
2. Quantify the dissimilarity between two DMSs, based on their **spatiotemporal persistent topology**. For this, we make use of a slight generalization of [the erosion distance](https://link.springer.com/article/10.1007/s41468-018-0012-6) by Amit Patel.

In addtion, using this code, users can generate various DMSs according to [Boids model](https://en.wikipedia.org/wiki/Boids). Find real time simulations at the bottom of [this webpage](https://research.math.osu.edu/networks/formigrams/).

<!--This code fill out the computational pipeline of taking in a dynamic metric space (DMS) generating a spatiotemporal Betti-0 function for it, and then compute erosion distance between these Betti-0 functions.--> 
The code is based on theoretical work by [Woojin Kim](https://wj-kim.com) and [Facundo Mémoli](https://people.math.osu.edu/memolitechera.1/) at Ohio State (see [*Spatiotemporal Persistent Homology for Dynamic Metric Spaces*](https://link.springer.com/article/10.1007/s00454-019-00168-w), published in Discrete & Computational Geometry, 2020).

<!--. The details can be found in the paper: [*Spatiotemporal Persistent Homology for Dynamic Metric Spaces*](https://link.springer.com/article/10.1007/s00454-019-00168-w), published in Discrete & Computational Geometry, 2020.-->
The code is authored by [Nate Clause](https://math.osu.edu/people/clause.15), a Math Ph.D. student at Ohio State, in collaboration with [Woojin Kim](https://wj-kim.com).

## Generating DMSs from Boids model
"boids_simulation.py" is code that uses a standard flocking boids model to generate DMSs, and is added to give the user an example of how to generate DMS, as well as to have easily-made DMS available to test with the other codes.
This code has two sets of initial parameters (again, this can be changed at the top of main):

```
python boids_simulation.py [num_points] [separation_force] [alignment_force] [cohesion_force] [dmsfile]
```

or:

```
python boids_simulation.py [num_points] [separation_force] [separation_radius] [alignment_force] [alignment_radius] [cohesion_force] [cohesion_radius] [dmsfile]
```

The forces and radii are parameters that alter the motion of the boids based on the rules based on the rules of the model. For a general explanation of the flocking boids model, see the [Wikipedia](https://en.wikipedia.org/wiki/Boids). Then [dmsfile] is the filename that the DMS data is saved out to.

The motion modeled in this scenario occurs on a 2-dimensional torus of width 500 units and height 250 units. The parameters should have the following ranges:

0 < [num_points] < infinity

0 <= [separation_force] < 10

0 <= [separation_radius] < infinity

0 <= [alignment_force] < 10

0 <= [alignment_radius] < infinity

0 <= [cohesion_force] < 10

0 <= [cohesion_radius] < infinity

We give a few remarks here. First, the number of points is the main factor deciding the runtime, so while we can have arbitrarily many points, the program can become very slow if we let the number of points go too high. Anything at 100 or fewer points should run just fine, more powerful infrastructures can handle hundreds of points (though for visual lag reasons, with more than 100 points, it is better to not show the plot as the program runs but instead to view the saved mp4 afterwards). 

Next, the force coefficients can be 0, this simply means that force is playing no role in the simulation. Technically, these coefficients could be input arbitrarily large, but making them more than 10 can make the simulation quite unstable. It is best overall to keep all of these coefficients less than 5 for "smoothness" of the output.

Lastly, all the radius parameters can be arbitrarily large, but there is no difference between having a radius parameter equal to the radius of the ambient space of motion and any larger parameter. With our current 500x250 window, this means making one of these parameters more than ~280 will yield no change from if it were 280. 


## Computing Spatiotemporal Betti-0 function
<!--All of the codes are currently set up to be run from a command line, for example, the current approach to run "betti_generator.py" would be to enter command line or the terminal in a python IDE and write something such as:
```
python betti_generator.py [dmsfile] [bettifile] 40 0 50 5
```
All of these codes have portions at the top of main where the user can change the approach to how these functions are called/executed to another method as desired. Later, code will also be posted that exemplifies a larger-scale version with parallelization for generating Betti-0 functions as well as computing erosion distance between a large number of Betti-0 functions. -->

"betti_generator.py" is code that takes in a DMS file and generates a Betti-0 function for that DMS. It currently has one set of initial parameters:

```
python betti_generator.py [dmsfile] [bettifile] [start_threshold] [end_threshold] [spacing] [time_samples]
```

[dmsfile] and [bettifile] are the names of the files containing the input DMS and the output Betti-0 function, respectively. The threshold and spacing inputs determine the discretization values to be used for the Vietoris-Rips filtration when generating the Betti-0 function. Consider the following example input:

```
python betti_generator.py DMS1.txt Betti.txt 0 50 5 20
```

Then, for computing the Betti-0 function of **DMS1.txt**, the Vietoris-Rips complexes will be computed at the thresholds  0 (=**start_threshold**), 5, 10, ..., 50 (=**end_threshold**) and the Betti-0 function will be saved into a file named **Betti.txt**. The Betti-0 function will use 20 evenly spaced time samples from the DMS. i.e. if the DMS consists of 100 points moving through 1000 timesteps, the Betti-0 function will be based on 20 evenly spaced timesteps. The boids simulation currently is setup to use 1000 timesteps, but creating a Betti-0 function based on 1000 timesteps will be very slow, and computing the erosion distance between to Betti-0 functions of this size is not computationally feasible. The code is currently built so that the DMS file it reads in is a dynamic point cloud and then it computes the time series of distance matrices using Euclidean distance in the function "get_dist". This "get_dist" function is the one to change should the DMS be based on other metrics. In the future, additional functionality will be added so the user could input either a DMS as a dynamic point cloud or as a time series of distance matrices already. Of note is that the "boids_simulation.py" output is a dynamic point cloud, as the current "betti_generator.py" code is setup to handle.

## Computing Erosion distance between Betti-0 functions. 
There are currently two files codes which can be used to perform this task, "erosion_distance.py" and "same_size_erosion_distance.py". These will soon be merged into a single file. At this time, if the user has wants to compute erosion distance between two betti-0 functions of the same size (say, both 10 by 20 by 20 arrays), we recommend using "same_size_erosion_distance.py". The code in "erosion_distance.py" can break if the thresholds for the Rips parameter in generating the betti-0 functions did not go high enough so that the "top-level" of the betti-0 functions are filled with only 1s, whereas the recommended code can deal with this situation, but currently only when the betti-0 functions are of the same size. We expect by the end of April 30th, 2020 that this issue will be resolved, and these two erosion distance computation codes will be merged back together into just "erosion_distance.py", and with some further runtime optimizations implemented.

"same_size_erosion_distance.py" is code that takes in two betti-0 function files and computes the erosion distance between them. It currently has one set of initial parameters:

```
python same_size_erosion_distance.py [bettifile1] [bettifile2] 
```
With this input, the code would compute the erosion distance of the two betti-0 functions and display this result.
