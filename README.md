# koopman-cbf

Code for collision avoidance experiments with CBFs synthesized using learned Koopman operators.

This repository contains 7 main scripts that run different parts of experiments demonstrating the Koopman + CBF framework. The first 4 scripts are concerned with learning an approximate Koopman operator for a Dubin's car model and then using the learned model to run experiments in the Georgia Tech Robotarium. These files are (all can be run separately by utilizing learned models that are stored in the repository): 
- 'dubin_learning.m' generates training data and learns a model
- 'robotarium_obstacle_avoidance.m' executes a simulated experiment at the Robotarium where a single agent avoids a fixed obstacle
- 'robotarium_collision_avoidance.m' executes a simulated experiment at the Robotarium where 5 agents avoid collision with each other
- 'robotarium_collision_obstacle_avoidance.m' executes a simulated experiment at the Robotarium where 5 agents avoid collision with each other and a fixed obstacle 

The 3 remaining scripts are concerned with learning an approximate Koopman operator for an UAV and simulating 2 different experiments utilizing the learned model. These files are (all can be run separately by utilizing learned models that are stored in the repository):
- 'uav_learning.m' generates training data and learns a model
- 'uav_collision_avoidance.m' executes a simple experiment where 2 UAVs avoid collision with each other
- 'uav_collision_ge_exp.m' executes an experiment where 3 UAVs attempt to land on the same landing pad while avoiding collision with each other and hard impacts with the ground
