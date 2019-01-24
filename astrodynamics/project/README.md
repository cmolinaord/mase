# Computation of possible future launch windows to reach Jupiter from Earth using Poliastro

With this project, I want to develop an algorithm that could compute a flyby of a planet from the Earth, with the possibility of performing a gravity assist on it. As first try, the program was developed to depert from the Earth and flyby Jupiter, but the aim is to extend this capacity to the pair of planets desired.

## Earth_to_Jupiter
Calculates the manuever and the launch windows from Earth to Jupiter, assuming tangential impulsive maneuver from Earth orbit around the Sun

Earth_to_Jupiter(dv):

'dv' is given with velocity units (u.m/u.s)

Returns:
- ss_transfer = trasfer orbit object
- t = time to the first next launch window
- T = period between possible launch windows
- t_transfer = time of transfer from launch to Jupiter arrival
   
