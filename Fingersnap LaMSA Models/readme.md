# Finger snap LaMSA model

## Summary

The goal of this model is to capture the trends related to the effects of skin friction and compression on maximum velocity which were found expeirmentally. To see more details on the results of these experiments, as well as the development of this model, please refer back to the main paper. 

There are two major functions that can be used to run simulations on this model: simulation.m and fricVar.m. The first function allows the user to provide details the latch, load, spring, motor, and friction capabilities of a system of interest and receive in return details on the kinematic and dynamic evolution of the system with respect to time, as well as key information such as unlatch and take off times, positions, and velocities. To assist in the analysis of friction effects, we developed and provide the function fricVar.m, which takes similar inputs to simulation.m and calls it several times for many friction coefficients. This function creates four figures detailing the effects of friction and returns a single structure containing key kinematic and dynamic information for each coefficient of friction tested. Additional detail on these functions can be found in the Usage section. 

## Required Files and Setup

To ensure proper functioning of all parts of this model, the following functions must be within the same folder:

- Ff.m
- Fl.m
- Fs.m
- yLoad.m
- dyLoad.m
- d2yLoad.m
- latchODE.m
- findLoad.m
- findtUnlatch.m
- simulation.m
- fricVar.m

All the above files are written in MATLAB version R2020b. 

## Usage

As explained in the summary section, two main functions are provided to run the simulation. The first is titled simulation.m. By inputing information regarding the system, the function will return detailed information on the evolution of the system over time as well as other key points of information. The inputs and outputs are detailed below. The function can be called using a command as follows:

```sh
[t,y,dy,d2y,x,dx,d2x,Fspring,netF,unlatchTime,unlatchTimeError,loadUnlatch,timeTO,loadTO,s,Wf,N] = simulation(loadStr,latchStr,startStr,fricStr,detailsStr)
```

In this function call, the inputs are defined as follows:

- startStr: structure with info about t =0. fields:
    - .startTheta (double): initial angle between latch and load in degrees
    - .startVx (double): intial horizontal velocity of the latch

- loadStr: structure with info about spring and loading force. Fields:
    - .k (double):  spring constant
    - .m (double): mass of the load

- fricStr: structure with info about friction values/model. Fields:
    - .form (vector): contains info for the nonlinear friction model [F_f = mu*N(N/N_0)^n] Vector has:
        - n (double): exponential factor N is raised to. Must be -1<n<0
        - N_0 (double): N/N_0. Usually 1.
    - .mus (double): static coefficient of friction.
    - .muk (double): dynamic coefficient of friction. 

- latchStr: structure with info about the latch. Fields:
    - .R: double - radius of the circular latch
    - .motor: cell - information on unlatching motor. Contains:
        - Motor type: string - "linear_motor", "constant_force". Additional Information:
            - If motor type is "linear_motor", cell also needs:
                - Motor Definitions (vector): contains information about the motor including: 
                    - max_Force (double): maximum motor force
                    - max_velocity (double): maximum velocity of unlatching motor
                    - range_of_motion (double): motor range
                - Braking (boolean): If true, motor stops at end of range of motion. If false, motor continues regardless of range of motion

- detailsStr: structure with additional information. Fields:
    - .overShoot (double): for each mu, percentage of time after t_t0 to continue calculating position, velocity, acceleration, and force data.


simulation.m outputs the following:

- t (n x 1 vector): Each time point for which a position, velocity, etc. is calculated.
- y (n x 1 vector): The vertical position of the load at each time in vector t.
- dy (n x 1 vector): The velocity of the load at each time in vector t.
- d2y (n x 1 vector): The acceleration of the load at each time in vector t.
- x (n x 1 vector): The horizontal position of the latch at each time in vector t.
- dx (n x 1 vector): The velocity of the latch at each time in vector t.
- d2x (n x 1 vector): The acceleration of the latch at each time in vector t.
- Fspring (n x 1 vector): The force the spring is exerting on the load at each time in vector t. 
- netF (n x 1 vector): The net force on the load at each time in vector t. 
- unlatchTime (double): The calculated time at which the load and latch are no longer in contact; when netF = Fspring.
- unlatchTimeError (double): The error associated with the calculated time. 
- loadUnlatch (2 x 1 vector): The vertical position and velocity of the load at the unlatch time. 
- timeTO (double): The calculated time at which the load disengages from the spring and has no force acting on it; when Fspring = 0. 
- loadTO (2 x 1 vector): The vertical position and velocity of the laod at the take off time. 
- s (m x 1 vector): The arc length travelled by the load on the circular region of the latch until the unlatch time.
- Wf (m x 1 vector): The work done by friction between latch and load until the unlatch time. 
- N (m x 1 vector): The normal force between the latch and load until the unlatch time. 

To better test the variation caused by friction, we provide a second function entitled fricVar.m. This function runs the simulation.m function for several provided values of friction coefficient and provides the results and graphical representations of the results. The function can be called as follows:

```sh
[kineticsa] = fricVar(startStr,loadStr,fricStr,latchStr,detailsStr)
```

The inputs startStr, loadStr, and latchStr follow the same structure as those required by the simulation.m function. The fricStr, loadStr, and detailsStr inputs however, are now defined as:

- fricStr: structure with info about friction values/model. Fields:
    - .form (vector): contains info for the nonlinear friction model [F_f = mu*N(N/N_0)^n] Vector has:
        - n (double): exponential factor N is raised to. Must be -1<n<0
        - N_0 (double): N/N_0. Usually 1.
    - .mus (vector): contains all mu_s which should be tested.
    - .mu_rat (double): ratio between mu_s and mu_k.
- loadStr: structure with info about spring and loading force. Fields:
    - .k (double):  spring constant
    - .loading (string): type of loading - "nonlinear" or "phenomenological"
    - .loadingInfo (cell): additional information
        - If loading is phenomenological, cell requires:
            - alpha (double): linear relationship between loaded force and mu
            - F0 (double): loaded force at mu=0
        - If loading is nonlinear, cell can be empty.
    - .loadMaxF (double): maximum force spring can be loaded to 
    - .m (double): mass of the load
- detailsStr: structure with additional information. Fields:
    - .overShoot (double): for each mu, percentage of time after t_t0 to continue calculating position, velocity, acceleration, and force data.
    - .backVec (vector): indexes of mu_vec to graph on the 'Kinetics' figure.       

This function then performs the following:
1. Creates a folder within the directory numbered based on how many times fricVar.m has been run.
2. Creates a txt file within this folder with a summary of the parameters tested.
3. Runs simulation.m with every provided friction coefficient.
4. Creates and saves the following four figures to the created folder:
    a. Kinetics: Selects three friction coefficients as well as those specified in detailsStr.backVec and plots the vertical position and velocity of the load and the normal force between latch and load for each. 
    b. Energetic Evolution: Selects two friction coefficients and plots the evolution of spring potential energy into kinetic energy, dissipated friction energy, and energy dissipated by the structure for each coefficient side by side.
    c. Trends: Plots the initial spring force, unlatch time, and take off velocity found for each friction coefficient. 
    d. Energetics: Plots the relationship between the initial spring potential energy and the maximum kinetic energy. Also plots the derivatives of both energy values with respect to the friction coefficient. A smaller axis is also created within this plot which focuses on the friction coefficients at which the crossover in energy derivatives occurs. 
5. Returns to the original folder.

## Function Details

To assist in understanding, use, and improvement of this model, here we will discuss in greater detail the individual functions developed for this model. 

_yLoad.m_
This function takes as inputs the latch radius and the current horizontal position of the latch. Using the equation of a circle, this function calculates the vertical position of the load if the load and the latch were still in contact. 

_dyLoad.m_
Building off the previous function, dyLoad takes the latch radius, horizontal position of the latch, and velocity of the latch as inputs. By taking the derivative of the equation used in function _yLoad.m_ with respect to horizontal latch position, the function calculates the expected vertical load velocity. 

_d2yLoad.m_
Similarily to _dyLoad.m_, by taking the derivative with respect to horizontal latch position of the equation used in _dyLoad.m_, we developed an equation which can calculate the expected vertical acceleration of the load from the latch radius and the horizontal position, velocity, and acceleration of the latch. 
One major shortcoming of this model is that the compressibility of the latch and load (finger pads) is accounted for in the friction model and is decoupled from the actual geometric changes of the system. These previous three functions, which dictate the relationship between the positions, velocities, and accelerations of the latch and load based on the geometry of these structures, is likely one major area where improvements can be made for future iterations of this model. 

_Fs.m_
This function takes a description of the spring and the current vertical position of the load to calculate the spring force acting on the load. Currently, the spring description is set up to contain the equilibrium vertical position of the spring and the spring constant. However, this can be adjusted simply by adding aditional variables to the springDescription cell array and adjusting the spring force equation in this function. 

_Fl.m_
This function takes in a description of the latch as well as the latch radius, the friction coefficient, and the current positions, velocities, and accelerations of the latch and load to calculate the force that the modelled unlatching motor exerts on the latch. Currently, only the latch description, latch position, and latch velocity variables are used but the other variables are included in the function definition to allow for modifications for different motor types or unlatching methods. Two types of unlatching are currently supported, which can be defined by the first cell in the latch description cell array: a constant force motor, and a linear force_velocity motor. Additional details for each motor type can be included in the latch description cell array, allowing this function to be modified as needed. 

_latchODE.m_
The _latchODE.m_ function calls on all previously defined functions to create a differential equation that can be solved using ordinary differential equation (ODE) solvers. By combining equations 3.1, 3.3, 3.4, and 3.6 described in the main paper, a differential equation can be found which determines the latch acceleration from the inputs: the time, latch position and velocity, friction coefficients and properties, latch and load mass, latch and spring properties, and the calculated unlatch time. As in most cases, the friction model used is not linear with respect to the normal force, the calculation of the latch acceleration requires the nonlinear solving of the used equations. Additionally, we have implemented checks which will throw errors if the function detects the latch getting stuck during unlatching. 

_find_tUnlatch.m_
For performance reasons, it is preferable to calculate an approximate unlatch time prior to solving the ODE. In order to do this, the ODE solver ode45 is used. An exit flag is specified for when the net force becomes equivalent to the spring force, which defines the unlatch time. This method allows for a rapid determination of an approximate unlatch time.

_simulation.m_
This function combines all previous functions to properly analyze the kinematics and dynamics of a given system. The inputs and outputs of this function are described in detail in the Usage section. First, the approximate unlatch time is calculated. This is then used with the _latchODE.m_ function to create a vector of time and latch positions and velocities. Using the functions described above, the latch acceleration; load position, velocities, and accelerations; normal, spring, and net forces; arc length travelled; and work due to friction are all calculated at each point in time above. From this, the position and velocity of the load at the unlatch point is determined and basic harmonic motion equations are used to determine the motion of the load and the point in time when the spring force reaches zero, which is defined as the take off time. After this point, the load has no external forces acting on it and simply moves at constant velocity for an amount of time determined by the overshoot field in the detailsStr structure. 

_findLoad.m_
When using a nonlinear model of friction, the maximum spring force that can be stored varies with the friction coefficient. We created this function to perform this calculation. During loading, the maximum force that can be stored is dependent upon the friction force and the spring force, which results in some normal force between the latch and the load. However, due to the nonlinear relationship between friction and normal force in the model chosen, this solution cannot be found analytically. Therefore, a nonlinear solver is implemented to solve it. 

_fricVar.m_
As mentioned in the Usage section, this function calls the _simulation.m_ function for each provided friction coefficient to assist in analysis of the effect of friction on this system. The inputs and outputs are described in detail in the Usage section. For each given friction coefficient, the loaded force is first determined. Based on whether the specified loading regime is "nonlinear" or "phenomenological", the loaded force is calculated either by calling the _findLoad.m_ function or by the phenomonelogical model of loading, defined in the Supplementary Information. Once the loaded force is calculated, the _simulation.m_ function is called for this friction coefficient and the results are stored in the output kinematicsa variable. This process is repeated for each friction coefficient. Once all values are analyzed, a new folder is created to store the results of the analysis. The four figures described in the Usage section are created and saved in this folder. 


