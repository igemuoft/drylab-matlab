iGEM_Toronto_MATLAB_model.m - uses square(). Works with the symbols toolbox (original by Andrew West)

iGEM_Toronto_MATLAB_model_edited.m - doesn't use square() function. Works without the additional toolbox. Also made easier to customize the total simulation time and light switch toggle time.

iGEM_Toronto_MATLAB_model_edited2.m - displays the equilibrium values for x2, theta, lambda, under both light on and off conditons.

iGEM_Toronto_MATLAB_model_edited4.m - displays graph for for 2D graph, parasweep, eqn_time

parasweep_kG.m - sweepts the ODE for theta over parameters (k,gamma_theta) for fixed amount of time

eqm_time.m - Determines equilibrium time from light off to light on, defined as the time required for lambda to increase by 85% of difference between light_on equilibrium and light_off equilibrium

rk.m - original runge-kutta function from the iode package

rk2.m - modified runge-kutta for use in main file

