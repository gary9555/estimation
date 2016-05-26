function trackErrorNorm=run(simConst,estConst,doplot)
% run()
%
% Main function for Extended Kalman Filter programming exercise.
% It is used to simulate the true model, call the estimator and show the 
% results.
%
% Inputs:
%   simConst    constants used for the simulation (as in SimulationConstants.m)
%   estConst    estimator constants (as in EstimatorConstants.m), constants
%               used by the estimator
%   doplot      if doplot=true      plots are generated
%               if doplot=false     no plots are generated
%
% Output:
%   trackErrorNorm  the tracking error e, as defined in the exercise sheet
%
% You can run the function simply with "run()" if you want to use the
% default parameters as provided by EstimatorConstants.m and
% SimulationConstants.m, and generate the plots.
%
% Class:
% Recursive Estimation
% Spring 2016
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Michael Muehlebach
% michaemu@ethz.ch
%
% --
% Revision history
% [15.04.11, ST]    first version by Sebastian Trimpe
% [30.04.12, PR]    2012 version, added unknown wheel radius
% [06.05.13, MH]    2013 version
% [24.04.15, MM]    2015 version
% [14.04.16, MM]    2016 version

% clear command window, close figures
clc;
close all;

if nargin==0
   % Define the simulation constants that are used in simulation, but not 
   % accessible to the estimator.  
   simConst = SimulationConstants();
   
   % Define the constants accessible to the estimator.
   estConst = EstimatorConstants();
   
   % Generate plots by default.
   doplot=true;
end

%% Setup

% Set the random number generator state.
% Uncomment to make results reproducable. This setting was used to generate
% the plot in the problem description.
% rand('seed',1);
% randn('seed',1);


%% Simulation
% The function 'Simulator' simulates the robot kinematics and generates
% measurements.
[tm, wheelRadius, loc, input, sense ] = Simulator( simConst );

%% Run the Estimator

N=simConst.N;

% Initialize the estimator.  
estState = [];
posEst = zeros(N,2);
oriEst = zeros(N,1);
radiusEst = zeros(N,1);
posVar = zeros(N,2);
oriVar = zeros(N,1);
radiusVar = zeros(N,1);

[posEst(1,:),oriEst(1),radiusEst(1),posVar(1,:),oriVar(1),radiusVar(1),estState] = ...
    Estimator(estState,zeros(1,2),zeros(1,2),0,estConst);

% Call the estimator for each time step.
for n = 2:N
    [posEst(n,:),oriEst(n),radiusEst(n),posVar(n,:),oriVar(n),radiusVar(n),estState] = ...
        Estimator(estState,input(n,:),sense(n,:),tm(n),estConst);
end



%% The results
% Plots of the results.

% Calculate the total tracking error.
% Replace the following:
errorVecX = loc(:,1) - posEst(:,1);
errorVecY = loc(:,2) - posEst(:,2);


trackErrorNorm = sqrt(sum(errorVecX.'*errorVecX + errorVecY.'*errorVecY) / N)

if doplot
    % Add your plots here to debug the estimator and verify your
    % implementation.
    %%% plot x %%%%
    figure(1);
    plot(loc(1:N,(1)),'-g');
    hold;
    plot(posEst(1:N,(1)),'-r');
      
    %%% plot r %%%
%     figure(2);
%     plot(loc(1:N,(3)),'-g');
%     hold;
%     plot(oriEst,'-r');
    
    %%% plot 2d position %%%
    figure(3);
    plot(loc(1:N,(1)),loc(1:N,(2)),'-g');
    hold;
    plot(posEst(1:N,(1)),posEst(1:N,(2)),'-r');

    %%% plot x variance %%%
      figure(4);
      plot(posVar(1:N,(1)),'-g');
end
    
return;
