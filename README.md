# estimation

function [posEst,oriEst,radiusEst, posVar,oriVar,radiusVar,estState] = Estimator(estState,actuate,sense,tm,estConst)
% [posEst,oriEst,posVar,oriVar,baseEst,baseVar,estState] =
% 	Estimator(estState,actuate,sense,tm,knownConst,designPart)
%
% The estimator.
%
% The function will be called in two different modes:
% If tm==0, the estimator is initialized; otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k), [1x2]-vector
%                   actuate(1): u_V, drive wheel angular velocity
%                   actuate(2): u_R, drive wheel angle
%   sense           sensor measurements z(k), [1x2]-vector, INF if no
%                   measurement
%                   sense(1): z_d, distance measurement
%                   sense(2): z_r, orientation measurement
%   tm              time, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConstants.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): x position estimate
%                   posEst(2): y position estimate
%   oriEst          orientation estimate (time step k), scalar
%   radiusEst       estimate of wheel radius W (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   radiusVar       variance of wheel radius estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
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
% [19.04.11, ST]    first version by Sebastian Trimpe
% [30.04.12, PR]    adapted version for spring 2012, added unknown wheel
%                   radius
% [06.05.13, MH]    2013 version
% [23.04.15, MM]    2015 version
% [14.04.16, MM]    2016 version


%% Mode 1: Initialization
if (tm == 0)
    % expectation vector
    estState.est = [0;                            % estimated orientation
                    0;                            % estimated x position 
                    0;                            % estimated y position
                    estConst.NominalWheelRadius]; % estimated wheel radius 
    
    % previous timestamp              
    estState.prev_t = 0; 
    
    % variance vector
    estState.var = [ 0.4;       % orientation variance 
                     0.33;      % x variance
                     0.33;      % y variance
                     0.01 ] ;    % radius variance
    % Output
    posEst = [estState.est(2) estState.est(3)];
    oriEst = estState.est(1);
    posVar = [estState.var(2) estState.var(3)];
    oriVar = estState.var(1);
    radiusEst = estConst.NominalWheelRadius;
    radiusVar = estState.var(4);
    return;
end


%% Mode 2: Estimator iteration.
% If we get this far tm is not equal to zero, and we are no longer
% initializing.  Run the estimator.

% prior update
% be aware that gamma is bounded
delta_t = tm - estState.prev_t;
r = estState.est(1);
x = estState.est(2);
y = estState.est(3);
W = estState.est(4);
u_v = actuate(1);
u_r = actuate(2);
B = estConst.WheelBase;

estState.est = estState.est + [ -delta_t * W * u_v * sin(u_r) / B;
                                delta_t * W * u_v * cos(u_r) * cos(r);
                                delta_t * W * u_v * cos(u_r) * sin(r);
                                0 ];
% Let A_t = dq / dx ,  L_t = dq / dv                    
A_t = [ 0                               0    0      -u_v * sin(u_r) / B;
        -W * u_v * cos(u_r) * sin(r)    0    0      u_v * cos(u_r) * cos(r);
        W * u_v * cos(u_r) * cos(r)     0    0      u_v * cos(u_r) * sin(r);
        0                               0    0              0               ];
% integrate x'=q(x, u, 0, t), time step is T
%Xp = estState.prev_Xm + []

% measurement update



% Replace the following:
posEst = [0 0];
oriEst = 0;
posVar = [0 0];
oriVar = 0;
radiusEst = 0;
radiusVar = 0;

end
