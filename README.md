function [posEst,oriEst,radiusEst, posVar,oriVar,radiusVar,estState] = Estimator(estState,actuate,sense,tm,estConst,n)
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
%     estState.var = [ 0.4    0   0     0;       % orientation variance 
%                       0   0.33  0     0;      % x variance
%                       0     0  0.33   0;      % y variance
%                       0     0   0     0.1 ] ;    % radius variance
    estState.var = [ 0.05    0   0     0;       % orientation variance 
                      0   0.4  0     0;      % x variance
                      0     0  0.4   0;      % y variance
                      0     0   0     0.0005 ] ;    % radius variance
    
    % Output
    posEst = [estState.est(2) estState.est(3)];
    oriEst = estState.est(1);
    posVar = [estState.var(2,2) estState.var(3,3)];
    oriVar = estState.var(1,1);
    radiusEst = estConst.NominalWheelRadius;
    radiusVar = estState.var(4,4);
    return;
end


%% Mode 2: Estimator iteration.
% If we get this far tm is not equal to zero, and we are no longer
% initializing.  Run the estimator.

%%%%%% prior update %%%%%%%%
% be aware that gamma is bounded
delta_t = tm - estState.prev_t;
r = estState.est(1);
x = estState.est(2);
y = estState.est(3);
W = estState.est(4);
u_v = actuate(1);
u_r = actuate(2);
Q_v = estConst.VelocityInputPSD;
Q_r = estConst.AngleInputPSD;
B = estConst.WheelBase;

% integrate x'=q(x, u, 0, t), time step is T
estState.est = estState.est + delta_t * [ -W * u_v * sin(u_r) / B;
                                W * u_v * cos(u_r) * cos(r);
                                W * u_v * cos(u_r) * sin(r);
                                0 ];
                            
estState.est(1) = mod(estState.est(1), 2*pi);

% Let A_t = dq / dx ,  L_t = dq / dv   , P_t the variance matrix of the states                
A_t = [ 0                               0    0      -u_v * sin(u_r) / B;
        -W * u_v * cos(u_r) * sin(r)    0    0      u_v * cos(u_r) * cos(r);
        W * u_v * cos(u_r) * cos(r)     0    0      u_v * cos(u_r) * sin(r);
        0                               0    0              0               ];

L_t = [ -W * u_v * sin(u_r) / B         -W * u_v * cos(u_r) / B;
        W * u_v * cos(u_r) * cos(r)     -W * u_v * sin(u_r) * cos(r);
        W * u_v * cos(u_r) * sin(r)     -W * u_v * sin(u_r) * sin(r);
                   0                                0               ]; 

P_p = estState.var;
P_p = P_p  +  delta_t * (A_t * P_p  +  P_p * transpose(A_t)  +  L_t * [Q_v 0; 0 Q_r] * transpose(L_t));
              
               
%%%%%%% measurement update %%%%%%%%
r_p = estState.est(1);
x_p = estState.est(2);
y_p = estState.est(3);
W_p = estState.est(4);
temp = sqrt(x_p^2 + y_p^2);
R_r = estConst.CompassNoise;
R_d = estConst.DistNoise;

% Let H_k = dh / dx , M_k = dh / dw
H_k = [ 1           0           0           0;
        0        x_p/temp    y_p/temp       0 ];

M_k = [ 1 0; 0 1];

K_k = P_p * transpose(H_k) / (H_k * P_p * transpose(H_k) + M_k * [R_r 0; 0 R_d] * transpose(M_k));

if sense(2)==Inf
    sense(2) = r_p;
end
if sense(1) ==Inf
    sense(1)=temp;
end
estState.est = estState.est + K_k * ([sense(2) - r_p; sense(1) - temp]);

estState.var = (eye(4) - K_k*H_k)*P_p;
% make sure the W_p doesnt't go out of bounds
estState.est(4) = min(estState.est(4),0.15);
estState.est(4) = max(estState.est(4),0.05);

% change previous timestamp
estState.prev_t = tm;

% Replace the following:
posEst = [estState.est(2) estState.est(3)];
oriEst = estState.est(1);
posVar = [estState.var(2,2) estState.var(3,3)];
oriVar = estState.var(1,1);
radiusEst = estState.est(4);
radiusVar = estState.var(4,4);

end
