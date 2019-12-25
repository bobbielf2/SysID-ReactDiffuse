% --- system setup ----------------------
T = 15;                             % time domain = [0,T]
n = 10;  t = linspace(0,T,n+1);     % snapshots at which the QoI's are computed
L = 40;                             % space domain = [-L,L]
m = 400; x = linspace(-L,L,m+1);    % spatial grid
u0 = rand(2,m+1)*0.2+0.4;           % fixed initial condition: random number in [0.4, 0.6]

% --- reaction-diffusion system ---------
G = @(theta) RDsystem(theta,t,x,u0);
% EXPLANATION:
% qoi = G(theta)
%   Input: theta = [D1, D2, R10, R11, R12, R13, R20, R21, R22, R23]
%          theta contains 10-dim prefactors.
%          REQUIRE: theta(1) > 0 and theta(2) > 0
%   Output: qoi
%       qoi(:,i) = QoI's for the i-th snapshot
%       qoi(1,:) = mean of the local minima of the field
%       qoi(2,:) = mean of the local maxima of the field
%       qoi(3,:) = mean of the bump-size distribution
%       qoi(4,:) = std of the bump-size distribution

% --- generate QoI ------------
%theta =[ D1,  D2, R10, R11, R12, R13, R20, R21, R22, R23]
%theta =[.01, 0.4, 0.1,  -1,   0,   1, 0.9,   0,   0,  -1];
%theta =[.03, 1.2, 0.1,  -1,   0,   1, 0.9,   0,   0,  -1];
%theta =[ .1,   4, 0.1,  -1,   0,   1, 0.9,   0,   0,  -1];
% Can use the above three example theta's as the "true" theta.
theta = [.01,  .4, 0.1,  -1,   0,   1, 0.9,   0,   0,  -1];
qoi = G(theta);

