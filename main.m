close all
clear all

% --- system setup ----------------------
rng(2020)
T = 15;                             % time domain = [0,T]
n = 2;  t = linspace(0,T,n+1);     % snapshots at which the QoI's are computed
L = 40;                             % space domain = [-L,L]
m = 300; x = linspace(-L,L,m+1);    % spatial grid
u0 = rand(2,m+1)*0.2+0.4;           % fixed initial condition: random number in [0.4, 0.6]

% --- reaction-diffusion system ---------
G = @(theta) RDsystem(theta,t,x,u0);
% EXPLANATION:
% qoi = G(theta)
%   Input: theta = [log(D1), log(D2), R10, R11, R12, R13, R20, R21, R22, R23]
%          theta contains 10-dim prefactors.
%          theta(1) & theta(2) are the log of diffusivities
%   Output: qoi
%       qoi(:,i) = QoI's for the i-th snapshot
%       qoi(1,:) = mean of the local minima of the field
%       qoi(2,:) = mean of the local maxima of the field
%       qoi(3,:) = mean of the bump-size distribution
%       qoi(4,:) = std of the bump-size distribution

% --- generate QoI ------------
%theta =[      D1,       D2, R10, R11, R12, R13, R20, R21, R22, R23]
%theta =[log(.01), log(0.4), 0.1,  -1,   0,   1, 0.9,   0,   0,  -1];
%theta =[log(.03), log(1.2), 0.1,  -1,   0,   1, 0.9,   0,   0,  -1];
%theta =[log( .1), log(  4), 0.1,  -1,   0,   1, 0.9,   0,   0,  -1];
% Can use the above three example theta's as the "true" theta.
theta = [log(.03),  log(1.2), 0.1,  -1,   0,   1, 0.9,   0,   0,  -1];
qoi = G(theta);
nparams = length(theta);

% --- MCMC ------------

% select the method used
method = 'dram'; % 'mh','am','dr', or 'dram', see below

% data
clear model data params options
data.G = G;
data.ydata = qoi;

rel_QoI_sigmas = [0.2 0.2 0.2 0.2]'; % Relative noise std, ie 0.1 = 10% of signal.
% Computes the individual sigma values based on the relative noise std.
for i=1:length(rel_QoI_sigmas)
  full_QoI_sigmas(i,:) = rel_QoI_sigmas(i) * abs(qoi(i,:));
end
data.full_QoI_sigmas = full_QoI_sigmas;

% parameters for mcm
switch method
 case 'mh'
   nsimu    = 3000;
   drscale  = 0;
   adaptint = 0;
 case 'dr'
  nsimu    = 3000;
  drscale  = 2; 
  adaptint = 0;
 case 'am'
  nsimu    = 3000;
  drscale  = 0; 
  adaptint = 100;
 case 'dram'
  nsimu    = 30000;
  drscale  = 2; 
  adaptint = 100;
end

% create input arguments for the dramrun function

model.lnLikelihood    = @my_lnLikelihood;
model.lnPrior = @my_lnPrior;

%params.par0    = ones(1,nparams)*0.5; % initial parameter values
params.par0 = theta;
params.n       = -1;    % number of observations % DEACTIVATE FOR NOW
params.sigma2  = 1;  % prior for error variance sigma^2 % DEACTIVATE FOR NOW
params.n0      = 1;    % prior accuracy for sigma^2

params.parmu0   = zeros(1,nparams);                % prior mean of theta
%params.parmu0 = [log(.01),  log(.4), 0.1,  -1,   0,   1, 0.9,   0,   0,  -1];
params.parsig0  = ones(1,nparams)*2;            % prior std of theta

options.nsimu    = nsimu;               % size of the chain
options.adaptint = adaptint;            % adaptation interval
options.drscale  = drscale;
%options.qcov     = eye(nparams,nparams)*1e-2;      % initial proposal covariance 
options.qcov = [
    2.0424e-03   8.8944e-04  -1.6086e-03   5.6531e-04   8.2412e-04  -1.9139e-03   7.1007e-04  -1.0130e-03   3.7535e-04   1.3850e-03
   8.8944e-04   1.2063e-03  -9.7559e-04   5.1325e-04   3.8824e-04  -2.4652e-03  -1.1815e-03   9.3982e-04   2.2758e-03   5.8860e-04
  -1.6086e-03  -9.7559e-04   2.1601e-03  -7.5744e-04  -9.7002e-04   1.4490e-03   1.2561e-04   7.5514e-04  -1.0151e-03  -1.4357e-03
   5.6531e-04   5.1325e-04  -7.5744e-04   3.9994e-04   3.3254e-04  -1.0868e-03  -4.6759e-04   2.4861e-04   8.8563e-04   4.3719e-04
   8.2412e-04   3.8824e-04  -9.7002e-04   3.3254e-04   4.8037e-04  -5.6714e-04   1.4810e-04  -5.1613e-04   1.7573e-04   6.6442e-04
  -1.9139e-03  -2.4652e-03   1.4490e-03  -1.0868e-03  -5.6714e-04   6.4959e-03   2.6233e-03  -2.7536e-03  -5.2813e-03  -9.2467e-04
   7.1007e-04  -1.1815e-03   1.2561e-04  -4.6759e-04   1.4810e-04   2.6233e-03   3.5430e-03  -3.1639e-03  -4.2236e-03   4.0791e-04
  -1.0130e-03   9.3982e-04   7.5514e-04   2.4861e-04  -5.1613e-04  -2.7536e-03  -3.1639e-03   3.6649e-03   3.8109e-03  -9.7155e-04
   3.7535e-04   2.2758e-03  -1.0151e-03   8.8563e-04   1.7573e-04  -5.2813e-03  -4.2236e-03   3.8109e-03   6.4720e-03   3.8776e-04
   1.3850e-03   5.8860e-04  -1.4357e-03   4.3719e-04   6.6442e-04  -9.2467e-04   4.0791e-04  -9.7155e-04   3.8776e-04   1.2908e-03
   ];

options.printint = 100; % display frequency

% run the chain
tic
[results,chain,~,chainll] = dramrun(model,data,params,options);
toc

save w3.mat
%tau=iact(chain); % Integrated Autocorrelation Time
%%
figure(2);clf
plot(chain(:,1),chain(:,2),'.')
xlabel('k_1'); ylabel('k_2');title('MCMC chain');
% add 95% ellipses of the proposal to the plot
hold on;axis manual
%ellipse(kopt+[100 100],cmat*6,'Linewidth',1,'Color','black')
ellipse(mean(chain(:,1:2)),results.R(:,1:2)'*results.R(:,1:2)*6,'Linewidth',1,'Color','red')
hold off; axis normal

for i=1:nparams
  figure(100+i)
  plot(chain(:,i),'-')
  xlabel('Iter')
  ylabel(['\theta_{',num2str(i),'}'])
end
%%
% tune for speed
% nsimu = 300 steps; adaptint = 10 steps
% 1) mesh size = 400, no halt, time = 304.56s ~ 1.015s/step
% 2) mesh size = 400, halt solver if soln > 100, time = 321.43s ~ 1.07s/step
% 3) mesh size = 300, no halt, time = 206.03s ~ 0.687s/step, iact = [33.3962, 28.9054, 26.9823, 12.2057, 26.5132, 27.9502, 17.9670, 29.9875, 31.1914, 32.0952]
% 4) mesh size = 300, no halt, no DR, 134.17s ~ 0.447s/step, iact = [31.2557, 33.4225, 13.6346, 36.5620, 29.8325, 31.5109, 29.5755, 29.8422, 20.2711, 28.5394]
% nsimu = 30000; adaptint = 100 steps
% * mesh size = 300, halt if soln > 100, time = 22222.45s ~ 0741s/step
