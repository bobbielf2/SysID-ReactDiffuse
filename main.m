close all
clear all

% --- system setup ----------------------
rng(2020)
T = 15;                             % time domain = [0,T]
n = 2;  t = linspace(0,T,n+1);     % snapshots at which the QoI's are computed
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
params.par0 = [.01,  .4, 0.1,  -1,   0,   1, 0.9,   0,   0,  -1];
params.n       = -1;    % number of observations % DEACTIVATE FOR NOW
params.sigma2  = 1;  % prior for error variance sigma^2 % DEACTIVATE FOR NOW
params.n0      = 1;    % prior accuracy for sigma^2

params.parmu0   = zeros(1,nparams);                % prior mean of theta
%params.parmu0 = [.01,  .4, 0.1,  -1,   0,   1, 0.9,   0,   0,  -1];
params.parsig0  = ones(1,nparams)*2;;            % prior std of theta

options.nsimu    = nsimu;               % size of the chain
options.adaptint = adaptint;            % adaptation interval
options.drscale  = drscale;
%options.qcov     = eye(nparams,nparams)*1e-2;      % initial proposal covariance 
options.qcov = [    9.6279e-06    1.283e-05   0.00011859   -0.0004412   -2.166e-05 ...
                    0.00040131  -0.00011012   0.00036685  -3.8164e-05   -0.00025168;
    1.283e-05    0.0020935  -0.00016269    0.0023743  -0.00012774   -0.0058966    0.0028876 ...
                    -0.0028077   0.00022461     0.0030311;
   0.00011859  -0.00016269     0.010236    -0.027527   -0.0025137     0.021304   -0.0062778 ...
                    0.021252   -0.0027936     -0.014034;
   -0.0004412    0.0023743    -0.027527     0.090705    0.0056354    -0.079086 ...
                    0.023473    -0.073306    0.0085081      0.051631;
   -2.166e-05  -0.00012774   -0.0025137    0.0056354   0.00074959   -0.0036176   0.00098428 ...
                    -0.0042627    0.0004458     0.0026801;
   0.00040131   -0.0058966     0.021304    -0.079086   -0.0036176     0.085881    -0.029933 ...
                    0.065416   -0.0073013     -0.047901;
  -0.00011012    0.0028876   -0.0062778     0.023473   0.00098428    -0.029933 ...
                    0.012035    -0.019366    0.0021255      0.014583;
   0.00036685   -0.0028077     0.021252    -0.073306   -0.0042627     0.065416    -0.019366 ...
                    0.061729   -0.0062406     -0.045181;
  -3.8164e-05   0.00022461   -0.0027936    0.0085081    0.0004458   -0.0073013    0.0021255 ...
                    -0.0062406    0.0018947     0.0039184;
  -0.00025168    0.0030311    -0.014034     0.051631    0.0026801    -0.047901 ...
                    0.014583    -0.045181    0.0039184       0.03518];

options.printint = 100; % display frequency

% run the chain
[results,chain,~,chainll] = dramrun(model,data,params,options);

%tau=iact(chain); % Integrated Autocorrelation Time

% figure(2);clf
% plot(chain(:,1),chain(:,2),'.')
% xlabel('k_1'); ylabel('k_2');title('MCMC chain');
% % add 95% ellipses of the proposal to the plot
% %hold on;axis manual
% %ellipse(kopt+[100 100],cmat*6,'Linewidth',1,'Color','black')
% %ellipse(mean(chain)-[50,0],results.R'*results.R*6,'Linewidth',1,'Color','red')
% %hold off; axis normal

for i=1:nparams
  figure(100+i)
  plot(chain(:,i),'-')
  xlabel('Iter')
  ylabel('\theta_1')
end