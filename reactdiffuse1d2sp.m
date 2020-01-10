function [sol, t, x, u0, sole, te, ie] = reactdiffuse1d2sp(t,x,u0,th_val,th_ind)
% Reaction-Diffusion equation 1-D, 2 species.
% u_t = D1*u_xx + R1(u,v)
% v_t = D2*v_xx + R2(u,v)
%   where D1 & D2 are diffusivities, R1 & R2 are the reaction terms
% This system is solved by the built-in "pdepe" solver

% th=[ D1, D2, R10, R11, R12, R13, R20, R21, R22, R23]
th = [.01, .4, 0.1,  -1,   0,   1, 0.9,   0,   0,  -1]; % default param
if nargin >= 4 && ~isempty(th_val)
    if nargin < 5 || isempty(th_ind)
        th_ind = 1:numel(th_val);
    elseif numel(th_val)~=numel(th_ind)
        error('th_val and th_ind need to have same length!')
    end
    th(th_ind) = th_val; % modify the default param
end

% define space and time domains
if nargin < 1 || isempty(t)
    T = 15; n = 20;
    t = linspace(0,T,n+1);      % time = [0,T]
else
    T = max(t); n = numel(t)-1;
end
if nargin < 2 || isempty(x)
    L = 40; m = 400;
    x = linspace(-L,L,m+1);     % domain = [-L,L]
else
    L = max(x); m = numel(x)-1;
end


% define PDE
mypde = @(x,t,u,DuDx) RDpde(x,t,u,DuDx,th);

% define initial condition
if nargin < 3 || isempty(u0)
    u0 = rand(2,m+1)*0.2+0.4; % random number in [0.4, 0.6]
end
myic = @(xi) RDic(xi,x,u0);

% define boundary condition
mybc = @RDbc;

% define event function
opt = odeset('Events',@myEvent);

%tic
[sol,t,sole,te,ie] = pdepe(0,mypde,myic,mybc,x,t,opt);
% fprintf('ie = %d\n ',ie)
%sol = pdepe(0,mypde,myic,mybc,x,t);
%sole = []; te = []; ie = [];
%toc

% ie indicates which species's magnitude has exceeded 100 
% te is the time when the above event happens
% sole are the solutions at te
% The program terminates at te, to estimate the magnitude of each solution
% at time t_final, we assume exponential growth/decay overtime, then
%   max(u1(t_final)) ~ max(u1(te))^(t_final/te)
%   max(u2(t_final)) ~ max(u2(te))^(t_final/te)

% plot result
if 0 % change to 1 to plot the solution
    dt = T/n;
    for i = 1:numel(t)
        plot(x, sol(i,:,1)); hold on
        plot(x, sol(i,:,2)); hold off
        legend({'u','v'})
        ylim([0,4])
        title(['time = ',num2str(dt*(i-1))])
        drawnow
    end
   keyboard
end


% --------------------------------------------------------------
function [c,f,s] = RDpde(x,t,u,DuDx,th)
c = [1; 1]; 
f = [th(1); th(2)] .* DuDx; 
s = [th(3) + th(4) * u(1) + th(5) * u(2) +  th(6) * u(1).^2 .* u(2);
     th(7) + th(8) * u(1) + th(9) * u(2) + th(10) * u(1).^2 .* u(2)];
% --------------------------------------------------------------
function u0 = RDic(x,X,UV0)
u0 = UV0(:,X == x);
% --------------------------------------------------------------
function [pl,ql,pr,qr] = RDbc(xl,ul,xr,ur,t)
pl = [0; 0]; 
ql = [1; 1]; 
pr = [0; 0]; 
qr = [1; 1]; 
% --------------------------------------------------------------
function [value,isterminal,direction] = myEvent(m,t,xmesh,umesh)
% Don't let magnitude of solutions go beyond 100.
n = numel(xmesh);
k = numel(umesh)/n;
value = zeros(k,1);
for i = 1:k
    value(i) = max(abs(umesh(i:k:end)),[],'all')>1e2;
end
isterminal = true(size(value));
direction = zeros(size(value));
