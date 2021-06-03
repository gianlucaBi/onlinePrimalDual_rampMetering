% This script simulates three control policies for a two-link network
% modeled using the CTM. 
%
% The control policies can be selected by settign the 'ctr_type' variable
% 1=no control, 
% 2=ALINEA, 
% 3=primal-dual
%
% Author: Gianluca Bianchin
% Department of Electrical, Computer, and Energy Engineering 
% University of Colorado Boulder 
% gianluca.bianchin@colorado.edu
%
% Date: June 3, 2021
%


clear all; close all; clc;
global n m r G B phi beta dMax sMax xJam out_neighb R in_neighb uRef ctr_type xCrt Kp xRef Qu Qx Ky e nu alpha wt ALINEA_ctr gain


%% Select controller type:
ctr_type    = 3;        % 1=no control, 2=ALINEA, 3=primal-dual

%% Initializations        
tMax    = 100;
n       = 64;                   % Number of highways
m       = 17;                   % Number of control inputs
r       = n;
off_ramps   = 47:64;            % list of off-ramps
load networkModel               % loads the matrices in_neighb, in_neighb, and R
B       = eye(n,m);
gain    = 1;                    % When larger than 1, makes traffic system faster 
tAux    = 0:.1:tMax;

uRef    = 5*ones(m,1);                          % Traffic demand on onramps
% Define exogenous input w_t
wDmean  = .5*square(.1*ones(r,1)*tAux);         % Mean value of the signal
% wDmean   = .1*sin(.1*ones(r,1)*tAux);
wD          = 0*(wDmean + .1*randn(r,length(tAux)));    % Pick a realization of the signal w_t
wt          = @(t) interp1(tAux, wD', t)';
figure; 
plot(tAux, wDmean, 'linewidth',2)
legend('mean(w)')

% Demand and supply parameters
phi     = 1*ones(n,1);                      % Free-flow velocity
beta    = 1*ones(n,1);                      % Congested velocity
dMax    = 3*ones(n,1);                      % Demand saturation value
sMax    = 3*ones(n,1);                      % Supply saturation value
xJam    = 6*ones(n,1);                      % Jam density
xCrt_d  = dMax./phi;                        % Critical density demand 
xCrt_s  = xJam-sMax./beta;                  % Critical density supply 
xCrt    = min(xCrt_d,xCrt_s);

% Plot demand-supply diagram 
x       = [0:.1:xJam].*ones(n,1);
for ii=1:size(x,2)
    s(:,ii)   = supply(x(:,ii));
    d(:,ii)   = demand(x(:,ii));
end
figure; plot(x(1,:),d(1,:),'linewidth',2)
hold on; plot(x(1,:),s(1,:),'linewidth',2)
legend('flow-density diagram')


% Parameters of the CTM in the free-flow regime
R_T         = R';
A           = (R_T-eye(n))*diag(phi);
C           = eye(n);
G           = -C*inv(A)*B;


%% Controllers parameters
% ALINEA
ALINEA_ctr  = [25; 26; 27; 29; 30; 31; 34; 35; 36; 37; 40; 41; 42; 43; 44; 45; 46];     % for each input this vector contains the downstream link to be controlled by ALINEA
Kp          = 50*ones(m,1);

% Optimization problem parameters
nu      = .5;                   % Regularization parameter
Ky      = eye(n);
e       = xCrt;
Qx      = 1*eye(n);             % Output cost matrix
Qu      = 1*eye(m);             % Input cost matrix
xRef    = xCrt;
alpha   = 10;





%% Solve for states using ODE solver     
x0          = 0*rand(n,1);          % Initial condition state
u0          = 0*rand(m,1);          % Initial condition control input
lm0         = 0*rand(r,1);          % Initial condition Lagrange multipliers
[tODE,zODE] = ode45(@(t,x) diffEq(t,x),[0 tMax],[x0; u0; lm0]);                 % Solve ODE 


% Extract interesting quantities for plot purposes
for jj=1:length(tODE)
    x                       = zODE(jj,1:n)';
    [fIn(:,jj), fOut(:,jj)] = getflows(x);
    Phi(jj)                 = sum(fOut(off_ramps,jj));
    cstrViol(jj)            = norm(max(Ky*x-e,0));
    if jj>1
        avgPhi(jj)              = trapz(tODE(1:jj),Phi)./tODE(jj);
    end
end

% States at final time
xODE    = zODE(:,1:n);
uODE    = zODE(:,n+1:n+m);
lmODE   = zODE(:,n+m+1:n+m+r);
x       = zODE(end,1:n)';
u       = zODE(end,n+1:n+m)';
lm      = zODE(end,n+m+1:n+m+r)';





%% Outputs and figures
figure; plot(tODE,cstrViol,'linewidth',2)
legend('||Kx-e|| (Constraint violation) ')


plt_list    = [39 43 46];        % list of links you would like to visualize
figure;
subplot(4,1,1)
plot(tODE,xODE(:,plt_list),'linewidth',2)
legend('x(t) (Densities)')
subplot(4,1,2)
plot(tODE, uRef.*ones(1,length(tODE)),'-.','linewidth',2); hold on
plot(tODE, uODE,'linewidth',2)
legend({'Input uRef','Input u'})
subplot(4,1,3)
plot(tODE,fOut(plt_list,:),'linewidth',2)
legend({'Link out-flows'})
subplot(4,1,4)
plot(tODE,Phi,'linewidth',2)
hold on

% Compare with MPC 
% Note:     'workspace_MPC.m' has been generated from a previous run of 
%           the same code with 'ctr_type = 3'
load workspace_MPC
plot(tMem,Phi,'linewidth',2)
legend({'GF','MPC'})










%% Functions
function [dz,y] = diffEq(t,z)
    global n m r G B uRef ctr_type xCrt Kp xRef Qu Qx Ky e nu alpha wt ALINEA_ctr gain
     
    x       = z(1:n);
    u       = z(n+1:n+m);
    lm      = z(n+m+1:n+m+r);
    y       = x + wt(t);
    
        %% Controller dynamics
    switch ctr_type
        case 1      % No control
            u       = uRef;
            du      = zeros(m,1);
            dlm     = zeros(r,1);
        case 2      % ALINEA
            du      = Kp.*(xCrt(ALINEA_ctr)-y(ALINEA_ctr));         %[xCrt(3)-x(3); xCrt(4)-x(4)];
            u       = max(min(u, uRef),0);         % Inflow cannot be larger than demand
            dlm     = zeros(r,1);
        case 3      % Saddle-flow
            Lu      = Qu*(u-uRef) + G'*Qx*(y-xRef) + G'*Ky'*lm; 
            Pu      = max(u-alpha*Lu,0);
            du      = Pu-u;                  
    
            Llm     = Ky*y-e-nu*lm;
            Plm     = max(lm+alpha*Llm,0);
            dlm     = Plm-lm;               
    end
    
    
    
    %% Plant dynamics
    [fIn, fOut] = getflows(x);
    dx          = gain*(-fOut + fIn + B*min(u,uRef));
    
    
    
    %% Updates
    dz      = [dx; du; dlm];
end

% CTM demand funtion
function s = demand(x)
    global phi dMax 
    s = min(phi.*x,dMax);
end

% CTM supply funtion
function s = supply(x)
    global beta sMax xJam
    s = max(min(beta.*(xJam-x),sMax),0);
end

% Compute outflow from a link, uses FIFO allocation policy
function [fIn, fOut] = getflows(x)
    global n out_neighb R in_neighb 
    d       = demand(x);
    s       = supply(x);
    % Compute outflows
    for ii=1:n
        sr = [];
        for jj=1:n
            if out_neighb(ii,jj)==1     % if j is an out-neighbor of i
                sr(jj)      = s(jj)/R(ii,jj);
            else 
                sr(jj)      = inf;
            end
        end
    fOut(ii)    = min([d(ii), sr]);
    end
    % Compute inflows
    for ii=1:n
        ff = 0;
        for jj=1:n          % if j is an out-neighbor of i
            if in_neighb(ii,jj)==1   
                ff = ff + R(jj,ii)*fOut(jj);
            end
        end
        fIn(ii)    = ff;
    end
    fIn         = fIn';
    fOut        = fOut';
end

