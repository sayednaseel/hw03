% run consecutive steps for rimless wheel

clear all;close all;clc

% Set a tolerance for testing
abstol = 1e-4;

% Use these parameter values
parms = struct('alpha', 0.3, 'rgyr', 0, 'gamma', 0.04, 'tmax', 5);
thetadotplusstar = -sqrt(4*sin(parms.alpha)*parms.gamma/sin(2*parms.alpha)^2)*cos(2*parms.alpha);% angular velocity at limit cycle by analytical derivation
x0 = [0.3; -0.5]; % a sample initial state
%x0 = [0.3; -0.2]; % a sample initial state

N = 20;
x = [];  % store a big vector of states for multiple steps
t = [];  % big time vector
e = [];  % big energies vector
x0s = x0'; % store the initial states
xprev = x0;
tprev = 0; % cumulative time of previous step end

for i = 1:N
    [xnext, te, xs, ts, energies] = rimlesswheelstep(xprev,parms);
    x = [x; xs]; % states are stored as a column of thetas, col of thetadots
    x0s = [x0s; xnext']; % only the step initial conditions
    t = [t; ts+tprev];
    e = [e; energies];
    tprev = t(end); % update cumulative time
    xprev = xnext;
end
xends=x0s;


figure('name','states v time')
plot(t,x)
xlabel('time (dimensionless)');
ylabel('states (dimensionless)');
title('states vs. time for rimless wheel');
legend('theta', 'thetadot');
grid on;
%% Energy plot
%clf
figure('name','Energies')
subplot(121)
plot(t,e, '.')
xlabel('time (dimensionless)');
ylabel('energy (dimensionless)');
title('energies vs. time for rimless wheel');

subplot(122)
plot(ts,energies) % zoom in on energy for the last step
xlabel('time (dimensionless)');
ylabel('energy (dimensionless)');
title('energies vs. time for final step of rimless wheel');

%% Thetadot per step 
%clf
figure('name','Angular velocity ICs')
plot(0:length(xends)-1, xends(:,2),'.')
xlabel('step number');
ylabel('initial thetadot (dimensionless)');
title('initial thetadot vs. step number');
hold on
plot([0,length(xends)-1], thetadotplusstar*[1 1],':')
legend('thetadot', 'analytical limit cycle')

%% Test energy conservation of final step
% We expect energy to be conserved for rimless wheel
assert(std(energies) < abstol, 'Rimless wheel simulation does not conserve energy accurately');

%% State plot; not really a test, but demonstration plot
%clf
figure('name','states v time')
plot(t,x)
xlabel('time (dimensionless)');
ylabel('states (dimensionless)');
title('states vs. time for rimless wheel');
legend('theta', 'thetadot');
grid on;
%% Test computational fixed point solution

parms = struct('alpha', 0.3, 'rgyr', 0, 'gamma', 0.04, 'tmax', 5);
x0 = [0.3; -0.5]; % a sample initial state

% make a one-argument function for root-finding, including custom parms
rwerror = @(x) rimlesswheelerror(x, parms);

xstar = findroot(rwerror, x0); % solve for limit cycle

% get data from a limit cycle
[xnext, te, xs, ts, energies] = rimlesswheelstep(xstar);

%% Test against analytical solution
% whether fixed point matches
% analytical solution

% Analytical solution predicts thetaplus
thetadotplusstar = -sqrt(4*sin(parms.alpha)*parms.gamma/sin(2*parms.alpha)^2)*cos(2*parms.alpha);
thetadotfinal = xstar(2); % thetadotplus for limit cycle 
% simulated steps

assert(abs(thetadotplusstar - thetadotfinal) < abstol,'error');

%% Parameter study of speed vs slope
gammas = linspace(0.03, 0.06, 15);
thetadot = -0.3; % first try for thetadot

xstars = zeros(length(gammas),2);
for i = 1:length(gammas)
    parms = struct('alpha', 0.3, 'rgyr', 0, 'gamma', gammas(i), 'tmax', 10);
    x0 = [0.3; thetadot]; % a sample initial state
    rwerror = @(x) rimlesswheelerror(x, parms);
    xstar = findroot(rwerror, x0); % solve for limit cycle
    xstars(i,:) = xstar';
    thetadot = xstar(2);
    [xnext,te,xs,ts]=rimlesswheelstep(xstar,parms);
    speeds(i) = 2*sin(parms.alpha)/te;
end

%clf;
figure('name','speed vs slope')
plot(gammas, speeds);
xlabel('Slope gamma (rad)');
ylabel('Speed (normalized)');
title('Rimless wheel speed vs. slope');
axis([0 Inf 0 Inf]) % helpful to include 0 in axes
grid on;
%% Parameter study of speed vs alpha

alphas = linspace(0.01, 0.33, 15);
thetadot = -0.3; % first try for thetadot

xstars = zeros(length(alphas),2);
for i = 1:length(alphas)
    parms = struct('alpha', alphas(i), 'rgyr', 0, 'gamma', 0.04, 'tmax', 10);
    x0 = [alphas(i); thetadot]; % a sample initial state
    rwerror = @(x) rimlesswheelerror(x, parms);
    xstar = findroot(rwerror, x0); % solve for limit cycle
    xstars(i,:) = xstar';
    thetadot = xstar(2);
    [xnext,te,xs,ts]=rimlesswheelstep(xstar,parms);
    speeds(i) = 2*sin(parms.alpha)/te;
end

%clf;
figure('name','speed vs alpha')
plot(alphas, speeds);
xlabel('Interleg angle alpha (rad)');
ylabel('Speed (normalized)');
title('Rimless wheel speed vs. Interlegangle');
%axis([0 Inf 0 Inf]) % helpful to include 0 in axes
grid on;
