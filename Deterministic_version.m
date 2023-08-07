clc;
clear;

%This Program is written by Arman Behrad 

%% Parameters
tau = 10;
k = 0.75;
w = 3.6;
tau2 = 2;
sigma = 1;
T_end = 20*tau;
dt = 0.05*tau2;
t = 0:dt:T_end;

theta = @(x) x.^2./(k.^2 + x.^2); %Non-linearity term
%% simulate the time-evolution of the system

% Initialization
E1 = zeros(size(t));
E2 = zeros(size(t));
N1 = zeros(size(t));
N2 = zeros(size(t));
E1(1) = 0.2;
E2(1) = 0.6;
N1(1) = 0;
N2(1) = 0;

% Iterative integration
for i = 2:length(t)
    dN1 = -N1(i-1) + sigma*sqrt(tau2*dt)*randn();
    dN2 = -N2(i-1) + sigma*sqrt(tau2*dt)*randn();
    dE1 = (-E1(i-1) + theta(w*E1(i-1) - E2(i-1) - 0.5) + N1(i-1))*dt/tau;
    dE2 = (-E2(i-1) + theta(E1(i-1) - E2(i-1) - 0.5) + N2(i-1))*dt/tau;
    E1(i) = E1(i-1) + dE1;
    E2(i) = E2(i-1) + dE2;
    N1(i) = N1(i-1) + dN1;
    N2(i) = N2(i-1) + dN2;
end

% Plotting
figure;
plot(t, E1, 'r', 'LineWidth', 2);
hold on;
plot(t, E2, 'b', 'LineWidth', 2);
title('Time-Evolution of E1 and E2');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
legend('E1', 'E2');


%% Time-evolution in the State-State
% Initialization
E1 = zeros(length(t), 20);
E2 = zeros(length(t), 20);
N1 = zeros(length(t), 20);
N2 = zeros(length(t), 20);
E1(1, :) = 0.2;
E2(1, :) = 0.6;
N1(1, :) = 0;
N2(1, :) = 0;

% Iterative integration
for i = 2:length(t)
    dN1 = -N1(i-1, :) + sigma*sqrt(tau2*dt)*randn(1, 20);
    dN2 = -N2(i-1, :) + sigma*sqrt(tau2*dt)*randn(1, 20);
    dE1 = (-E1(i-1, :) + theta(w*E1(i-1, :) - E2(i-1, :) - 0.5) + N1(i-1, :))*dt/tau;
    dE2 = (-E2(i-1, :) + theta(E1(i-1, :) - E2(i-1, :) - 0.5) + N2(i-1, :))*dt/tau;
    E1(i, :) = E1(i-1, :) + dE1;
    E2(i, :) = E2(i-1, :) + dE2;
    N1(i, :) = N1(i-1, :) + dN1;
    N2(i, :) = N2(i-1, :) + dN2;
end

% Plotting
Colours = hsv(size(E1, 2));
figure;
plot(E1(:, 1), E2(:, 1), 'r', 'LineWidth', 1);
hold on;
for i = 2:size(E1, 2)
    plot(E1(:, i), E2(:, i), 'r', 'LineWidth', 1, 'Color',Colours(i,:));
end
title('State-Space Evolution of E1 and E2');
xlabel('E1 (mV)');
ylabel('E2 (mV)');

%% Isocline of the System

% Define range for E1 and E2
E1_range = linspace(-5, 10, 100);
E2_range = linspace(-5, 10, 100);

% Define function handles for dE1/dt and dE2/dt
theta = @(x) x^2/(k^2+x^2);
dN1dt = @(N1) (-N1 + sigma*sqrt(tau2)*randn())/tau2;
dN2dt = @(N1) (-N1 + sigma*sqrt(tau2)*randn())/tau2;
dE1dt = @(E1, E2) (-E1 + theta(w*E1 - E2 - 0.5) + N1)/tau;
dE2dt = @(E1, E2) (-E2 + theta(E1 - E2 - 0.5) + N2)/tau;

% Compute isoclines
E1_isocline = zeros(length(E1_range), length(E2_range));
E2_isocline = zeros(length(E1_range), length(E2_range));

for i = 1:length(E1_range)
    for j = 1:length(E2_range)
        E1_isocline(i,j) = dE1dt(E1_range(i), E2_range(j));
        E2_isocline(i,j) = dE2dt(E1_range(i), E2_range(j));
    end
end

% Plot isoclines
figure;
contour(E1_range, E2_range, E1_isocline, [0,0], 'LineWidth', 1.5, 'Color', 'blue');
hold on;
contour(E1_range, E2_range, E2_isocline, [0,0], 'LineWidth', 1.5, 'Color', 'red');
title('Isoclines of E1 and E2');
xlabel('E1 (mV)');
ylabel('E2 (mV)');
legend('dE1/dt = 0', 'dE2/dt = 0');

%% Isocline with state space

%%Calculation of E1 and E2
%First elements of matrix E is for values of each trajectory
%Second elements of matrix E is for number of trajectories
num_trajec = 50;
E1 = zeros(length(t),num_trajec); % E1 pre-allocation
E2 = zeros(length(t),num_trajec); % E2 pre-allocation

%defining different random initial point for plotting state space
for i = 1:num_trajec %For each trajectory we need different initial point
    E1(1,i) = 4*rand-2; % E1 starting point
    E2(1,i) = 4*rand-2; % E2 starting point
end

%Stochastic term
N1 = zeros(length(t), num_trajec);
N2 = zeros(length(t), num_trajec);
N1(1, :) = 0;
N2(1, :) = 0;

% Iterative integration
for j = 1:num_trajec

    for i = 2:length(t)
        dN1 = -N1(i-1, j) + sigma*sqrt(tau2*dt)*randn();
        dN2 = -N2(i-1, j) + sigma*sqrt(tau2*dt)*randn();
        dE1 = (-E1(i-1, j) + theta(w*E1(i-1,j) - E2(i-1,j) - 0.5) + N1(i-1,j))*dt/tau;
        dE2 = (-E2(i-1,j) + theta(E1(i-1,j) - E2(i-1, :) - 0.5) + N2(i-1,j))*dt/tau;
        E1(i, j) = E1(i-1, j) + dE1;
        E2(i, j) = E2(i-1, j) + dE2;
        N1(i, j) = N1(i-1,j) + dN1;
        N2(i, j) = N2(i-1, j) + dN2;
    end

end

%Plotting
figure
plot(E1(:,1),E2(:,1), 'k', 'Linewidth', 1)
hold on
for i=2:num_trajec
    plot(E1(:,i),E2(:,i), 'k', 'Linewidth', 1)
end
title('Isoclines of E1 and E2 + State Space');
xlabel('E1 (mV)');
ylabel('E2 (mV)');
% Isocline Calculation
iso_1 = @(E_1,E_2)-E_1/tau + (w*E_1 - E_2 - 0.5)^2/((0.75^2) + (w*E_1 - E_2 - 0.5)^2)/tau;
iso_2 = @(E_1,E_2)-E_2/tau + (E_1 - E_2 - 0.5)^2/((0.75^2) + (E_1 - E_2 - 0.5)^2)/tau;
fimplicit(iso_1,'LineWidth',2) %Plotting nullcline for first E1
hold on
fimplicit(iso_2,'LineWidth',2) %Plotting nullcline for first E2
hold off

%% First passage time from high to low

% Parameters
sigma = 1;
T_end = 50*tau;
threshold = 0.1;

% Initial conditions
E1 = 0.9;
E2 = 0.1;
N1 = 0;
N2 = 0;

% Number of trajectories
n_traj = 100;

% Initialize first-passage time vector
T_fpt = zeros(n_traj, 1);

% Loop over trajectories
for i = 1:n_traj
    
    % Initialize time and state vectors
    t = 0:dt:T_end;
    E1_traj = zeros(size(t));
    E2_traj = zeros(size(t));
    
    % Set initial state
    E1_traj(1) = E1;
    E2_traj(1) = E2;
    
    % Loop over time steps
    for j = 2:length(t)
        
        % Compute noise terms
        dN1 = sigma*sqrt(tau2*dt)*randn();
        dN2 = sigma*sqrt(tau2*dt)*randn();
        zeta = randn();
        
        % Compute derivatives
        dE1 = (-E1_traj(j-1) + theta(w*E1_traj(j-1) - E2_traj(j-1) - 0.5) + N1 + dN1)/tau;
        dE2 = (-E2_traj(j-1) + theta(E1_traj(j-1) - E2_traj(j-1) - 0.5) + N2 + dN2)/tau;
        dN1dt = -N1/tau2 + sigma*sqrt(tau2)*zeta;
        dN2dt = -N2/tau2 + sigma*sqrt(tau2)*zeta;
        
        % Update state variables
        E1_traj(j) = E1_traj(j-1) + dt*dE1;
        E2_traj(j) = E2_traj(j-1) + dt*dE2;
        N1 = N1 + dt*dN1dt;
        N2 = N2 + dt*dN2dt;
        
        % Check if threshold has been crossed
        if E2_traj(j) <= threshold
            T_fpt(i) = t(j);
            break;
        end
        
    end
    
end

% Plot cumulative distribution of first-passage time
figure;
[f, x] = ecdf(T_fpt);
plot(x, f, 'LineWidth', 1.5);
title('Cumulative distribution of first-passage time');
xlabel('Time (ms)');
ylabel('Cumulative probability');

%% First-passage time from low to high!

% Parameters
sigma = 1;
T_end = 50*tau;
threshold = 0.9;

% Initial conditions
E1 = 0.2;
E2 = 0.6;
N1 = 0;
N2 = 0;

% Number of trajectories
n_traj = 100;

% Initialize first-passage time vector
T_fpt = zeros(n_traj, 1);

% Loop over trajectories
for i = 1:n_traj
    
    % Initialize time and state vectors
    t = 0:dt:T_end;
    E1_traj = zeros(size(t));
    E2_traj = zeros(size(t));
    
    % Set initial state
    E1_traj(1) = E1;
    E2_traj(1) = E2;
    
    % Loop over time steps
    for j = 2:length(t)
        
        % Compute noise terms
        dN1 = sigma*sqrt(tau2*dt)*randn();
        dN2 = sigma*sqrt(tau2*dt)*randn();
        zeta = randn();
        
        % Compute derivatives
        dE1 = (-E1_traj(j-1) + theta(w*E1_traj(j-1) - E2_traj(j-1) - 0.5) + N1 + dN1)/tau;
        dE2 = (-E2_traj(j-1) + theta(E1_traj(j-1) - E2_traj(j-1) - 0.5) + N2 + dN2)/tau;
        dN1dt = -N1/tau2 + sigma*sqrt(tau2)*zeta;
        dN2dt = -N2/tau2 + sigma*sqrt(tau2)*zeta;
        
        % Update state variables
        E1_traj(j) = E1_traj(j-1) + dt*dE1;
        E2_traj(j) = E2_traj(j-1) + dt*dE2;
        N1 = N1 + dt*dN1dt;
        N2 = N2 + dt*dN2dt;
        
        % Check if threshold has been crossed
        if E2_traj(j) >= threshold
            T_fpt(i) = t(j);
            break;
        end
        
    end
    
end

% Plot cumulative distribution of first-passage time
figure;
[f, x] = ecdf(T_fpt);
plot(x, f, 'LineWidth', 1.5);
title('Cumulative distribution of first-passage time');
xlabel('Time (ms)');
ylabel('Cumulative probability');

%% 95% of all trajectories remain in the initial state
% Parameters
tau = 10;
k = 0.75;
w = 3.6;
tau2 = 2;
dt = 0.05*tau2;
threshold = 0.9;

% Initial conditions
E1 = 0.1;
E2 = 0.1;

% Number of trajectories
n_traj = 1000;

% Initialize fraction of trajectories in initial state
frac_init = zeros(21, 21);

% Loop over T_end values
for i = 1:21
    T_end = i*tau;
    
    % Loop over sigma values
    for j = 1:21
        sigma = j/10;
        
        % Count fraction of trajectories in initial state
        n_init = 0;
        for k = 1:n_traj
            % Initialize time and state vectors
            t = 0:dt:T_end;
            E1_traj = zeros(size(t));
            E2_traj = zeros(size(t));

            % Set initial state
            E1_traj(1) = E1;
            E2_traj(1) = E2;

            % Loop over time steps
            for l = 2:length(t)

                % Compute noise terms
                dN1 = sigma*sqrt(tau2*dt)*randn();
                dN2 = sigma*sqrt(tau2*dt)*randn();
                zeta = randn();

                % Compute derivatives
                dE1 = (-E1_traj(l-1) + theta(w*E1_traj(l-1) - E2_traj(l-1) - 0.5))/tau + dN1/tau;
                dE2 = (-E2_traj(l-1) + theta(E1_traj(l-1) - E2_traj(l-1) - 0.5))/tau + dN2/tau;

                % Update state variables
                E1_traj(l) = E1_traj(l-1) + dt*dE1;
                E2_traj(l) = E2_traj(l-1) + dt*dE2;

            end

            % Check if at least 95% of trajectories are in initial state
            if sum(E2_traj >= threshold) == 0
                n_init = n_init + 1;
            end

        end
        
        % Compute fraction of trajectories in initial state
        frac_init(i, j) = n_init/n_traj;
        
    end
end

% Plot fraction of trajectories in initial state
figure;
x = 0.1:0.1:2.1;
y = 1:21;
surf(x, y, frac_init);
xlabel('sigma');
ylabel('T\_end (ms)');
zlabel('Fraction of trajectories in initial state');

%% Spectral analysis

% Parameters
T_end = 12;
dt = 0.01;
sigma = 1;
num_realizations = 100;

% System equations
tau = 10;
k = 0.75;
w = 3.6;
tau2 = 2;

% Initial conditions
initial_states = [0.9 0.1; 0.2 0.6];

% Preallocate arrays
%On each rows located each realization, on columns values at the specefic
%time (t)
E1_fft_high = zeros(num_realizations, T_end/dt); 
E2_fft_high = zeros(num_realizations, T_end/dt);
E1_fft_low = zeros(num_realizations, T_end/dt);
E2_fft_low = zeros(num_realizations, T_end/dt);
N1_fft = zeros(num_realizations, T_end/dt);
N2_fft = zeros(num_realizations, T_end/dt);

% Loop over initial states
for i = 1:size(initial_states, 1) %we have two initial state
    % Loop over realizations
    for j = 1:num_realizations
        % Initial conditions
        E = initial_states(i,:);
        N = [randn() randn()];

        % Preallocate arrays for time series
        E1 = zeros(1, T_end/dt);
        E2 = zeros(1, T_end/dt);
        N1 = zeros(1, T_end/dt);
        N2 = zeros(1, T_end/dt);
        %Assigning the initial value to vectors
        E1(1,1) = E(1,1);
        E2(1,1) = E(1,2);
        N1(1,1) = N(1,1);
        N2(1,2) = N(1,2);

        % Time-evolution loop
        for t = 1:(T_end/dt)-1
            zeta = randn();
            E1(t+1) = E1(t) + dt/tau*(-E1(t) + theta(w*E1(t) - E2(t) - 0.5) + N(1));
            E2(t+1) = E2(t) + dt/tau*(-E2(t) + theta(E1(t) - E2(t) - 0.5) + N(2));
            N1(t+1) = N1(t) + dt/tau2*(-N1(t) + sigma*sqrt(tau2)*zeta);
            N2(t+1) = N2(t) + dt/tau2*(-N2(t) + sigma*sqrt(tau2)*zeta);
            N = [N1(t+1) N2(t+1)];
        end

        % Compute and store Fourier transforms
        if initial_states(i,1) == 0.9
            E1_fft_high(j,:) = abs(fft(E1));
            E2_fft_high(j,:) = abs(fft(E2));
        else
            E1_fft_low(j,:) = abs(fft(E1));
            E2_fft_low(j,:) = abs(fft(E2));
        end
        N1_fft(j,:) = abs(fft(N1));
        N2_fft(j,:) = abs(fft(N2));
    end
end

% Compute and plot average spectra
E1_spectrum_high = mean(E1_fft_high, 1);
E2_spectrum_high = mean(E2_fft_high, 1);
E1_spectrum_low = mean(E1_fft_low, 1);
E2_spectrum_low = mean(E2_fft_low, 1);
N1_spectrum = mean(N1_fft, 1);
N2_spectrum = mean(N2_fft, 1);

figure()

plot(E1_spectrum_high,'LineWidth',2)
title('E1/2 spectrum starting in high/low state')
xlabel('Frequency (cycles per sample)')
ylabel('Magnitude')
hold on
plot(E2_spectrum_high,'LineWidth',2)
plot(E1_spectrum_low,'LineWidth',2)
plot(E2_spectrum_low,'LineWidth',2)
legend('E1(high)','E2(high)','E1(Low)','E2(Low)')

%figure for brownian motion
% N1_spectrum = fftshift(N1_spectrum);
% N2_spectrum = fftshift(N2_spectrum);
figure()
plot(N1_spectrum,'LineWidth',2)
hold on
plot(N2_spectrum,'LineWidth',2)
legend('N1','N2')
title('Experimental Spectrum of Brownian Noise N(t)');
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density')
xticks([0 300 600 900 1200])
xticklabels({'-10','-5','0','5','10'})
yticks([0 45 90 135 180])
yticklabels({'0','0.5','1','1.5','2'});

% Diffferntial spectrum
E1_differential = E1_spectrum_low - E1_spectrum_high;
E2_differential = E2_spectrum_low - E2_spectrum_high;
figure()
plot(E1_differential,'LineWidth',2)
hold on
plot(E2_differential,'LineWidth',2)
legend('E1 differential','E2 differential')
title('E1/2 differential spectrum')
xlabel('Frequency (cycles per sample)')
ylabel('Magnitude')

%% Comparison of Theoretical and Simulated Power Spectra of Brownian Noise
% Define parameters
tau_n = 2; % timescale of noise
sigma_n = 1; % amplitude of noise

% Load simulated noise data 
% N1 is a vector of noise values sampled at a fixed time interval
dt = 0.01; % time interval between samples
Fs = 1/dt; % sampling frequency
L = length(N1); % length of signal
t = (0:L-1)*dt; % time vector

% Compute power spectrum of simulated noise
Y = fft(N1);
P2 = abs(Y/L).^2;
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

% Compute theoretical spectrum of Brownian noise
w = linspace(0, 10, 1000);
S = ((tau_n)*(sigma_n^2)) ./ (1 + (tau_n^2)*(w.^2));

% Plot both spectra on the same graph
loglog(f, P1, 'b', w, S, 'r','LineWidth',2);
xlabel('Frequency (Hz)(log)');
ylabel('Power Spectral Density(log)');
title('Comparison of Theoretical and Simulated Power Spectra of Brownian Noise');
legend('Simulated Noise', 'Theoretical Brownian Noise');

