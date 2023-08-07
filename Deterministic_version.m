clc;
clear;

%This Program is written by Arman Behrad for Exercise No.4 of the Course
%"Spiking Network". Mart.No = 239526

%% Time evolution of System

%time increments and time axis

dt = 0.0001;
N_i = 100000;  %T_max
t_i = (1:N_i)*dt; 

% replicate stimuli along time axis

n = round(1/dt);  %number of replications

tau_E = 0.1;      % time-constant of E_cells activity

w = 2:0.25:4; % w used for Time Evolution

% obtain time-course of K-cell activities
efactor = exp( - dt/tau_E );

%Initial point
%E_0 = 0;
%E_ini = E_0*ones(2,1);
E_ini = [0 ; 0];

E_i = zeros(2,N_i); %peallocating
E_i(:,1) = E_ini ;
%E_matrix = [w , -1 ; 1 , -1];
%connection to matrix [E1_E1 E1_E2 E2_E1 E2_E2]
Constant = [-0.5 ; -0.5];

figure(1)

for j=1:length(w)
%     init_1 = -0.2 + (1.4 + 0.2)*rand; %Capturing two different initial value randomly between -0.2 and 1.4
%     init_2 = -0.2 + (1.4 + 0.2)*rand;
%     E_i(:,1) = [init_1 ; init_2];
E_ini = [0 ; 0];
for i=2:N_i
    %obtain activity of E cells
    E_matrix = [w(j) , -1 ; 1 , -1];
    E_input_i = E_matrix * (E_i(:,i)); 
    E_input_i = E_input_i + Constant;
    
    E_ss_i = ActivationFunction( E_input_i , w(j) ); %Result of E functions at steady states
    
    E_i(:,i) = E_ss_i + (E_ini - E_ss_i) * efactor;
   
    
    E_ini = E_i(:,i);
end

%Plotting the activity of populations in respect to time
subplot(3,3,j) 
plot(t_i,E_i(1,:),t_i,E_i(2,:),'LineWidth',2);
axis([0 10 0 1.2]);
title('Time-evolution of System Activaty');
ylabel('free response')
xlabel('Time')
legend('E_1','E_2')

%Time-evolution in the state space (E1,E2)
%  plot(E_i(1,:)',t_i,t_i,E_i(2,:)','LineWidth',2);
%  axis([-0.5 1.5 -0.5 1.5])
%  %title(['Initial value, X =  ', num2str(init_1), ' Y = ',num2str(init_2),"w = " + w(j)])
%  title("w = " + w(j))
end

 %% Calculation fix points using solve()

x_0 = [0,0 ; 0.5,0.5 ; 1,0.5]; %defining three different initial step which determines three different root for w>2.25 and remain same for w<2.25
tau = 10;
kappa = 0.75;
w = 2:0.25:4; % w used for state space analysis

% %Now we want to calculate fixed point for different values of w
fix_x = zeros(length(w),3); % pre-allocating E1 value of fix point
fix_y = zeros(length(w),3); % pre-allocating E2 value of fix point

for i=1:length(w)
       W=w(i);
       for j=1:3 %evaluating root for each w with three different initial point
         x0 = x_0(j,:);
         fun = @(x)exercise4(x,W); 
         z = fsolve(fun,x0);%solve the equation based on assigned initial point (x_0) and weight (W)
         fix_x(i,j) = z(1,1);
         fix_y(i,j) = z(1,2);
       end
end

%It is not that everythings happens as you expected!
fix_x(2,3) = fix_x(2,2);
fix_y(2,3) = fix_y(2,2);


% State Space Analysis and Plots
t = 0:0.1:10000; % time evolution vector
num_trajec = 50; % number of E trajectories

dt = t(2)-t(1);

%%Calculation of E1 and E2
%Here first elements of matrix E is for 9 different w
%Second elements of matrix E is for values of each trajectory
%Second elements of matrix E is for number of trajectories
%So at the end of day we have 50 trajectory for each w
E1 = zeros(length(w),length(t),num_trajec); % E1 pre-allocation
E2 = zeros(length(w),length(t),num_trajec); % E2 pre-allocation

%defining different random initial point for plotting state space
for i = 1:num_trajec %For each trajectory we need different initial point
    E1(1:length(w),1,i) = 4*rand-2; % E1 starting point
    E2(1:length(w),1,i) = 4*rand-2; % E2 starting point
end

%Now calculating different E based on the different initial points
for i = 1:length(t)-1
    E1(:,i+1,:) = (1-dt/tau)*E1(:,i,:) + dt*phi(w'.*E1(:,i,:) - E2(:,i,:) -0.5, kappa)/tau; % E1 time evolvement
    E2(:,i+1,:) = (1-dt/tau)*E2(:,i,:) + dt*phi(E1(:,i,:) - E2(:,i,:) -0.5, kappa)/tau; % E2 time evolvement
end


% Flow Field Calculation

[E1_q, E2_q] = meshgrid(-6:0.1:6,-6:0.1:6); % set up a grid for quiver

F1 = zeros(size(E1_q,1), size(E1_q,2), length(w)); % gradient field F1 pre-allocation
F2 = zeros(size(E2_q,1), size(E2_q,2), length(w)); % gradient field F2 pre-allocation

for i = 1:length(w)
F1(:,:,i) = gradient_for_E1(E1_q, E2_q, w(i), kappa); % gradient field F1 calculation
F2(:,:,i) = gradient_for_E1(E1_q, E2_q, w(i), kappa); % gradient field F2 calculation
end

figure(2)
for i = 1:9
    subplot(3,3,i) % create a good subplot formatting
    title("w = " + w(i))
    hold on
    for j=1:size(E1,3) % plotting of the E trajectories
        plot(E1(i,:,j),E2(i,:,j), 'k', 'Linewidth', 1)
    end

    % Isocline Calculation
    iso_1 = @(E_1,E_2)-E_1/tau + (w(i)*E_1 - E_2 - 0.5)^2/((0.75^2) + (w(i)*E_1 - E_2 - 0.5)^2)/tau;
    iso_2 = @(E_1,E_2)-E_2/tau + (E_1 - E_2 - 0.5)^2/((0.75^2) + (E_1 - E_2 - 0.5)^2)/tau;
    fimplicit(iso_1,'LineWidth',2) %Plotting nullcline for first E1
    fimplicit(iso_2,'LineWidth',2) %Plotting nullcline for first E2

    plot(fix_x(i,:),fix_y(i,:),'*r','LineWidth',2); %Plottinf fixed point
    quiver(E1_q, E2_q, F1(:,:,i), F2(:,:,i),3) % gradient field
    hold off
    xlabel('E1')
    ylabel('E2')
    xlim([-0.5 1.25])
    ylim([-0.5 1.25])
    set(gca,'FontSize',12)
    %legend('E1 trajectory','E2 trajectory','E1 isocline','E2 isocline','fixed point');
end

%% Bifurcation Diagram

x_0 = [0,0 ; 0.5,0.5 ; 1,0.5]; %defining three different initial step which determines three different root for w>2.25 and remain same for w<2.25
tau = 10;
kappa = 0.75;
w = 2:0.01:4;

for i=1:length(w)
       W=w(i);
       for j=1:3 %evaluatinf root for each w with three different initial point
         x0 = x_0(j,:);
         fun = @(x)exercise4(x,W); 
         z = fsolve(fun,x0);
         fix_x(i,j) = z(1,1);
         fix_y(i,j) = z(1,2);
       end
       plot(w(i),fix_x(i,:),'*b','LineWidth',2); %Plottinf fixed point
       hold on
       xlabel('W');
       ylabel('Es')
end

%% Jacobian Matrix
x_0 = [0.1847,0.5962 ; 0.3162,0.2548 ; 0.9461,0.1415]; %you have to pick up one of this initial step to rech to the identical steady state or you can just pick steady state instead 
tau = 10;
kappa = 0.75;
w =4; %determine your weight
x0 = x_0(3,:); %Pick one of the rows of matrix of x_0
fun = @(x)exercise4(x,w);
[x,fval,exitflag,output,jacobian] =fsolve(fun,x0);
B = [-0.5;-0.5]; %B is the constant input vector size(2,1)
% First run the LinearOrder2 function
LinearOrder2(jacobian,B,x0)
%% Functions
% Activation function of E1/E2
function E_s = ActivationFunction( E, w )   % activation function

E_1 =  ((w*E(1,:) - E(2,:) - 0.5)^2)/((0.75^2) + (w*E(1,:) - E(2,:) - 0.5)^2);
E_2 =  ((E(1,:) - E(2,:) - 0.5)^2)/((0.75^2) + (E(1,:) - E(2,:) - 0.5)^2);

E_s = [E_1 ; E_2];
end

% function of activation function in singular phi
function phi_val = phi(x, kappa)

      phi_val = x.^2./(kappa^2+x.^2);

end

% function to calculate gradient field for E1
function grad = gradient_for_E1(E1, E2, w, kappa)
   grad = -E1 + phi(w*E1- E2 - 0.5, kappa);
end

% function to calculate gradient field for E2
function grad = gradient_for_E2(E1, E2, w, kappa)
   grad = -E2 + phi(E1 - E2 -0.5, kappa);
end

function dxdt = exercise4(x,w)

dxdt = [-x(1) + (((w*x(1) - x(2) - 0.5)^2)/((0.75^2)+(w*x(1) - x(2) - 0.5)^2)); -x(2) + (((x(1) - x(2) - 0.5)^2)/((0.75^2)+((x(1) - x(2) - 0.5)^2)))];

end