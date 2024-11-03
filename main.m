%{
Author: Charlie Bailey & Sovann Bonini  
Assignment: Coding Project 1 - TCAS (Plane)
Creation Date: 10/14/2024
Inputs: Data_TCAS_A.csv & Data_TCAS_B.csv
Outputs: Figures for recorded positions (for both planes), figures of lines of best fit (for both planes)
figures for distance over time. Answers structure to make it easier for
storing resulting values
Purpose: Model a fundamental aviation problem which is "how to make sure
aircrafts avoid collisions in mid-air?". The student learn to apply their
skills to a real-life problem, making them understand the true value behind
computational methods and modeling.
%}

clear;
clc;
close all;

%% FETCH DATA
dataA = readmatrix('Data_TCAS_A.csv'); %read the data from the Data_TCAS_A.csv file
og_dataA = dataA; %"dummy array" to extract the Data_TCAS_A.csv, a way for us to modify the data without modifying the original file

dataB = readmatrix('Data_TCAS_B.csv'); %read the data from the Data_TCAS_B.csv file
og_dataB= dataB; %"dummy array" to extract the Data_TCAS_B.csv, a way for us to modify the data without modifying the original file

%% CLEAN UP THE DATA SET & STORE IT IN ANSWERS

%DATA SET A
og_dataA(isnan(og_dataA)) = []; %find any NaN values in the data set and delete it
t_A = og_dataA(:,1); %extract time column
posX_A = og_dataA(:,2); %extract X column
posY_A = og_dataA(:,3); %extract Y column
n_A = length(t_A); %fetch the number of elements in the array

pos_A = [t_A, posX_A, posY_A]; %position over time array

%DATA SET B
og_dataB(isnan(og_dataB)) = []; %find any NaN values in the data set and delete it
t_B = og_dataB(:,1); %extract time column
posX_B = og_dataB(:,2); %extract X column
posY_B = og_dataB(:,3); %extract Y column
n_B = length(t_B); %fetch the number of elements in the array

pos_B = [t_B, posX_B, posY_B]; %position over time array

%% LINE OF BEST FIT
%PLOT A
[coeffXA,S_XA] = polyfit(t_A,posX_A,1); %find the coefficients of the linear line of best fit for X_A
[coeffYA,S_YA] = polyfit(t_A,posY_A,1); %find the coefficients of the linear line of best fit for Y_A
[e_xA,d_XA] = polyval(coeffXA,t_A,S_XA); %the point of this line is to fetch the error of the line without extrapolation for X_A
[e_yA,d_YA] = polyval(coeffYA,t_A,S_YA); %the point of this line is to fetch the error of the line without extrapolation for Y_A

%coeff for f(x,t)
bXA = coeffXA(2); %vertical intercept (initial condition xA_0)
mXA = coeffXA(1); %slope of the line (u_A)

%coeff for f(y,t)
bYA = coeffYA(2); %vertical intercept (initial condition yA_0)
mYA = coeffYA(1); %slope of the line (v_A)

f_XA = bXA + mXA .* t_A; %compute line of best fit for posX_A
f_YA = bYA + mYA .* t_A; %compute line of best fit for posY_A

%PLOT B
[coeffXB,S_XB] = polyfit(t_B,posX_B,1); %find the coefficients of the linear line of best fit for X_B
[coeffYB,S_YB] = polyfit(t_B,posY_B,1); %find the coefficients of the linear line of best fit for Y_B
[e_xB,d_XB] = polyval(coeffXB,t_B,S_XB); %the point of this line is to fetch the error of the line without extrapolation for X_B
[e_yB,d_YB] = polyval(coeffYB,t_B,S_YB); %the point of this line is to fetch the error of the line without extrapolation for Y_B

%coeff for f(x,t)
bXB = coeffXB(2); %vertical intercept (initial condition xB_0)
mXB = coeffXB(1); %slope of the line (u_B)

%coeff for f(y,t)
bYB = coeffYB(2); %vertical intercept (initial condition yB_0)
mYB = coeffYB(1); %slope of the line (v_B)

%t_fitB = linspace(min(t_B),max(t_B),length(t_B));

f_XB = bXB + mXB .* t_B; %compute line of best fit for posX_B
f_YB = bYB + mYB .* t_B; %compute line of best fit for posY_B

%% EXTRAPOLATE
upper_Time = max(t_A) + max(t_B); %extend the timeframe using the respective maximum values of time from A & B

predictionTime = (min(t_A):upper_Time); %set up a new domain for extrapolation
[extrap_fitXA,deltaXA] = polyval(coeffXA, predictionTime,S_XA); %extract the extrapolated line of best fit with extrapolated values with the associated error for X_A
[extrap_fitYA,deltaYA] = polyval(coeffYA, predictionTime,S_YA); %extract the extrapolated line of best fit with extrapolated values with the associated error for Y_A

[extrap_fitXB,deltaXB] = polyval(coeffXB, predictionTime,S_XB); %extract the extrapolated line of best fit with extrapolated values with the associated error for X_B
[extrap_fitYB,deltaYB] = polyval(coeffYB, predictionTime,S_YB); %extract the extrapolated line of best fit with extrapolated values with the associated error for Y_B

%% COMPUTE DISTANCE

dU = mXB - mXA; %compute the difference in X velocities
dV = mYB - mYA; %compute the difference in Y velocities

%for loop to ensure the distance is compute at each instance, not only at the last instance of the extended time frame
for i = 1:length(predictionTime)
    dX(i) = extrap_fitXB(i) - extrap_fitXA(i); %xB - xA at each instance
    dY(i) = extrap_fitYB(i) - extrap_fitYA(i); %yB - yA at each instance
    D_t(i) = norm([dX(i), dY(i)]); %array of all distances within the extended time frame

    dD(i) = (dX(i).*(dU) + dY(i).*(dV)) ./ (D_t(i)); %compute the derivative of distance at each instance
end

%% Compute T_CA

%re-name the initial condition for X_A and Y_A within the code (semantics)
xA_0 = bXA;
yA_0 = bYA;

%re-name the initial condition for X_B and Y_B within the code (semantics)
xB_0 = bXB;
yB_0 = bYB;

%compute the difference between initial conditions for A and B
dX_0 = xB_0 - xA_0;
dY_0 = yB_0 - yA_0;

%compute t_CA using the formula
t_CA = (-(dX_0 .* dU) - (dY_0 .* dV)) ./ ((dU .^2) + (dV .^ 2));
disp(t_CA); %display t_CA
answers.tCA = t_CA; %store in a structure for easier check

%% Compute position at T_CA
[xA_tCA,deltaXAtCA] = polyval(coeffXA,t_CA,S_XA); %extrapolate the X position at t_CA for plane A
[yA_tCA,deltaYAtCA] = polyval(coeffYA,t_CA,S_YA); %extrapolate the Y position at t_CA for plane A

[xB_tCA,deltaXBtCA] = polyval(coeffXB,t_CA,S_XB); %extrapolate the X position at t_CA for plane B
[yB_tCA,deltaYBtCA] = polyval(coeffYB,t_CA,S_YB); %extrapolate the Y position at t_CA for plane B

%difference between the respective positions at t_CA
dx_tCA = xB_tCA - xA_tCA;
dy_tCA = yB_tCA - yA_tCA;

%compute the distance between planes at t_CA
D_tCA_unnormed = [dx_tCA dy_tCA];
D_tCA = norm(D_tCA_unnormed);

disp(D_tCA); %display D at t_CA

%check if the derivative at t_CA is infinitely close to zero (or even equal to zero)
dD_tCA = (dx_tCA .*(dU) + dy_tCA.* (dV)) ./ (D_tCA); %compute the derivative of D at t_CA
answers.dDtCA = dD_tCA; % store value in the structure

%% Computing Uncertainty With General Method

%reassigning variables for easier reading (because partials are really hard to analyze when written in text format in MATLAB
%initial positions
xA0 = xA_0;
yA0 = yA_0;
xB0 = xB_0;
yB0 = yB_0;

%Velocities
uA = mXA;
uB = mXB; 
vA = mYA;
vB = mYB;

Delta = length(predictionTime)*sum(predictionTime.^2) - sum(predictionTime) ^2; %computing the least linear square fit Delta

%uncertainty for the velocities using the formula from the least linear squares slide
delta_uA = deltaXA(1) * sqrt((length(predictionTime))/(Delta));
delta_vA = deltaYA(1) * sqrt((length(predictionTime))/(Delta));
delta_uB = deltaXB(1) * sqrt((length(predictionTime))/(Delta));
delta_vB = deltaYB(1) * sqrt((length(predictionTime))/(Delta));

% 1. Partial derivative with respect to xA,0
dxA0 = (uB - uA) ./ ((uB - uA).^2 + (vB - vA) .^2);
% 2. Partial derivative with respect to yA,0
dyA0 = (vB - vA) ./ ((uB - uA).^2 + (vB - vA) .^2);
% 3. Partial derivative with respect to xB,0
dxB0 = -(uB - uA) ./ ((uB - uA).^2 + (vB - vA) .^2);
% 4. Partial derivative with respect to yB,0
dyB0 = -(vB - vA) ./ ((uB - uA).^2 + (vB - vA) .^2);
% 5. Partial derivative with respect to uA
duA = ((xB0 - xA0) .* ((vB - vA).^2 -(uA - uB).^2) - 2 .* (uB - uA) .* (yB0 - yA0) .* (vB - vA)) ./ (((uB - uA).^2 + (vB - vA).^2).^2);
% 6. Partial derivative with respect to vA
dvA = ((yB0 - yA0) .* ((vB - vA).^2 -(uA - uB).^2) - 2 .* (uB - uA) .* (xB0 - xA0) .* (vB - vA)) ./ (((uB - uA).^2 + (vB - vA).^2).^2);
% 7. Partial derivative with respect to uB
duB = ((xB0 - xA0) .* ((uA - uB).^2 - (vB - vA).^2) + 2 .* (uB - uA) .* (yB0 - yA0) .* (vB - vA)) ./ (((uB - uA).^2 + (vB - vA).^2).^2);
% 8. Partial derivative with respect to vB
dvB = ((yB0 - yA0) .* ((uA - uB).^2 - (vB - vA).^2) + 2 .* (uB - uA) .* (xB0 - xA0) .* (vB - vA)) ./ (((uB - uA).^2 + (vB - vA).^2).^2);

%delta_t_cA using the error propagation formula
delta_tcA = sqrt((dxA0 .* deltaXA(1)).^2 + (dyA0 .* deltaYA(1)).^2 + (dxB0 .* deltaXB(1)).^2 + (dyB0 .* deltaYB(1)).^2 + (duA .* delta_uA).^2 + (dvA .* delta_vA).^2 + (duB .* delta_uB).^2 + (dvB .* delta_vB).^2); 
answers.uncTCA = delta_tcA; %store value in structure

%% DISTANCE UNCERTAINTY  

%compute the uncertainty for the distance at each instance using the general method
dxA = ((-1)*dX) ./ D_t;
dxB = (dX) ./ D_t;
dyA = ((-1)*dY) ./ D_t;
dyB = (dY) ./ D_t;

%uncertainty of distance at each point
delta_D = sqrt((dxA .* deltaXA).^2 + (dxB .* deltaXB).^2 + (dyA .* deltaYA).^2 + (dyB .* deltaYB).^2);

%re-calculate the partials for a specific point (t = t_ca) -> we re-did the formulas once again to avoid potential coding mistakes. It is hard-coding in a way, but rather be safe.
dxA = (xA_tCA - xB_tCA) ./ D_tCA;
dxB = (xB_tCA - xA_tCA) ./ D_tCA;
dyA = (yA_tCA - yB_tCA) ./ D_tCA;
dyB = (yB_tCA - yA_tCA) ./ D_tCA;

%re-compute the uncertainty at t = t_ca
delta_D_tCA = sqrt((dxA .* deltaXAtCA).^2 + (dxB .* deltaXBtCA).^2 + (dyA .* deltaYAtCA).^2 + (dyB .* deltaYBtCA).^2);

disp(delta_D_tCA); %display 
answers.uncDTCA = delta_D_tCA; %store in structure

%% WARNING

%compute the advisories 
if D_tCA < 3.3 && D_tCA > 2.0
    disp('Caution: Aircraft within 20-48s. Escape Maneuver recommended');
elseif D_tCA < 2.0
    disp('Warning: Collision expected within 15-35s. Escape maneuver needed in the Z direction');
elseif D_tCA > 3.3
    disp('Current Path is clear');
end

upBound_DtCA = D_tCA + delta_D_tCA;
lowBound_DtCA = D_tCA - delta_D_tCA;

stringUpBound = ['Maximum Closest Distance: ', num2str(upBound_DtCA), 'nmi'];
stringLowBound = ['Minimum Closest Distance: ', num2str(lowBound_DtCA), 'nmi'];

disp([newline stringUpBound]); %display maximum closest distance
disp(stringLowBound); %display minimum closest distance
disp('Lowest chance of collision with unchanged course: 68%'); %display chance of collision with unchanged course within +/- 1 standard deviation from D at t=t_ca

%% PLOTTING FIGURES 

%Plot Plane A position
fig1 = figure(); 

subplot(2,1,1); %create a subplot for the X position of Plane A
hold on;
scatter(t_A,posX_A, 'blue');
scatter(t_A,posY_A, 'red');
title('X & Y vs. Time of Plane B');
legend('p(x) vs. t','p(y) vs. t','Location','best');
xlabel('Position (nmi)');
ylabel('time (s)');

subplot(2,1,2); %create a subplot for the Y position of Plane A
hold on;
scatter(posX_A,posY_A, 'black');
legend('p(y) vs. p(x)');
title('Y vs. X of Plane B');
xlabel('p(X) (nmi)');
ylabel('p(Y) (nmi)');
saveas(fig1,"figure 1.png"); %save figure

%Plot Plane B position
fig2 = figure();
hold on;

subplot(2,1,1); %create a subplot for the X position of Plane B
hold on;
scatter(t_B,posX_B, 'blue');
scatter(t_B,posY_B, 'red');
legend('p(x) vs. t','p(y) vs. t','Location','best');
title('X & Y vs. Time of Plane B');
xlabel('Position (nmi)');
ylabel('time (s)');

subplot(2,1,2); %create a subplot for the Y position of Plane B
hold on;
scatter(posX_B,posY_B, 'black');
legend('p(y) vs. p(x)');
title('Y vs. X of Plane B');
xlabel('p(X) (nmi)');
ylabel('p(Y) (nmi)');
saveas(fig2,"figure 2.png"); %save figure

%Plot Plane A Line of best fit
fig3 = figure();

subplot(3,1,1); %create a subplot X(t) line of best fit of Plane A
hold on;
plot(t_A, f_XA, 'k');
%1 standard deviation away from the line of best fit
plot(t_A, f_XA + d_XA, '--b');
plot(t_A, f_XA - d_XA, '--b');
title('p(x)_a vs. t_{fitA}');
xlabel('time (s)');
ylabel('X Position (nmi)');
legend('p(x)','p(x) + \sigma_x','p(x) - \sigma_x','Location','best');

subplot(3,1,2); %create a subplot Y(t) line of best fit of Plane A
hold on;
plot(t_A, f_YA, 'k');
%1 standard deviation away from the line of best fit
plot(t_A, f_YA + d_YA, '--b');
plot(t_A, f_YA - d_YA, '--b');
title('p(y)_a vs. t_{fitA}');
xlabel('time (s)');
ylabel('Y Position (nmi)');
legend('p(y)','p(x) + \sigma_y','p(y) - \sigma_y','Location','best');

subplot(3,1,3); %create a subplot Y(t) vs. X(t) lines of best fit of Plane A
hold on;
plot(f_XA,f_YA); %plot of B position
title('p(y)_a vs. p(x)_a');
xlabel('X Position (nmi)');
ylabel('Y Position (nmi)');
saveas(fig3,"figure 3.png"); %save figure

%Plot Plane B Line of best fit
fig4 = figure();

subplot(3,1,1); %create a subplot X(t) line of best fit of Plane B
hold on;
plot(t_B, f_XB, 'k');
%1 standard deviation away from the line of best fit
plot(t_B, f_XB + d_XB, '--b'); 
plot(t_B, f_XB - d_XB, '--b');
title('p(x)_b vs. t_{fitB}');
xlabel('time (s)');
ylabel('X Position (nmi)');
legend('p(x)','p(x) + \sigma_x','p(x) - \sigma_x','Location','best');

subplot(3,1,2); %create a subplot Y(t) line of best fit of Plane B
hold on;
plot(t_B, f_YB, 'k');
%1 standard deviation away from the line of best fit
plot(t_B, f_YB + d_YB, '--b');
plot(t_B, f_YB - d_YB, '--b');
title('p(y)_b vs. t_{fitB}');
xlabel('time (s)');
ylabel('Y Position (nmi)');
legend('p(y)','p(x) + \sigma_y','p(y) - \sigma_y','Location','best');

subplot(3,1,3); %create a subplot Y(t) vs. X(t) lines of best fit of Plane B
hold on;
plot(f_XB,f_YB); %plot of B position
title('p(y)_b vs. p(x)_b');
xlabel('X Position (nmi)');
ylabel('Y Position (nmi)');
saveas(fig4,"figure 4.png"); %save figure

%Plot Distance
fig5 = figure();
hold on;
plot(predictionTime, D_t, 'r','LineWidth',2);
%1 standard deviation away from D(t)
plot(predictionTime, D_t + delta_D, '--b'); 
plot(predictionTime, D_t - delta_D, '--b');
xline(t_CA, 'k', 'LineWidth',2); %plot a vertical line at D(t_CA)
%1 standard deviation away from D(t_CA)
xline(t_CA + delta_tcA, '--p');
xline(t_CA - delta_tcA, '--p');
yline(D_tCA, '-c'); %horizontal line at D_tCA to show that D(t_CA) is when dD_tCA approximately equals 0
legend('D(t)','D(t) + \sigma_d','D(t) - \sigma_d','x = t_{CA}', 't_{CA} + \sigma_{t_CA}', 't_{CA} - \sigma_{t_CA}','D_tCA','Location','best');
title('Distance vs. Time');
xlabel('Time (s)');
ylabel('Distance (nmi)');
saveas(fig5,"figure 5.png"); %save figure
