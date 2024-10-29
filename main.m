%{
Author: Charlie Bailey & Sovann Bonini  
Assignment: Coding Project 1: TCAS
Creation Date: 10/14/2024
Inputs: Data_TCAS_A.csv & Data_TCAS_B.csv
Outputs: 
Purpose:
%}

clear;
clc;
close all;

%% FETCH DATA
dataA = readmatrix('Data_TCAS_A.csv'); %read the data from the XLSX file
og_dataA = dataA; %"dummy array" to extract the flightspeed.xlsx, a way for us to modify the data without modifying the original file

dataB = readmatrix('Data_TCAS_B.csv'); %read the data from the XLSX file
og_dataB= dataB; %"dummy array" to extract the flightspeed.xlsx, a way for us to modify the data without modifying the original file

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

%% PLOT DATA

figure();
%plot plane A scatter position
scatter(posX_A,posY_A);

figure();
%plot plane B scatter position
scatter(posX_B,posY_B);

%% LINE OF BEST FIT
%PLOT A
[coeffXA,S_XA] = polyfit(t_A,posX_A,1); %find the coefficients of the linear line of best fit
[coeffYA,S_YA] = polyfit(t_A,posY_A,1); %find the coefficients of the linear line of best fit

%coeff for f(x,t)
bXA = coeffXA(2);
mXA = coeffXA(1);

%coeff for f(y,t)
bYA = coeffYA(2);
mYA = coeffYA(1);

%t_fitA = linspace(min(t_A),max(t_A),length(t_A));

f_XA = bXA + mXA .* t_A; %compute line of best fit for posX_A
f_YA = bYA + mYA .* t_A; %compute line of best fit for posY_A

figure();
hold on;
plot(t_A, f_XA, '--k');
plot(t_A, f_YA, '--b');
plot(f_XA,f_YA); %plot of A position
legend('p(x)_a vs. t_{fitA}','p(y)_a vs. t_{fitA}','p(y)_a vs. p(x)_a','Location','best');

%PLOT B
[coeffXB,S_XB] = polyfit(t_B,posX_B,1); %find the coefficients of the linear line of best fit
[coeffYB,S_YB] = polyfit(t_B,posY_B,1); %find the coefficients of the linear line of best fit

%coeff for f(x,t)
bXB = coeffXB(2);
mXB = coeffXB(1);

%coeff for f(y,t)
bYB = coeffYB(2);
mYB = coeffYB(1); 

%t_fitB = linspace(min(t_B),max(t_B),length(t_B));

f_XB = bXB + mXB .* t_B; %compute line of best fit for posX_B
f_YB = bYB + mYB .* t_B; %compute line of best fit for posY_B

figure();
hold on;
plot(t_B, f_XB, '--k');
plot(t_B, f_YB, '--b');
plot(f_XB,f_YB); %plot of B position
legend('p(x)_b vs. t_{fitB}','p(y)_b vs. t_{fitB}','p(y)_b vs. p(x)_b','Location','best');

%% EXTRAPOLATE

figure(); %before extrapolation
hold on;
plot(f_XA,f_YA);
plot(f_XB,f_YB);

upper_Time = max(t_A) + max(t_B); %extend the timeframe

predictionTime = (min(t_A):upper_Time); %set up a new domain for extrapolatiomn
[extrap_fitXA,deltaXA] = polyval(coeffXA, predictionTime,S_XA); %extract the extrapolated line of best fit with extrapolated values with the associated error
[extrap_fitYA,deltaYA] = polyval(coeffYA, predictionTime,S_YA); %extract the extrapolated line of best fit with extrapolated values with the associated error

[extrap_fitXB,deltaXB] = polyval(coeffXB, predictionTime,S_XB); %extract the extrapolated line of best fit with extrapolated values with the associated error
[extrap_fitYB,deltaYB] = polyval(coeffYB, predictionTime,S_YB); %extract the extrapolated line of best fit with extrapolated values with the associated error

%Extrapolated data
figure();
hold on;
plot(extrap_fitXA,extrap_fitYA); %plot A extrapolated position
plot(extrap_fitXB,extrap_fitYB); %plot B extrapolated position

%% COMPUTE DISTANCE

dU = mXB - mXA;
dV = mYB - mYA;

for i = 1:length(predictionTime)
    dX(i) = extrap_fitXB(i) - extrap_fitXA(i); %xB - xA
    dY(i) = extrap_fitYB(i) - extrap_fitYA(i); %yB - yA 
    D_t(i) = norm([dX(i), dY(i)]); %"dummy" array

    dD(i) = (dX(i).*(dU) + dY(i).*(dV)) ./ (D_t(i));
end

figure();
plot(predictionTime, D_t);

figure();
hold on;
plot(predictionTime, dD);
yline(0);

%% Compute T_CA

xA_0 = bXA;
yA_0 = bYA;

xB_0 = bXB;
yB_0 = bYB;

dX_0 = xB_0 - xA_0;
dY_0 = yB_0 - yA_0;

t_CA = (-(dX_0 .* dU) - (dY_0 .* dV)) ./ ((dU .^2) + (dV .^ 2));
disp(t_CA);


%% Compute position at T_CA

xA_tCA = polyval(coeffXA,t_CA);
yA_tCA = polyval(coeffYA,t_CA);

xB_tCA = polyval(coeffXB,t_CA);
yB_tCA = polyval(coeffYB,t_CA);

dx_tCA = xB_tCA - xA_tCA;
dy_tCA = yB_tCA - yA_tCA;

unnormed = [dx_tCA dy_tCA];
d_tCA = norm(unnormed);

disp(d_tCA);
dD_tCA = (dx_tCA .*(dU) + dy_tCA.* (dV)) ./ (d_tCA);
disp(dD_tCA);

