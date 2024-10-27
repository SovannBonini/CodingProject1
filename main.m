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
bXA = coeffXA(1);
mXA = coeffXA(2);

%coeff for f(y,t)
bYA = coeffYA(1);
mYA = coeffYA(2);

t_fitA = linspace(min(t_A),max(t_A),length(t_A));

f_XA = bXA + mXA .* t_fitA; %compute line of best fit for posX_A
f_YA = bYA + mYA .* t_fitA; %compute line of best fit for posY_A

figure();
plot(f_XA,f_YA); %plot of A position

%PLOT B
[coeffXB,S_XB] = polyfit(t_B,posX_B,1); %find the coefficients of the linear line of best fit
[coeffYB,S_YB] = polyfit(t_B,posY_B,1); %find the coefficients of the linear line of best fit

%coeff for f(x,t)
bXB = coeffXB(1);
mXB = coeffXB(2);

%coeff for f(y,t)
bYB = coeffYB(1);
mYB = coeffYB(2); 

t_fitB = linspace(min(t_B),max(t_B),length(t_B));

f_XB = bXB + mXB .* t_fitB; %compute line of best fit for posX_B
f_YB = bYB + mYB .* t_fitB; %compute line of best fit for posY_B

figure();
plot(f_XB,f_YB); %plot of B position

%% EXTRAPOLATE

figure(); %before extrapolation
hold on;
plot(f_XA,f_YA);
plot(f_XB,f_YB);

upper_Time = max(t_A)+max(t_B); %extend the timeframe

predictionTimeA = (min(t_A):upper_Time); %set up a new domain for extrapolatiomn
[extrap_fitXA,deltaXA] = polyval(coeffXA, predictionTimeA,S_XA); %extract the extrapolated line of best fit with extrapolated values with the associated error
[extrap_fitYA,deltaYA] = polyval(coeffYA, predictionTimeA,S_YA); %extract the extrapolated line of best fit with extrapolated values with the associated error

predictionTimeB = (min(t_B):upper_Time); %set up a new domain for extrapolatiomn
[extrap_fitXB,deltaXB] = polyval(coeffXB, predictionTimeB,S_XB); %extract the extrapolated line of best fit with extrapolated values with the associated error
[extrap_fitYB,deltaYB] = polyval(coeffYB, predictionTimeB,S_YB); %extract the extrapolated line of best fit with extrapolated values with the associated error

%Extrapolated data
figure();
hold on;
plot(extrap_fitXA,extrap_fitYA); %plot A extrapolated position
plot(extrap_fitXB,extrap_fitYB); %plot B extrapolated position

%% COMPUTE DISTANCE
