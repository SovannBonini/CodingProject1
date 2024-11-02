%bernoulli s -> P0 = P1 + 1/2 rho V^2;

clear;
clc;
close all;

%% FETCH DATA
data = readmatrix('desktopWT_output_data.csv'); %read the data from the XLSX file
og_data = data; %"dummy array" to extract the flightspeed.xlsx, a way for us to modify the data without modifying the original file

%% CLEAN UP THE DATA SET & STORE IT IN ANSWERS

og_data(isnan(og_data)) = []; %find any NaN values in the data set and delete it
time = og_data(:,1);

%delete unwanted data (very low ouliers)
clean = []; %empty array
threshold = 125; %N -> we want to keep all values of Lift over 25N

row = size(og_data); % we want a single value to represent the size of the dataLift array

%go through each row of the array and check if there are any values under the threshold.
for i = 1:row
    if(og_data(i,:) >= threshold)
        clean = [clean; i]; %keep the index when the value is over the threshold
    end
end

time = time(clean); %new time array corresponding to the cleanLift array (updated dataset)
og_data_clean = og_data(clean, :); %new dataLift array corresponding to the cleanLift array (updated dataset)

pressure = og_data_clean(:,2); %copy the modified dataset into a new array "age" to represent the age column
N = length(pressure); %fetch the number of elements in the array

P0 = 84.2; %assume atmospheric pressure kPa
R = 0.2780;
T = 22+273.15; %assume room temperature
rho = P0 ./ (R .* T); %assume negligible change in rho

unc_P0 = 49.79 / 1000; %uncertainty of pressure atm
unc_P = (unc_P0/P0) * pressure; %uncertainty of Delta P
unc_T = (0.05/22) * T; %fractional uncertainty of temperature;

vi = sqrt((2 .* pressure) ./ rho );

v_max = max(vi); %max velocity
disp(v_max);

%% uncertainty
dpressure = (R .* T) ./ (sqrt(2 .* R .* T .* pressure ./ P0) .* P0);
dR = (T .* pressure) ./ (sqrt(2 .* R .* T .* pressure ./ P0) .* P0);
dT = (R .* pressure) ./ (sqrt(2 .* R .* T .* pressure ./ P0) .* P0);
dP0 = (R .* T .* pressure) ./ (sqrt(2 .* R .* T .* pressure ./ P0) .* P0 .^2);


dv_unnormed = [dpressure .* unc_P, dT .* unc_T, dP0 .* unc_P0];

for i = 1:length(dv_unnormed)
    dV(i) = norm(dv_unnormed(i)); %uncertainty for velocity
end

dV_max = norm(dv_unnormed);
disp(dV_max);

%% Plots
figure();
plot(time, pressure);
figure();
hold on;
plot(pressure,vi);
plot(pressure,vi + 2 * dV, '--k');
plot(pressure,vi - 2 * dV, '--k');
figure();
hold on;
plot(time,vi);
yline(mean(vi) + dV_max, '--r');
yline(mean(vi) - dV_max, '--r');
yline(mean(vi), '--b', 'LineWidth',2);
%plot(time,vi + dV_max, '--k');
%plot(time,vi - dV_max, '--k');