%% Clear everything from workspace
clear all
close all
clc

%% Read in csv file and plot raw data
filename = 'beam.csv';
Data = csvread(filename,1,0);
t  = Data(:,1);
x  = Data(:,2);
z  = Data(:,3);
ke = Data(:,4);
do_plot(0,t,x,z,ke);

%% Clean Data to when we release the stem + first return
% From the input file for liggghts
deflection_distance = 20.0e-3;
deflection_speed = 1.0;
wait_time = 25.0e-3;
t_release = deflection_distance/deflection_speed + wait_time;
id = find(t>=t_release);

t  = t(id);
x  = x(id);
z  = z(id);
ke = ke(id);

if z(1) < 0
    [~,id_start] = max(z);
else
    [~,id_start] = min(z);
end

id = id_start:length(z);
t  = t(id)-t(id(1));
x  = x(id);
z  = z(id);
ke = ke(id);
do_plot(0,t,x,z,ke);

%% Get frequency information
freq = find_freq(t,z);

%% Normalize deflection data in z direction
z = z./z(1);
t = freq.*t;
do_plot(1,t,z);

%% Find global damping coefficient 
two_PI = 2*pi;
fun = @(b,t)b(1).*exp(-b(2).*t).*cos(b(3)*t);
b00 = [1 0.001 two_PI];
lb = [0 0  0];
ub = [5 1 10];
b0 = lsqnonlin(@(b)fun(b,t)-z,b00,lb,ub);
NLM = fitnlm(t,z,fun,b0);
Global_Bond_Damping = NLM.Coefficients.Estimate(2);
do_plot(2,NLM);