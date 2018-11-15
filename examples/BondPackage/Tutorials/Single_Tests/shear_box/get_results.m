close all
clear
clc

%% Load data

data1 = csvread('results1.csv',1,0);
data2 = csvread('results2.csv',1,0);

data = [data1;data2];
time = data(:,1);
kinetic_energy = data(:,2);
bulk_den = data(:,3);
x_dis = data(:,4);
z_dis = data(:,5);
x_for = data(:,6);
z_for = data(:,7);

%% Plot all data
figure
subplot(2,1,1)
plot(x_dis,z_dis)
subplot(2,1,2)
plot(x_dis,x_for)

%% Only look at data from shear
clearvars -except data2
time = data2(:,1);
kinetic_energy = data2(:,2);
bulk_den = data2(:,3);
x_dis = data2(:,4);
z_dis = data2(:,5)-data2(1,5);
x_for = medfilt1(data2(:,6),25);
z_for = medfilt1(data2(:,7),55);

figure
subplot(2,1,1)
plot(x_dis,z_dis)
subplot(2,1,2)
plot(x_dis,x_for)