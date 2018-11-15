close all
clear
clc

%% What were the DOE variables and levels

Y = [1.00e9 1.00e9 1.00e9 5.00e9 5.00e9 5.00e9 1.1e10 1.1e10 1.1e10]';
B = [0.0000 0.0005 0.0010 0.0000 0.0005 0.0010 0.0000 0.0005 0.0010]';

%% What were the remaining material properties
Diam = 2.83e-3;
rho  = 125;

%% Extract all the data
NumRuns = length(Y);
Data = cell(NumRuns,1);
f_name = @(x)fullfile('results',['run_',num2str(x)],'beam.csv');
for k = 1:NumRuns
    Data{k} = csvread(f_name(k),1,0);
end

%% Calculate the frequency and damping ratio per run
freq = zeros(NumRuns,1);
Bglo = freq;
for k = 1:NumRuns
    [freq(k), Bglo(k)] = extract_properties(Data{k});
end

%% Calculate the undamped natural frequency 
lb = Diam;
area = pi*(Diam/2)^2;
K = Y*area/lb;
vol = 4 * pi * (Diam/2)^3 / 3;
m = vol*rho;
w2 = K/m;
w = sqrt(w2);

%% Run linear regression
LM_1 = stepwiselm([w,B],freq);
LM_2 = stepwiselm([w,B],Bglo);

LM_3 = stepwiselm([freq,Bglo],w);
LM_4 = stepwiselm([freq,Bglo],B);

figure
subplot(2,2,1)
LM_1.plot
subplot(2,2,2)
LM_2.plot
subplot(2,2,3)
LM_3.plot
subplot(2,2,4)
LM_4.plot