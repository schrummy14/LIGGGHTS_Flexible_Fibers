clear all
close all
clc

%% Read in data
data = csvread('compression.csv',1,0);
t =  data(:,1);
z = -data(:,2);
f =  data(:,3);
E =  data(:,4);

figure
plotyy(z,f,z,log(E))
drawnow
%% Fit a quadratic to the tail of the data
z_new = z;
f_new = f;
p = 2;
r2_tol = 0.99;
while 1
    A = [ones(length(z_new),1),z_new,z_new.^2];
    b = A\f_new;
    f_prime = A*b;
    mean_f = mean(f_new);
    n = length(f_new);
    top_part = sum((f_new-f_prime).^2)/(n-p-1);
    bot_part = sum((f_new - mean_f).^2)/(n-1);
    adj_r_squared = 1 - top_part/bot_part;
    
    if adj_r_squared < r2_tol || abs(b(1)) > 100
        z_new(1:10) = [];
        f_new(1:10) = [];
    else
        done = 1;
    end
end
figure
plot(z_new,f_new,z_new,f_prime)