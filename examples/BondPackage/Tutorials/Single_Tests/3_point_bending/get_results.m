% clear all
% close all
clc

%% Read in data and plot it
data = csvread(fullfile('post','3-point.csv'),1,0);
x = 1000*abs(data(:,2));
f = data(:,3);

id_start = find(f==0);
x(id_start) = [];
f(id_start) = [];
x = x - x(1);

[~,id_max] = max(x);

data_2 = f(1:id_max);

figure
plot(x(1:id_max),f(1:id_max),x(id_max+1:end),f(id_max+1:end))
xlabel ('Displacement (mm)')
ylabel ('Force (N)')
legend ('Compression','Release','Location','NW')
title 'Raw Data'

%% Smooth the data a litte
f = medfilt1(f,15);
figure
plot(x(1:id_max),f(1:id_max),x(id_max+1:end),f(id_max+1:end))
xlabel ('Displacement (mm)')
ylabel ('Force (N)')
legend ('Compression','Release','Location','NW')
title 'Smoothed Data'

%% Find linear portion of the curve

r2_tol = 0.99;
p = 1;
done = 0;
x_new = x(1:id_max);
f_new = f(1:id_max);

while ~done
    A = zeros(length(x_new),p+1);
    A(:,1) = 1;
    for k = 1:p
        A(:,k+1) = x_new.^k;
    end
    b = A\f_new;
    
    f_prime = b(1);
    for k = 1:p
        f_prime = f_prime + b(k+1)*x_new.^k;
    end
    mean_f = mean(f_new);
    n = length(f_new);
    top_part = sum((f_new-f_prime).^2)/(n-p-1);
    bot_part = sum((f_new - mean_f).^2)/(n-1);
    adj_r_squared = 1 - top_part/bot_part;
    if adj_r_squared < r2_tol || b(1) > 0.125
        x_new(end) = [];
        f_new(end) = [];
    else
        done = 1;
    end
end

figure 
plot(x_new,f_new,'k-',x_new,f_prime,'r--')