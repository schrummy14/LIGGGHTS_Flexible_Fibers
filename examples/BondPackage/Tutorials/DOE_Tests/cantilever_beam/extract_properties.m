function [freq,damp] = extract_properties(Data)

    t  = Data(:,1);
    z  = Data(:,3);

    %% Clean Data to when we release the stem + first return
    % From the input file for liggghts
    deflection_distance = 20.0e-3;
    deflection_speed = 1.0;
    wait_time = 25.0e-3;
    t_release = deflection_distance/deflection_speed + wait_time;
    id = find(t>=t_release);

    t  = t(id);
    z  = z(id);

    if z(1) < 0
        [~,id_start] = max(z);
    else
        [~,id_start] = min(z);
    end

    id = id_start:length(z);
    t  = t(id)-t(id(1));
    z  = z(id);

    %% Get frequency information
    freq = find_freq(t,z);

    %% Normalize deflection data in z direction
    z = z./z(1);
    t = freq.*t;

    %% Find global damping coefficient 
    two_PI = 2*pi;
    fun = @(b,t)b(1).*exp(-b(2).*t).*cos(b(3)*t);
    b00 = [1 0.001 two_PI];
    lb = [0 0  0];
    ub = [5 1 10];
    b0 = lsqnonlin(@(b)fun(b,t)-z,b00,lb,ub);
    NLM = fitnlm(t,z,fun,b0);
    damp = NLM.Coefficients.Estimate(2);
    
end