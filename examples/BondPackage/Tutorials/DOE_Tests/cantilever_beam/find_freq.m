function freq = find_freq(t,x,DO_GRAD)

    if nargin < 3
        DO_GRAD = true;
    end
    
    x_org = x;
    if DO_GRAD
        x = gradient(x,t);
    end
    x1 = x(1:end-1);
    x2 = x(2:end);
    id = find(x1.*x2 < 0);
    t_id_1 = t(id);
    t_id_2 = t(id+1);
    x_id_1 = x(id);
    x_id_2 = x(id+1);
    
    m = (x_id_2 - x_id_1)./(t_id_2 - t_id_1);
    t_zero = t_id_1 - x_id_1./m;
    
    std_diff = std(diff(t_zero));
    if std_diff > 0.001
        freq = find_freq(t,x_org,false);
        return
    end
    T = 2*mean(diff(t_zero));
    freq = 1./T;

end