function flag = make_DOE_txt(var_names,varargin)

    Num_Vars = length(varargin);
    Num_Runs = (1:length(varargin{1}))';
    
    fid = fopen('DOE_Params.txt','w');
    
    flag = write_var2file(fid,Num_Runs,var_names(1),'int');
    for k = 1:Num_Vars
        flag = flag && write_var2file(fid,varargin{k},var_names(k+1));
    end
    
    fclose(fid);
    
end

function flag = write_var2file(fid,var,var_name,class_type)
    
    if nargin < 4
        class_type = 'double';
    end

    fprintf(fid,'variable %s universe &\n',var_name);
    
    switch class_type
        case 'int'
            for k = 1:length(var)-1
                fprintf(fid,'%i &\n',var(k));
            end
            fprintf(fid,'%i\n\n',var(end));
        
        case 'double'
            for k = 1:length(var)-1
                fprintf(fid,'%25.20e &\n',var(k));
            end
            fprintf(fid,'%25.20e\n\n',var(end));
    end
    
    flag = true;
end
    