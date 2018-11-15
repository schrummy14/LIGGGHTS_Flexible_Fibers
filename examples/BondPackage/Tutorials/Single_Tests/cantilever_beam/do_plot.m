function h = do_plot(type, varargin)

    switch type
        case 0
            t  = varargin{1};
            x  = varargin{2};
            z  = varargin{3};
            ke = varargin{4};
            h = cell(3,1);
            figure
            subplot(3,1,1)
            h{1} = plot(t,x);
            xlabel '(a)'
            ylabel 'x-Displacement'
            grid on
            subplot(3,1,2)
            h{2} = plot(t,z);
            ylabel 'y-Displacement'
            xlabel('(b)')
            grid on
            subplot(3,1,3)
            h{3} = plot(t,ke);
            ylabel 'Kinetic Energy'
            xlabel({'Time (seconds)';'(c)'})
            grid on
        case 1
            t = varargin{1};
            z = varargin{2};
            figure
            h = plot(t,z);
            xlabel 'Normalized Time (t/T)'
            ylabel 'z-Displacement (z/max(|z|))'
            grid on
        case 2
            NLM = varargin{1};
            t = NLM.Variables.x1;
            z = NLM.Variables.y;
            figure
            h = cell(2,1);
            subplot(2,1,1)
            h{1} = plot(t,z,'b-',t,NLM.Fitted,'r--');
            ylabel({'Normalized','Displacement (z/max(|z|))'})
            xlabel '(a)'
            legend('Data','Fitted Model')
            grid on
            subplot(2,1,2)
            h{2} = semilogy(t,abs(NLM.Residuals.Raw));
            ylabel({'Normalized', 'Displacement Error'})
            xlabel({'(b)','Normalized Time (t/T)'});
            grid on
    end
    
end