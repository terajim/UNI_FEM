function [x] = dspace( start, stop , number, ratio )
    if ratio == 1
        x = linspace(start,stop,number);
    else
        a = (stop-start)/((1-ratio.^(number-1))/(1-ratio));
        x(1) = start;
        x(2) = start+a;
        for i = 3:number
            x(i) = x(i-1)+ a*ratio^(i-2);
        end
    end

end

