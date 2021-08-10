function [t, y] = euler_expl(f, tspan, y0, opts, param)
    if length(tspan) == 2
        t0 = tspan(1);
        tf = tspan(2);
        h = (tf-t0)/N;
        t = t0:h:tf;
    else
        h = tspan(2)-tspan(1);
        t = tspan;
    end
    
    y=zeros(length(t),length(y0));
    y(1,:) = y0';
    
    for k=1:length(t)-1
        y(k+1,:)=y(k,:)+h*f(t(k), y(k, :)', param)';
    end
end