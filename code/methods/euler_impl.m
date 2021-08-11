function [t,y] = euler_impl(f,tspan,y0,opts, param)
% Explizites Euler-Verfahren für Systeme von DGLs erster Ordnung
% INPUT :   f - Funktion für rechte Seite
%           tspan - Intervall [t_start , t_end]
%           y0 - Anfangswert für y
%           N - Anzahl der Integrationschritte im Zeitintervall [t_start , t_end]
% OUTPUT :  t - Vektor, der diskretisierte Zeit in tspan enthält
%           y - numerische Lösung mit y_i = y(t_i)
    if length(tspan) == 2
        t0 = tspan(1);
        tf = tspan(2);
        h = (tf-t0)/N;
        t = t0:h:tf;
    else
        h = tspan(2)-tspan(1);
        t = tspan;
    end
    N = length(t)-1;
    y=zeros(N+1,length(y0));
    y(1,:)=y0(:)';
    for k = 1 : N

        Y = transpose(y(k,:)) + h * feval( f, t(k), transpose(y(k,:)), param); % Expliziter Euler => initial estimate
        [Y,isConverged]= newton4euler(f,t(k+1),transpose(y(k,:)),Y,h,param);

        if ~ isConverged
            error(['Explizites Euler-Verfahren unterbrach bei Schritt ',num2str(k)])
        end

        y(k+1,:) = transpose(Y);
    end
end

function [Y,isConverged]=newton4euler(f,x,ytranspose,Y,h,param)
    % Funktion für das Newton Verfahren, angepasst für implizites Euler-Verfahren

    RelTol = 1E-3;
    maxit = 10e4;

    isConverged= (1==0);  % startet mit FALSE
    for n=1:maxit
        fValue = feval(f, x, Y,param);
        F = ytranspose + h * fValue - Y;
        dFdY = h * fValue - eye(length(Y));
        increment=dFdY\F;
        Y = Y - increment;
        if norm(increment,2) < RelTol*norm(Y,2)
            isConverged= (1==1);  % wird TRUE hier
            return
        end
    end
end