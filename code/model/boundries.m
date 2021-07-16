function r = boundries(Za, Zb, Z0, ZT)
%
    
    r =[Za(1)-Z0(2); 
        Za(3)-Z0(4); 
        Za(2)-Z0(1); 
        Za(4)-Z0(3); 
        Zb(5)+1;        % lambda_1(tf)+lambda_0,x
        Zb(2)-ZT(1);
        Zb(7);
        Zb(4)-ZT(2)];

end
