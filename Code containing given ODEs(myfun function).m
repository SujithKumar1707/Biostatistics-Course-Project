function yo = myfun(a,b,c,d,h,delx,maxx)
    temp = maxx/delx;
    yo = zeros(ceil(temp),1);
    xt = 0.0 ; yt =-0.00684124; i = 1;
    while xt < maxx
        yo(i,1) = yt;
        % The below lines represent the functions given
        %yt=yt+ (d*sin(xt^2)-a*sin(h*xt*yt)-b*cos(xt*(yt^2))-c*sin(xt*(yt^3)))*delx;
        yt=yt+(d*sin(xt)-a*xt*sin(h*xt*yt)-b*sin(yt^2)-c*xt*(yt^3))*delx;
        %yt = yt + (d*cos(xt) - a*sinh(h*yt) - b*(yt^2) - c*(yt^3))*delx;
        %yt=yt+ (d*cos(xt^3)-a*cosh(h*yt)-b*xt*(yt)^2-c*sin(xt*(yt)^3))*delx;
        %yt=yt+(d*sin(xt)-a*xt*sin(h*xt*yt)-b*sin(yt^2)-c*xt*(yt^3))*delx;
        %yt = yt + d * cos(xt) -  a * sinh(h*yt) * delx - b * yt.^2 * delx - c * yt.^3 * delx;
        %yt=yt+ (d*cos(xt)-a*cos(h*xt*yt)-b*xt*(yt^2)-c*cos(xt*(yt^3)))*delx;
        %yt=yt+ (d*sin(xt^2)-a*sin(h*xt*yt)-b*cos(xt*(yt^2))-c*sin(xt*(yt^3)))*delx;
        %yt=yt+ (d*cos(xt^3)-a*cosh(h*yt)-b*xt*(yt)^2-c*sin(xt*(yt)^3))*delx;
        %yt=yt+ (d*sin(xt^2)-a*sin(h*xt*yt)-b*cos(xt*(yt^2))-c*sin(xt*(yt^3)))*delx;
       % yt=yt+ (d*sin(xt)-a*sin(h*yt)-b*sin(xt*(yt)^2)-c*sin((yt)^3))*delx;
        %yt=yt+(d*sin(xt)-a*xt*sin(h*xt*yt)-b*sin(yt^2)-c*xt*(yt^3))*delx;
        %yt=yt+ (d*sin(xt)-a*sin(h*yt)-b*sin(xt*(yt)^2)-c*sin((yt)^3))*delx;
        xt = xt + delx;
        i = i + 1;
    end
end