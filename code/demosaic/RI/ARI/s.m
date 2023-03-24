function y = s(x)

a = -0.5;

if (2 > x) && (x > 1)
	y = a*abs(x)^3-5*a*abs(x)^2+8*a*abs(x)-4*a;
elseif a >= 2
	y = 0;
else
	y = (a+2)*abs(x)^3-(a+3)*abs(x)^2+1;
end
