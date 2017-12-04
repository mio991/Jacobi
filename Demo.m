h=1/(2^5)
x=h:h:(2^5-1)*h
f=pi^2*sin(pi*x)
u0=sin(32*x)
u=Helmholz(f',u0')