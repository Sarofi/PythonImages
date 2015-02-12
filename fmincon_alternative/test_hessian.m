function hq=test_hessian(l,lambda)
x1=[0.8147;0.9058];
x2 = [0.1270; 0.9134];
q=1/(4*pi^2*l(1)*l(2))*exp(-0.5*((x1(1)^2+x1(2)^2)/l(1)+(x2(1)^2+x2(2)^2)/l(2)));
hq=zeros(2,2);
hq(1,1)=q*((x1(1)^2+x1(2)^2)^2+8*l(1)^2-8*l(1)*(x1(1)^2+x1(2)^2))/(4*l(1)^2);
hq(2,2)=q*((x2(1)^2+x2(2)^2)^2+8*l(2)^2-8*l(2)*(x2(1)^2+x2(2)^2))/(4*l(2)^2);
hq(1,2)=q*(x1(1)^2+x1(2)^2-2*l(1))/(2*l(1)^2)*(x2(1)^2+x2(2)^2-2*l(2))/(2*l(2)^2);
hq(2,1)=hq(1,2);