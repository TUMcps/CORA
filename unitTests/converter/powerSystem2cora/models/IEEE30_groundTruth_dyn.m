function f = IEEE30_groundTruth_dyn(x,y,u)

f(1,1)=x(7) - 120*pi;
f(2,1)=x(8) - 120*pi;
f(3,1)=x(9) - 120*pi;
f(4,1)=x(10) - 120*pi;
f(5,1)=x(11) - 120*pi;
f(6,1)=x(12) - 120*pi;
f(7,1)=(2500*x(13))/53 - (100*x(7))/53 + (12000*pi)/53 - (12500*y(1)*sin(x(1)))/53;
f(8,1)=(2500*x(14))/53 - (100*x(8))/53 + (12000*pi)/53 - (12500*y(2)*sin(x(2) - y(31)))/53;
f(9,1)=(2500*x(15))/53 - (100*x(9))/53 + (12000*pi)/53 - (12500*y(3)*sin(x(3) - y(32)))/53;
f(10,1)=(2500*x(16))/53 - (100*x(10))/53 + (12000*pi)/53 - (12500*y(4)*sin(x(4) - y(33)))/53;
f(11,1)=(2500*x(17))/53 - (100*x(11))/53 + (12000*pi)/53 - (12500*y(5)*sin(x(5) - y(34)))/53;
f(12,1)=(2500*x(18))/53 - (100*x(12))/53 + (12000*pi)/53 - (12500*y(6)*sin(x(6) - y(35)))/53;
f(13,1)=u(1) - (1911387046407553*x(7))/36028797018963968 - x(13) + (28670805696113295*pi)/4503599627370496;
f(14,1)=u(2) - (1911387046407553*x(8))/36028797018963968 - x(14) + (28670805696113295*pi)/4503599627370496;
f(15,1)=u(3) - (1911387046407553*x(9))/36028797018963968 - x(15) + (28670805696113295*pi)/4503599627370496;
f(16,1)=u(4) - (1911387046407553*x(10))/36028797018963968 - x(16) + (28670805696113295*pi)/4503599627370496;
f(17,1)=u(5) - (1911387046407553*x(11))/36028797018963968 - x(17) + (28670805696113295*pi)/4503599627370496;
f(18,1)=u(6) - (1911387046407553*x(12))/36028797018963968 - x(18) + (28670805696113295*pi)/4503599627370496;