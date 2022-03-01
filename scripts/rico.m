function [cMA,cAR] = rico(z,B,fc,T)
% filtro ricorsivo reiezione 50 Hz

% input parameters
%	z  (0.01) attenuazione minima alla frequenza fc
%	B  (2-8) larghezza di banda corrispondente alla attenuazione 0.707
% 	fc (50) frequenza di centro banda
% 	T  intervallo di campionamento
% output parameters
%	cMA  filter coefficients (MA part)
%	cAR  filter coefficients (AR part)


b = pi*B*T;
a = b*z;
c1 = -2*(1-a)*cos(2*pi*fc*T);
c2 = (1-a)^2;
c3 = 2*(1-b)*cos(2*pi*fc*T);
c4 = -(1-b)^2;
cMA = [1 c1 c2];
cAR = [1 -c3 -c4];

end