clc
close all

f0 = 2400; % ������� �������
Vm = 600; % ������������� ��������
Vi = 1200; %�������������� ��������
T = 1/Vm; % ������ ���������� �������
T0 = 1/f0; %������ �������
q = 2 ^ (Vi / Vm); % ���������� ��������
A = sqrt(2/T); % ���������
n = 64;% �������
dT = T0/n;
t = 0:dT:T;


Sc = A*cos(2*pi*f0*t);
Ss = A*sin(2*pi*f0*t);
%sum(Sc.*Sc)*dT
%sum(Ss.*Ss)*dT
%sum(Sc.*Ss)*dT

alpha = 0;
beta = 0;
SNRdB = 15;%1:2:15;%������/��� � ��
SNR = 10.^(SNRdB/10);%������/���
Pe = zeros(1,length(SNRdB));
a = sqrt(SNR*(1-(sqrt(2)/2)));
b = sqrt(SNR*(1+(sqrt(2)/2)));
PeTheor = marcumq(a, b) - (1/2)*(besselj(0, a.*b).*exp(-(a.^2 + b.^2)/2));
nErrMax = 100;
for ns = 1:length(SNRdB) %���� �� ��������� ��������� ������/��� � ��
    nErr = 0; %��������� �������� �������� ����� ������
    nTest = 0; %��������� �������� �������� ����� ��������� 
    sigma = sqrt((sum(Sc.^2))/(2*SNR(ns)));
    thetaPrev = 0; 
    bigThetaPrev = 0;
    phiPrev = 0;%2*pi*rand;
    while nErr < nErrMax %���� ������������� ��� ����� �������� ��������� ������/���
       nTest = nTest + 1;  
       i = floor(q * rand); %����� ��������� i � ��������� [0, q-1]
       deltaTheta = (i*2*pi)/q;
       theta = mod(thetaPrev + deltaTheta, 2*pi);
       thetaPrev = theta;
       epsilon = (2*pi*rand)-pi;
       phi = mod((phiPrev + (alpha*epsilon) + beta), 2*pi);%�����
       phiPrev = phi;
       r = (cos(theta + phi)*Sc) + (sin(theta + phi)*Ss) + (sigma*randn(1, length(t)));
       rc = sum(r.*Sc)*dT;
       rs = sum(r.*Ss)*dT;
       bigTheta = atan2(rs, rc);
       if bigTheta < 0, bigTheta = bigTheta + 2*pi;end
       %bigTheta = mod(bigTheta, 2*pi);
       deltaThetaI = mod(bigTheta - bigThetaPrev, 2*pi);
       if deltaThetaI <= pi/4 && deltaThetaI >= -(pi)/4
           resI = 0;
       elseif deltaThetaI <= 3*pi/4 && deltaThetaI >= pi/4
           resI = 1;
       elseif deltaThetaI <= 5*pi/4 && deltaThetaI >= 3*pi/4
           resI = 2;
       else
           resI = 3;
       end

      
       if resI~=i
            nErr = nErr + 1;
       end
        bigThetaPrev = bigTheta;
    end
    Pe(ns) = nErr/nTest;
    disp([nTest,nErr,SNRdB(ns)]);
end

figure;

grid on
semilogy(SNRdB, PeTheor, 'black', SNRdB, Pe, 'o', 'LineWidth',2);
xlabel('SNRdB') 
ylabel('Pe')
grid on






