clc
close all
clear all

f0 = 2400; % частота несущей
Vm = 600; % модуляционная скорость
Vi = 1200; % информационная скорость
T = 1/Vm; % период следования сигнала
T0 = 1/f0; %период несущей
q = 2 ^ (Vi / Vm); % количество сигналов
A = sqrt(2/T); % амплитуда
n = 64;% отсчеты
dT = T0/n;
t = 0:dT:T;


s = zeros(q, length(t));

%постороение сигнального множества
for i = 0:q-1
    s(i + 1, :) = A*cos(2*pi*f0*t - 2*pi*i/q);
end


Sc = A*cos(2*pi*f0*t);
Ss = A*sin(2*pi*f0*t);
sum(Sc.*Sc)*dT
sum(Ss.*Ss)*dT
sum(Sc.*Ss)*dT

alpha = 0.002;
beta = 0.0005;
SNRdB = 2:2:15;%Сигнал/шум в дБ
SNR = 10.^(SNRdB/10);%Сигнал/шум
Pe = zeros(1,length(SNRdB));

SNRdBTheor = 0:15;
SNRTheor = 10.^(SNRdBTheor/10);
a = sqrt(SNRTheor*(1-(sqrt(2)/2)));
b = sqrt(SNRTheor*(1+(sqrt(2)/2)));
PeTheor = marcumq(a, b) - (1/2)*(besselj(0, a.*b).*exp(-(a.^2 + b.^2)/2));
nErrMax = 100;
for ns = 1:length(SNRdB) %Цикл по значениям отношения сигнал/шум в дБ
    nErr = 0; %Начальное значение счетчика числа ошибок
    nTest = 0; %Начальное значение счетчика числа испытаний 
    nt = sum(sum(s.*s));
    sigma = sqrt(nt/(2*q*SNR(ns)));
    thetaPrev = 0; 
    bigThetaPrev = 0;
    phiPrev = 2*pi*rand;
    while nErr < nErrMax %Цикл моделирования при одном значении отношения сигнал/шум
       nTest = nTest + 1;  
       i = randi(q)-1; %Выбор случайной i в интервале [0, q-1]
       deltaTheta = (i*2*pi)/q;
       theta = mod(thetaPrev + deltaTheta, 2*pi);%mod2pi
       thetaPrev = theta;
       epsilon = (2*pi*rand)-pi;
       phi = mod((phiPrev + (alpha*epsilon) + beta), 2*pi);%сдвиг
       phiPrev = phi;
       r = (cos(theta + phi)*Sc) + (sin(theta + phi)*Ss) + (sigma*randn(1, length(t)));
       rc = sum(r.*Sc)*dT;
       rs = sum(r.*Ss)*dT;
       bigTheta = atan2(rs, rc);
       bigTheta = mod(bigTheta, 2*pi);
       
       resI = mod(round((mod(bigTheta - bigThetaPrev,2*pi)*q) / (2 * pi)),4) ;
       if resI~=i
            nErr = nErr + 1;
       end
       bigThetaPrev = bigTheta;
    end
    Pe(ns) = nErr/nTest;
    disp([nTest,nErr,SNRdB(ns)]);
end

figure;

semilogy(SNRdBTheor, PeTheor, 'black', SNRdB, Pe, 'o', 'LineWidth',2);
xlabel('SNRdB') 
ylabel('Pe')
grid on




