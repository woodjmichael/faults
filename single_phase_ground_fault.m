%% user inputs

Rp = 271.7;
Rf = 1;%109;

c = 1;                          % compensation ratio
e = 1.1;                        % starting voltage
uk0tr = 0.12;                   % transformer short circuit voltage (pu) ?   
X0 = -642.727;                  % capacitive reactance to ground, from powerfactory

Zni = 0;%1e9                    % if non-zero, network is isolated
Znr = 0;%10                     % if non-zero, network is resistance grounded

%% fortescue (0,d,i)

a = j*sqrt(3)/2 - 0.5;
F = [   1   1   1   ;     ...
        1   a^2 a   ;     ...
        1   a   a^2 ];       

%% constants

Vn = 15 * 10^3 / sqrt(3);
In = 25;

e =  e*[1; a^2; a];           % e1, e2, e3 [pu]           
N = 0;
v = e + N;                      % prefault v1, v2, v3 [pu] 

e_odi = F\e;
e_d = e_odi(2);

w = 2 * pi * 50;
Ed = e_d * Vn;                  % for equivalent circuit                   
j = sqrt(-1);                   % prevents a matlab warning 

%% calculated constants

% calculate compensation reactance from known X0
Xe = -X0 / 3 * c;               % coil reactance = 1/3 of cap. reactance
C  = 1 / (w * abs(X0) );        % back-out total capacitance of 3 feeders
Zc = 1/(j*w*C);                 % capacitive impedance

% calculate transformer zero sequencey impedance
Ub = 15.5 * 10^3;               % 15.5 or 15
Sb = 25 * 10^6;                 % 25 or 100
Z0tr = uk0tr * Ub^2 / Sb;       % transformer zero seq impedance

% "zero sequence earthing impedance" (actually its in ABC domain)
Zn = (Rp * j*Xe) ./ (Rp + j*Xe); 

% special: can run as isolated or resistance grounded network
if Zni 
    Zn = Zni;
end 
if Znr 
    Zn = Znr;
end
    
%Zn = 5;    % low resistance grounding
%Zn = 1e15;   % isolated

%%  fault values in 0di components

Z0 = 3*Zn*Zc ./ (3*Zn + Zc) + Z0tr;

I0 = Ed ./ (3*Rf + Z0);

% V0 = I0*Z0
V0  = Zn*Ed ./ (Zn + Rf*(1 + j*w*C * 3*Zn));


%% fault values in ABC (123) components (easy way)

If = 3*I0;

Vres = -3*V0;

Ir = Vres ./ (3*Rp);

Ires_f = -3*V0*(j*w*(C*2/3) + 1/(3*Zn));     
Ires_f = -1* Ires_f;                        % to match Enel img axis                                       

Ires_h = 3*V0*j*w*(C*1/3);
Ires_h = -1*Ires_h;                         % to match Enel img axis

%% fault values in ABC components (Fortescue way)

Vd = 0;     % Zd=0
Vi = 0;     % Zi=0

Ii = I0;       
Id = I0;

% v_postfault = v_prefault + -1*(v_due_to_fault)
v = v - F*[V0; Vd; Vi]/Vn;
vres = sum(v);
v12 = v(1) - v(2);         
v23 = v(2) - v(3);         
v31 = v(3) - v(1);  

% I_postfault = 0 + -1*(I_due_to_fault)
I = -F*[I0; Id; Ii];
Ires = sum(I);

%% display values

% rotate 90º clockwise to match Enel axes
v = v*j;
vres = vres*j;
e = e*j;
Ires = Ires*j;
Ires_f = Ires_f*j;
Ires_h = Ires_h*j;

clc
disp(['Rp (chosen) = ',         num2str(abs(    Rp))])
disp(['Rf (chosen) = ',         num2str(abs(    Rf))])
disp(['Zn = ',                  num2str(abs(    Zn))])
disp([' '])
disp(['Calculated from CESI method_________________________________________'])
disp(['vres = ',                num2str(abs(    Vres)/Vn)])    
disp(['I0 = ',                  num2str(abs(    I0))])
disp(['If = Ir = ',             num2str(abs(    If))])
disp(['Ir = Vres / 3Rp = ',     num2str(abs(    Ir))])
disp(['Ires_f = ',              num2str(abs(    Ires_f))])    
disp(['Ires_h = ',              num2str(abs(    Ires_h))])    

disp([' '])
disp(['Calculated from Fortescue___________________________________________'])
disp(['Ires = ',    num2str(abs(    Ires)),' (',    num2str(angle(  Ires)*180/pi()),'º)'])
disp(['vres = ',    num2str(abs(    vres)),' (',    num2str(angle(  vres)*180/pi()),'º)'])
disp(['e1 = ',      num2str(abs(    e(1))),' (',    num2str(angle(  e(1))*180/pi()),'º), e2 = ', num2str(abs(    e(2))),' (',num2str(angle(  e(2))*180/pi()),'º), e3 = ', num2str(abs(    e(3))),' (',num2str(angle(  e(3))*180/pi()),'º)'])
disp(['v1 = ',      num2str(abs(    v(1))),' (',    num2str(angle(  v(1))*180/pi()),'º), v2 = ', num2str(abs(    v(2))),' (',num2str(angle(  v(2))*180/pi()),'º), v3 = ', num2str(abs(    v(3))),' (',num2str(angle(  v(3))*180/pi()),'º)'])
disp(['v12 = ',     num2str(abs(    v12)),' (',     num2str(angle(  v12)*180/pi()),'º), v23 = ', num2str(abs(    v23)),' (',num2str(angle(   v23)*180/pi()),'º), v31 = ', num2str(abs(    v31)),' (',num2str(angle(   v31)*180/pi()),'º)'])

%% plots

cla
figure(1)
hold on
quiver(0,0,real(vres),imag(vres))
quiver(0,0,real(Ires_f/In),imag(Ires_f/In))     % CESI method
quiver(0,0,real(Ires_h/In),imag(Ires_h/In))     % CESI method

axis([-3 3 -3 3])
legend('vres','ires_f','ires_h')
hold off

figure(2)
hold on
quiver(0,0,real(vres),imag(vres))
quiver(0,0,real(e(1)),imag(e(1)))
quiver(0,0,real(e(2)),imag(e(2)))
quiver(0,0,real(e(3)),imag(e(3)))
quiver(0,0,real(v(1)),imag(v(1)))
quiver(0,0,real(v(2)),imag(v(2)))
quiver(0,0,real(v(3)),imag(v(3)))

axis([-3 3 -3 3])
legend('vres','e1','e2','e3','v1','v2','v3')
hold off
