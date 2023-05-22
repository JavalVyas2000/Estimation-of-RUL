clc;
clear all;
close all;

%Initial Conditions
%Constants
S_n=0.0878;
alpha_n=0.5;
F=96485;
T=298.15;
S_p=0.0878;
R=8314;
alpha=0.5;
kp=0.6*10^(-10);
kn=1.1*10^(-9);
EOCV=4.2;
EODV_valid=3;
EODV_cycle=3.2;
EOCC=0.1;

% Positive Electrode
L_p=7.4*10^(-5);
epsilon_p=0.338;
epsilon_fl_p=0.142;
Rs=2*10^-6;
Cs_max_p=51555;
Ce_0_p=1000;
Ds_p=10^(-13);
De_p=2.5*10^(-10);
sigma_p=10;
kappa_p=2.5;
transfer_no_p=0.2;
n=1;

%Membrane Separator
L_ms=2.5*10^(-5);
epsilon_e_m=0.37;
Ce_0_m=1000;
De_m=2.5*10^(-10);
kappa_m=2.5;
transfer_no_m=0.2;

%Negative Electrode
Ln=7.5*10^(-5);
epsilon_n=0.440;
epsilon_fl_n=0.07;
Rs=2*10^-6;
Cs_max_n=30555;
Ce_0_n=1000;
Ds_n=3.8*10^(-14);
De_n=2.5*10^(-10);
sigma_n=100;
kappa_n=2.5;
transfer_no_n=0.2;
Curr=1.67;

asp=3*(1-epsilon_p-epsilon_fl_p)/Rs;
Jptot=Curr/(asp*S_p*L_p);
Cons_mass_1_p=De_p*epsilon_p^(1.5)/(epsilon_p);

asn=3*(1-epsilon_n-epsilon_fl_n)/Rs; 
Jntot=-Curr/(asn*S_n*Ln);
Cons_mass_1_n=De_n*epsilon_n^(1.5)/(epsilon_n);


kn1 = 1.1*10^(-6);
kn2 = 1.1*10^(-12);
Msei = 0.162; 
psei = 1690;
ksei =5.5*10^(-16); 
kappa_sei = 5*10^(-8);
Mplating =6.94*10^(-3); 
pplating = 0.535 * 106;
kplating = 1.1 * 10^(-18); 
kappa_plating=1.1*10^7;
y0=[51555*0.48 1000 30555*0.83];
% p0=[25361 1000 25361];
% y0=[30555 1000 30555];
Uref_pp=0;
Uref_np=0;
% tspan=1:1:000;
t0=1;
tf=7000;
delta=7000;
Voltage=[];
time=[];
cp_max=51555;
cn_max=30555;
C_sei_loss=0;
V=4.2;
count=0;
for i=1:3
    Voltage=[];
    time=[];
    p_time=count*delta;
    count=count+1;
    [t,y]=ode45(@(t,y) charge(t,y,p_time,cp_max,cn_max,Curr,C_sei_loss,epsilon_p,epsilon_n,F,asp,asn,Jntot,Jptot,kp,kn,alpha,Ln,kn1,kn2,Msei,psei,ksei,kappa_sei,R,T,S_n,Uref_pp,Uref_np,S_p,L_p),[t0 tf],y0);
    disp(y(3,end))
    Voltage=[];
    time=[];
    p_time=count*delta;
    count=count+1;
    [t,y]=ode45(@(t,y) discharge(t,y,p_time,epsilon_p,Curr,epsilon_n,F,asp,asn,Jntot,Jptot,kp,kn,alpha,Ln,kn1,kn2,Msei,psei,kappa_sei,R,T,S_n,Uref_pp,Uref_np),[t0 tf],y0);

end

function [dydt] = discharge(t,y,p_time,epsilon_p,Curr,epsilon_n,F,asp,asn,Jntot,Jptot,kp,kn,alpha,Ln,kn1,kn2,Msei,psei,kappa_sei,R,T,S_n,Uref_pp,Uref_np)
    dydt=zeros(3,1);
    if y(1)<=51555 && y(3)<=30555
        dydt(1) =asp*Jptot/(epsilon_p*F);
        dydt(2)=0;
        dydt(3)=asn*Jntot/(epsilon_n*F);
        theta_p=y(1)/51555;
        theta_n=y(3)/31555;
        if theta_p<1
            Uref_p= 1654107.79310*theta_p^10 - 12495115.2783*theta_p^9 + 42158126.8123*theta_p^8 - 83659025.2732*theta_p^7 + 108125643.253*theta_p^6 - 95100808.008*theta_p^5 + 57644387.325*theta_p^4 - 23776061.484*theta_p^3 + 6386329.45*theta_p^2 - 1008737.242*theta_p + 71156.162;
            Uref_pp=Uref_p;
        else
            Uref_p=Uref_pp;
        end
        if theta_n>0
            Uref_n=(1.2 + 118.2*theta_n^0.5 - 706.07*theta_n + 2217.65*theta_n^1.5 - 1675.13*theta_n^2)/(1.0 + 131.76*theta_n^0.5 - 32.14*theta_n - 746.85*theta_n^1.5 + 15502.95*theta_n^2 - 14213.075*theta_n^2.5);
            Uref_np=Uref_n;
        else
            Uref_n=Uref_np;
        end
        Jp0=12*F*kp*(y(1)*y(2))^0.5;
        Jn0=12*F*kn*(y(3)*y(2))^0.5;
        eta_act_p =(R*T/(alpha*F*100))*asinh(Jptot/(2*Jp0));
        eta_act_n =(R*T/(alpha*F*100))*asinh(Jntot/(2*Jn0));

        J_des_n0= 2*F*(kn1*kn2*y(3)*y(2))^0.5;
        J_des_n = 2*J_des_n0*sinh((0.5*F*eta_act_n)/(R*T));

        % 
        C_sei_loss =(Jntot -J_des_n0)/(Ln*F);
        del_sei_t=abs(Jntot-J_des_n)*Msei/(psei*F*S_n*Ln);
        Rsei=del_sei_t/(10000*kappa_sei);
        % Rsei=del_sei_t/(kappa_sei);
        V = Uref_p - eta_act_p - Uref_n + eta_act_n - Jntot*Rsei;
%         figure(1)
%         scatter(p_time+t,Curr);
%         hold on 
        figure(1)
        scatter(p_time+t,V);
        hold on 
    end
end

function [dydt] = charge(t,y,p_time,cp_max,cn_max,Curr,C_sei_loss,epsilon_p,epsilon_n,F,asp,asn,Jntot,Jptot,kp,kn,alpha,Ln,kn1,kn2,Msei,psei,ksei,kappa_sei,R,T,S_n,Uref_pp,Uref_np,S_p,L_p)
    dydt=zeros(3,1);
    if y(1)<=51555 && y(3)<=30555
        theta_p=y(1)/(cp_max-C_sei_loss);
        theta_n=y(3)/(cn_max-C_sei_loss);
        if theta_p<1
            Uref_p= 1654107.79310*theta_p^10 - 12495115.2783*theta_p^9 + 42158126.8123*theta_p^8 - 83659025.2732*theta_p^7 + 108125643.253*theta_p^6 - 95100808.008*theta_p^5 + 57644387.325*theta_p^4 - 23776061.484*theta_p^3 + 6386329.45*theta_p^2 - 1008737.242*theta_p + 71156.162;
            Uref_pp=Uref_p;
        else
            Uref_p=Uref_pp;
        end
        if theta_n>0
            Uref_n=(1.2 + 118.2*theta_n^0.5 - 706.07*theta_n + 2217.65*theta_n^1.5 - 1675.13*theta_n^2)/(1.0 + 131.76*theta_n^0.5 - 32.14*theta_n - 746.85*theta_n^1.5 + 15502.95*theta_n^2 - 14213.075*theta_n^2.5);
            Uref_np=Uref_n;
        else
            Uref_n=Uref_np;
        end

        Jp0=1*F*kp*(y(1)*y(2))^0.5;
        Jn0=1*F*kn*(y(3)*y(2))^0.5;
        Jsei_n0= 2*F*kn1*(kn2 + ksei)*(y(3)*y(2))^0.5;
        eta_act_p =(R*T/(alpha*F*1000))*asinh(Jptot/(2*Jp0));
        eta_act_n =(R*T/(alpha*F*1000))*asinh(Jntot/(2*Jsei_n0));

        J_des_n0= 2*F*(kn1*kn2*y(3)*y(2))^0.5;
        J_des_n = 2*J_des_n0*sinh((0.5*F*eta_act_n)/(R*T));

        % 
        C_sei_loss =(Jntot -J_des_n0)/(Ln*F);
        del_sei_t=abs(Jntot-J_des_n)*Msei/(psei*F*S_n*Ln);
        Rsei=del_sei_t/(100*kappa_sei);
        dydt(1) =-asp*Curr/(asp*S_p*L_p*10*epsilon_p*F);
        dydt(2)=0;
        dydt(3)=-asn*J_des_n0*Curr/(asp*S_p*L_p*10*epsilon_p*F);
        V = Uref_p - eta_act_p - Uref_n - eta_act_n + Jntot*Rsei;
        % Jntot=(4.2 - Uref_p + eta_act_p + Uref_n + eta_act_n)/Rsei;
        fprintf("%d %d %d %d %d %d\n",V,Uref_p,eta_act_p,Uref_n,eta_act_n,Jntot*Rsei);
            if V>4.2
            V=4.2;
            Curr=(asn*S_n*Ln)*(4.2 - Uref_p + eta_act_p + Uref_n + eta_act_n)/Rsei;
            disp(Curr)
            end
%         figure(1)
%         scatter(p_time+t,Curr);
%         hold on 
        figure(1)
        scatter(p_time+t,V);
        hold on 
    end
end

% function r