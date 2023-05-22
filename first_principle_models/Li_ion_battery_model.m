clc
clear all
close all 
param_battery;
cycles=3;
discap=ones(cycles);
no_cycles=1:1:cycles;
Csnmax=Csnmax0;
Cspmax=Cspmax0;
cap_x=1:1:cycles;
Rsei=0;
for k=1:1:cycles
    disp("Cycle_No")
    disp(k)
    I=1.67;
    Jn_tot = -I/(As(2)*S*L(2));        
    Jp_tot =I/(As(1)*S*L(1)); 
    Jtot=[Jp_tot Jn_tot];
    if k==1
        t0=0;
        tf=4000;
        y0=[48816.6103038211 1000 157.338639698784];
    else
        t0=t_discharge(end);
        tf=t0+4000;
        y0=p0;   
    end
    
    %% Initialting Charging with constant current
    options=odeset('Maxstep',4);
    [t, y]=ode45(@charge, [t0 tf],y0,options,I);
    
    for i=1:length(y)
        if (y(i,3)<=Csnmax || y(i,3)>0) && y(i,1)>=0
            y_used(i,1)=y(i,1);
            y_used(i,2)=y(i,2);
            y_used(i,3)=y(i,3);
            t_used(i,1)=t(i,1);
        else
            break;
        end
    end
  %% Calculating the voltages after the concentrations of CC charging are obtanied
  for j=1:1:length(y_used)
    theta_p(j) = y_used(j,1)/Csmax(1);
    Jp(j)=F*K(1)*sqrt(y_used(j,1)*Ce);
    f_sin_p(j)= Jtot(1)/(2*Jp(j));
    Urefp(j) = 1654107.79310*(((theta_p(j))^10))-12495115.2783*(((theta_p(j))^9))+42158126.8123*(((theta_p(j))^8))-83659025.2732*(((theta_p(j))^7))+108125643.253*(((theta_p(j))^6))-95100808.008*(((theta_p(j))^5))+57644387.325*(((theta_p(j))^4))-23776061.484*(((theta_p(j))^3))+6386329.45*(((theta_p(j))^2))-1008737.242*(theta_p(j))+71156.162;
    eta_act_p(j) = (R*T*asinh(f_sin_p(j)))/D;
    
    theta_n(j) = y_used(j,3)/Csmax(2);
    Jn(j)=F*K(2)*sqrt(y_used(j,3)*Ce);
    Jn_sei(j)=F*sqrt(kn1*(kn2+ksei)*y_used(j,3)*Ce);
    f_sin_n(j)= Jtot(2)/(2*Jn_sei(j));
  

    Urefn_num(j)=(1.2+118.2*((theta_n(j))^0.5)-706.07*theta_n(j)+2217.65*((theta_n(j))^1.5)-1675.13*((theta_n(j))^2));
    Urefn_den(j)=1+131.67*((theta_n(j))^0.5)-32.14*theta_n(j)-746.85*((theta_n(j))^1.5)+15502.95*((theta_n(j))^2)-14213.075*((theta_n(j))^2.5);
    Urefn(j)=Urefn_num(j)/Urefn_den(j);
    eta_act_n(j) = (R*T*asinh(f_sin_n(j)))/D;
    
    Vcell(j) = Urefp(j)-eta_act_p(j)-Urefn(j)+eta_act_n(j)-Jtot(2)*(1+10^(-9)*j);  %Rsei taken as such, to fit curve
    time_CCC(j)=t(j);
    if Vcell(j)>=EOCV
        break;
    end
  end
%   figure(1)
%   plot((time_CCC)/3600,Vcell,'Color','b','LineWidth',6);
%   grid on
%   hold on
%   CC_Curr=ones(1,j)*1.67;
%     figure(4)
%     plot((time_CCC/3600),CC_Curr,'Color','b','LineWidth',6);
%   grid on
%   hold on
    %% Initializing the CV mode charging
    x=9000;     
    time_CVC(1)=time_CCC(j);
    for n=1:1:x
        if n==1
            I(n)=1.67;
        end
        Jn_tot = -I(n)/(As(2)*S*L(2));        
        Jp_tot =I(n)/(As(1)*S*L(1)) ;
        Jtot=[Jp_tot Jn_tot];

        if n==1
            ti=time_CCC(j);
            tff=ti+.25;
            y0=[y_used(j,1) 1000 y_used(j,3)]; 
        else
            ti=t(end);
            tff=ti+.25;
            y0=y_cv;
        end
        options=odeset('Maxstep',0.5);
        [t, y_charge]=ode45(@charge_cv, [ti tff],[y0],options,I(n));
        y_cv=[y_charge(end,1) 1000 y_charge(end,3)];

        theta_p = y_cv(1)/Csmax(1);
        Jp=F*K(1)*sqrt(y_cv(1)*Ce);
        f_sin_p= Jtot(1)/(2*Jp);
        Urefp = 1654107.79310*(((theta_p)^10))-12495115.2783*(((theta_p)^9))+42158126.8123*(((theta_p)^8))-83659025.2732*(((theta_p)^7))+108125643.253*(((theta_p)^6))-95100808.008*(((theta_p)^5))+57644387.325*(((theta_p)^4))-23776061.484*(((theta_p)^3))+6386329.45*(((theta_p)^2))-1008737.242*(theta_p)+71156.162;
        eta_act_p = (R*T*asinh(f_sin_p))/D;

        theta_n = y_cv(3)/Csmax(2);
        Jn=F*K(2)*sqrt(y_cv(3)*Ce);
        Jn_sei=F*sqrt(kn1*(kn2+ksei)*y_cv(3)*Ce);
        f_sin_n= Jtot(2)/(2*Jn_sei);

        Urefn_num=(1.2+118.2*((theta_n)^0.5)-706.07*theta_n+2217.65*((theta_n)^1.5)-1675.13*((theta_n)^2));
        Urefn_den=1+131.67*((theta_n)^0.5)-32.14*theta_n-746.85*((theta_n)^1.5)+15502.95*((theta_n^2))-14213.075*((theta_n)^2.5);
        Urefn=Urefn_num/Urefn_den;
        eta_act_n = (R*T*asinh(f_sin_n))/D;

        Res(n) = Urefp-eta_act_p-Urefn+eta_act_n;
        I(n+1)=(((EOCV-(Res(n)))*(As(2)*S*L(2))))/(1+0.000001*n); 
        time_CVC(n+1)=t(end);

        if I(n+1)<=0.05
            break
        end
    end

%     figure(1)
%     m=ones(length(time_CVC),1)*4.2;   
%     plot(time_CVC/3600,m,'Color','b','LineWidth',6)  
%     hold on
%     
%      figure(4)
%     plot((time_CVC/3600),I,'Color','b','LineWidth',6);
%   grid on
%   hold on
  last_curr=I(n+1);
  
  %% Initializing CC mode discharging
    Rsei=1;
    I= -1.67/2;
    Jp_tot = I/(As(1)*S*L(1)); 
    Jn_tot = -I/(As(2)*S*L(2));
    Jtot=[Jp_tot Jn_tot];
    t0=time_CVC(end);
    disp(t0);
    tf=t0+3600*2.7; 
    p0=[y_cv(1)+1000 1000 y_cv(3)-1000];
    
    options=odeset('Maxstep',10);
    [t, y_discharge]=ode45(@discharge, [t0:0.1:tf],[p0],options,I);
    for i=1:length(y_discharge)
        if (y_discharge(i,1)<=Cspmax || y_discharge(i,1)>0) && y_discharge(i,3)>=0 
            y_used1(i,1)=y_discharge(i,1);
            y_used1(i,2)=y_discharge(i,2);
            y_used1(i,3)=y_discharge(i,3);
            t_out_final1=t(i,1);        
        else
            break;
        end
    end
        t_discharge=[];
  %% Calculating the voltages after the concentrations of CC discharging are obtanied
  for j=1:1:length(y_used1)
    theta_p(j) = y_used1(j,1)/Csmax(1);
    Jp(j)=F*K(1)*sqrt(y_used1(j,1)*Ce);
    f_sin_p(j)= Jtot(1)/(2*Jp(j));
    Urefp(j) = 1654107.79310*(((theta_p(j))^10))-12495115.2783*(((theta_p(j))^9))+42158126.8123*(((theta_p(j))^8))-83659025.2732*(((theta_p(j))^7))+108125643.253*(((theta_p(j))^6))-95100808.008*(((theta_p(j))^5))+57644387.325*(((theta_p(j))^4))-23776061.484*(((theta_p(j))^3))+6386329.45*(((theta_p(j))^2))-1008737.242*(theta_p(j))+71156.162;
    eta_act_p(j) = (R*T*asinh(f_sin_p(j)))/D;
    
    theta_n(j) = y_used1(j,3)/Csmax(2);
    Jn(j)=F*K(2)*sqrt(y_used1(j,3)*Ce);
    f_sin_n(j)= Jtot(2)/(2*Jn(j));
    Urefn_num(j)=(1.2+118.2*((theta_n(j))^0.5)-706.07*theta_n(j)+2217.65*((theta_n(j))^1.5)-1675.13*((theta_n(j))^2));
    Urefn_den(j)=1+131.67*((theta_n(j))^0.5)-32.14*theta_n(j)-746.85*((theta_n(j))^1.5)+15502.95*((theta_n(j))^2)-14213.075*((theta_n(j))^2.5);
    Urefn(j)=Urefn_num(j)/Urefn_den(j);
    eta_act_n(j) = (R*T*asinh(f_sin_n(j)))/D;
    
    Vcell1(j) = Urefp(j)-eta_act_p(j)-Urefn(j)+eta_act_n(j)-Jtot(2)*(1+10^(-10)*j);
    t_discharge(j)=t(j);
    if Vcell1(j)<=EODV(2)
        break;
    end
  end
   p0=[y_used1(j,1) 1000 y_used1(j,3)];
   x_axis=[t(1)/3600 t(2)/3600];
   y_axis=[EOCV Vcell1(1)];
   Csei_loss_per=-5.791*10^(-7)*k^2+0.00565*k-0.0001287;
   Csei_loss=Csei_loss_per*Cspmax0/100;
   Csnmax=Csnmax0-Csei_loss;
   Cspmax=Cspmax0-Csei_loss;

%    figure(1)
%     line(x_axis,y_axis,'Color','b','LineWidth',6);
%    hold on
%    plot((t_discharge)/3600,Vcell1,'Color','b','LineWidth',6);
%    hold on

%    if mod(k,50)==0 || k==1
       figure(5)
       plot((t_discharge)/3600,Vcell1,'Color','b','LineWidth',0.1);
       hold on 
%    end
   xlabel("Time (hrs)")
    ylabel("Voltage (V)")
   title("Voltage vs Time")
   
%     figure(4)
%    DisCurr=ones(1,length(t_discharge))*I;
%    plot([t(1)/3600 t(2)/3600],[last_curr I],'Color','b','LineWidth',6);
%    plot(t_discharge/3600,DisCurr,'Color','b','LineWidth',6);
%    plot([t_discharge(end)/3600 t_discharge(end)/3600],[I -2*I],'Color','b','LineWidth',6);
%    hold on
%    xlabel("Time (hrs)")
%     ylabel("Curr (Amp)")
%    title("Current vs Time")

 time_CCC=[];    
 Vcell=[];
 Vcell1=[];
 time_CVC=[]; 
 if k==1
     discap_max=0.835*(t_discharge(end)-t_discharge(1))/3600;
     discap(k,1)=discap_max;
     frac_cap=1;
 else
     discap(k,1)=0.835*(t_discharge(end)-t_discharge(1))/3600;
     frac_cap(k,1)=discap(k,1)/discap_max;
 end
end

figure(3)
plot(no_cycles,discap(:,1),'Color','b','LineWidth',6)
xlabel("No of Cycles")
ylabel("discap")
title("Capacity vs Cycle")

function [ dydt ] = charge_cv(t,y,I)
param_battery;
dydt=zeros(3,1);
Jp_tot =I/(As(1)*S*L(1));
Jn_tot = -I/(As(2)*S*L(2));       
J_deso_n=F*sqrt(kn1*kn2*y(3)*y(2));
J_seio_n=F*sqrt(kn1*(kn2+Ksei)*y(3)*y(2));
dydt(1)=(-As(1)*Jp_tot)/(Es(1)*F);
dydt(2)=0;
dydt(3)=(-As(2)*Jn_tot*J_deso_n)/(Es(2)*F*J_seio_n);
end

function [ dydt ] = charge(t,y,I)
param_battery;
dydt=zeros(3,1);
Jp_tot =I/(As(1)*S*L(1));
Jn_tot = -I/(As(2)*S*L(2));       
J_deso_n=F*sqrt(kn1*kn2*y(3)*y(2));
J_seio_n=F*sqrt(kn1*(kn2+Ksei)*y(3)*y(2));
dydt(1)=(-As(1)*Jp_tot)/(Es(1)*F);
dydt(2)=0;
dydt(3)=(-As(2)*Jn_tot*J_deso_n)/(Es(2)*F*J_seio_n);
end

function [ dydt ] = discharge(t,y,I)
param_battery;
dydt=zeros(3,1);
Jp_tot =I/(As(1)*S*L(1)) ;
Jn_tot = -I/(As(2)*S*L(2));        
dydt(1)=(-As(1)*Jp_tot)/(Es(1)*F);
dydt(2)=0;
dydt(3)=(-As(2)*Jn_tot)/(Es(2)*F);
end
