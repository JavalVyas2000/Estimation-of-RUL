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
        Jp0=10*F*kp*(y(1)*y(2))^0.5;
        Jn0=10*F*kn*(y(3)*y(2))^0.5;
        eta_act_p =(R*T/(alpha*F*100))*asinh(Jptot/(2*Jp0));
        eta_act_n =(R*T/(alpha*F*100))*asinh(Jntot/(2*Jn0));

        J_des_n0= 4*F*(kn1*kn2*y(3)*y(2))^0.5;
        J_des_n = 4*J_des_n0*sinh((0.5*F*eta_act_n)/(R*T));

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
        figure(2)
        scatter(p_time+t,Curr);
        hold on
    end
end