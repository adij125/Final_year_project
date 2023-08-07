
%% SImulink outputs 
% v_simu=out.relay_voltage
% seq_simu=out.relay_voltage_a
% i_simu=out.relay_current_a 
% i_grid=out.current_grid
% imp=(v_simu(1)-v_simu(2))/i_simu(1) 


% 
% fl=i_simu-i_grid
%% Base parameters 
Vb = 315e3;
Sb = 300e6;
Zb = Vb^2/Sb;
Ib = Sb/(sqrt(3)*Vb);
freq=60;
w = 2*pi*freq; 
%% Fault Params
m=0.75;  
Rf = 10;
Fault_Type = 1; % LG = 1, LL = 2; 

     % Line Parameters Model
    leng = 50; %km
    
    R10 = [ 0.02715    0.25155];
    L10 = [0.00100    0.00364];
    C10 = [ 1.16622e-08  6.99588e-09]; 

    ZL1 = (2.7362 + 1i*37.90822)*leng/100;
    ZL2 = ZL1;
    ZL0 = (25.711 + 1i*138.85)*leng/100;
    ZL = [ZL1; ZL2; ZL0];
    K0 = (ZL0/ZL1-1);
    


 

% Grid Impedance
Psc = 5000e6; % Short Circuit level
Zs_abs = Vb^2/Psc;
XR_ratio = 10;
Rs = Zs_abs/sqrt(1+XR_ratio^2);
Ls = (Rs*XR_ratio)/w;
Zs = Rs + 1i*Ls*w; 

% Transformer Winding Impedance
Rt_pu = 0;
Lt_pu = 0.04;
Zt = 2*(1i*Lt_pu)*Zb; 

% Mho Setting
reach = 0.8;


%phase

limit_LG=0; 
limit_LL=0;
V2=1;
% for z=1:361 
% q=-83.4;
%% Phase terms for simulink model
q=-79.3;
PhaseA=-30 + q-1; 
PhaseB=210 + q-1; 
PhaseC=90 + q-1;

%% Theory
for z=1:361
curr = 1*exp(1i*(z-1)*pi/180);
% This term is used in angle injection to distinguish injected current before a resistive current angle is found: curr_th=1;
% curr_th=1;
% for z=1:100
Rfault = Rf/Zb;
Transformer_imp=Zt/Zb; 
line_impedance1=ZL1/Zb; 
line_impedance0=ZL0/Zb; 
S2=Zs/Zb;

if Fault_Type == 1
%Line to ground fault
V_th=(curr*((1-m)*line_impedance1+S2))+V2; 
Z_1_thev=(S2+(1-m)*line_impedance1);
Z_0_thev=(1/(1/(S2+(1-m)*line_impedance0)+1/(Transformer_imp+m*line_impedance0))); 
Z_2_thev=Z_1_thev; 

I_a1=V_th/(Z_1_thev+Z_2_thev+Z_0_thev+3*Rfault); 
I_fault=3*I_a1;
%% The sets of current angles we inject when using zero sequence currents to make relay errors resistive
% curr(1)=1*exp(1i*(-41.9)*pi/180);
% curr(2)=1*exp(1i*(-40.4025)*pi/180);
% curr(3)=1*exp(1i*(-38.96)*pi/180);
% curr(4)=1*exp(1i*(-37.57)*pi/180);
% curr(5)=1*exp(1i*(-36.23)*pi/180); 
% curr(6)=1*exp(1i*(-34.95)*pi/180);
% curr(7)=1*exp(1i*(-33.71)*pi/180);
% curr(8)=1*exp(1i*(-32.521)*pi/180);
% curr(9)=1*exp(1i*(-31.386)*pi/180);
% curr=1*exp(1i*(70)*pi/180);

I_r0=I_a1*(S2+(1-m)*line_impedance0)/(Transformer_imp+S2+line_impedance0) 
%% Our relay voltage assuming a current angle of 0 
% V_r1_test=(1*m*line_impedance1)+(V_th-I_a1*Z_1_thev);


%% Code for angle injections mimicking synchronous generstor operation
% V_synch=1*(Transformer_imp +S2+line_impedance1)+V2-(V_r1_test); 

% curr=1*exp(1i*((angle(V_synch)-pi/2)));

I_r1=curr; 
I_r2=0; 

V_th=(curr*((1-m)*line_impedance1+S2))+V2; 
I_a1=V_th/(Z_1_thev+Z_2_thev+Z_0_thev+3*Rfault); 
I_r0=I_a1*(S2+(1-m)*line_impedance0)/(Transformer_imp+S2+line_impedance0);
I_AR=I_r1+I_r2+I_r0; 

V_r1=(I_r1*m*line_impedance1)+(V_th-I_a1*Z_1_thev); 
V_r2=-I_a1*Z_2_thev;
V_r0=I_r0*m*line_impedance0-I_a1*Z_0_thev; 
V_AR=V_r1+V_r2+V_r0;

Apparent_Impedance1(z)=V_AR/(I_AR+K0*I_r0); 
%% Code for Simulink relay impedance 
% true_imp=v_simu(1)/(i_simu(1)+K0*seq_simu(3));

%% Code that allows detection of fault resistance where relay impedance falls outside mho
% if angle(reach*line_impedance1-Apparent_Impedance1(z))-angle(Apparent_Impedance1(z))>=(pi/2) && limit_LG==0
% limit_LG=Apparent_Impedance1(z) 
% limit_Rf=Rfault*Zb
% end

%% If statements to allow plotting of multiple relay impedances for varying line parameters
% if p==1 
%     temp1(z)=Apparent_Impedance1(z); 
% elseif p==2 
%     temp2(z)=Apparent_Impedance1(z); 
% elseif p==3 
%     temp3(z)=Apparent_Impedance1(z); 
% elseif p==4 
%     temp4(z)=Apparent_Impedance1(z);
% 
% end 


end
%% Line to line 
if Fault_Type==2
Z_1_thev=(S2+(1-m)*line_impedance1);
Z_2_thev=Z_1_thev; 

V_th=(curr*((1-m)*line_impedance1+S2))+V2; 
 
I_a1=V_th/(Z_1_thev+Z_2_thev+Rfault);

I_fault=-i*sqrt(3)*I_a1; 



V_r2=(I_a1*(Rfault+Z_2_thev)*(Z_2_thev/(Z_2_thev+Rfault))); 

%% Angle injection using sequence relay voltage to estimate fault current angle
% curr = 1*exp(1i*(angle(V_r2)-pi/2));
% curr = 1*exp(1i*(-83.4*pi/180));

%% Angle injection using synchronous generator mimicry
% V_r1_test=(1*m*line_impedance1+(V_th-(I_a1*Z_1_thev))); 
% V_synch=1*(line_impedance1+Transformer_imp+S2)+V2-(V_r1_test);
% curr=1*exp(1i*((angle(V_synch))-pi/2));

V_th=(curr*((1-m)*line_impedance1+S2))+V2; 
 
I_a1=V_th/(Z_1_thev+Z_2_thev+Rfault);
V_r1=(curr*m*line_impedance1+(V_th-(I_a1*Z_1_thev))); 
V_r2=(I_a1*(Rfault+Z_2_thev)*(Z_2_thev/(Z_2_thev+Rfault))); 

I_fault=-i*sqrt(3)*I_a1; 
V_AR=V_r1+V_r2;

I_AR=curr; 

Apparent_Impedance2(z)=(V_r1-V_r2)/curr;

% if angle(reach*line_impedance1-Apparent_Impedance2(z))-angle(Apparent_Impedance2(z))>=(pi/2) && limit_LL==0
% limit_LL=Apparent_Impedance2(z) 
% R_limit=Rfault*Zb
% end

% if p==1 
%     temp1l(z)=Apparent_Impedance2(z); 
% elseif p==2 
%     temp2l(z)=Apparent_Impedance2(z); 
% elseif p==3 
%     temp3l(z)=Apparent_Impedance2(z); 
% elseif p==4 
%     temp4l(z)=Apparent_Impedance2(z);
% 
% end 

end






% end 
end
% error=imp-test_imp; 
% percent_error=abs(error)*100/abs(imp)
% end
hold on
%% Plotting true position of fault along line
true=plot(m*line_impedance1,'x','DisplayName','position of fault','Color','black');

%% Plotting the mho
theta=-2*pi:0.00001:2*pi;
x=(reach/2)*(abs(line_impedance1))*(exp(1i*(theta))+exp(1i*angle(line_impedance1))); 
%% Plotting the apparent impedances where 1 is the line to ground impedance and 2 is the line to line impedance

plot(x,'Color','black') 

plot([0, 1i*reach*abs(line_impedance1)-0.02,1i*reach*abs(line_impedance1)+abs(line_impedance1),abs(line_impedance1)-0.02,0],'Color','black');
% end 
test=[0,line_impedance1];
plot(test,'--');
imp1=plot(Apparent_Impedance1,'DisplayName','L-G fault relay impedance (Script)','Color','blue');
%% Relay impedances obtained from Simulink
% berlin=plot(0.0174+0.0417i,'*','DisplayName','Simulink relay Impedances','Color','magenta'); 
% plot(0.0219+0.0444i,'*','Color','red');  
% plot(0.0172+0.047i,'*','Color','#A2142F');  
% plot(0.0158+0.0442i,'*','Color','black'); 
%R=10;

berlin=plot(0.0716+0.0391i,'*','DisplayName','Simulink relay Impedances','Color','magenta'); 
plot(0.0984+0.0446i,'*','Color','red');  
plot(0.0772+0.0682i,'*','Color','#A2142F');  
plot(0.0659+0.0538i,'*','Color','black'); 

%% Line to line 
% berlin=plot(0.0115+0.0027i,'*','DisplayName','Simulink relay Impedances','Color','magenta'); 
% plot(-0.0341+0.0375i,'*','Color','red');  
% plot(0.0007+0.0831i,'*','Color','#A2142F');  
% plot(0.0463+0.0483i,'*','Color','black'); 

%R=0.5
% berlin=plot(0.0049+0.0327i,'*','DisplayName','Simulink relay Impedances','Color','magenta'); 
% plot(-0.0063+0.0417i,'*','Color','red');  
% plot(0.0028+0.0529i,'*','Color','#A2142F');  
% plot(0.014+0.0438i,'*','Color','black'); 


%% PLotting relay impedance values for current angles ranging from 0 to 360
% imp2=plot(Apparent_Impedance2,'--','DisplayName','L-L fault relay impedance (Script)','Color','magenta'); 



% % % % plot(imp,'x');
% Res=plot([m*line_impedance1,m*line_impedance1+0.1],'DisplayName','Purely resistive fault error')
% impMV1=plot(Apparent_Impedance2(1),'x','DisplayName','L-L  relay impedance R=1','Color','blue'); 
% impMV2=plot(Apparent_Impedance2(91),'x','DisplayName','L-L  relay impedance R=3','Color','magenta'); 
% impMV3=plot(Apparent_Impedance2(181),'x','DisplayName','L-L  relay impedance R=4','Color','Red'); 
% impMV4=plot(Apparent_Impedance1,'DisplayName','L-G relay impedance R=0.7','Color','green'); 

% impMV1l=plot(temp1l,'--','DisplayName','L-L  relay impedance X-R ratio:5','Color','blue'); 
% impMV2l=plot(temp2l,'--','DisplayName','L-L  relay impedance X-R ratio:10','Color','magenta'); 
% impMV3l=plot(temp3l,'--','DisplayName','L-L  relay impedance X-R ratio:15','Color','red'); 
% impMV4l=plot(temp4l,'--','DisplayName','Line to Line fault relay impedance:R=0.7','Color','green'); 

% 
%% Plots showing relay impedance values at specific current angles
imp0g=plot(Apparent_Impedance1(1),'+','DisplayName','phi=0 (Script)','Color','magenta') 
imp1g=plot(Apparent_Impedance1(91),'+','DisplayName','phi=90 (Script)','Color','red') 
imp2g=plot(Apparent_Impedance1(181),'+','DisplayName','phi=180(Script)','Color','#A2142F') 
imp3g=plot(Apparent_Impedance1(271),'+','DisplayName','phi=270 (Script)','Color','black') 
 
% imp0l=plot(Apparent_Impedance2(1),'+','DisplayName','phi=0 (Script)','Color','magenta') 
% imp1l=plot(Apparent_Impedance2(91),'+','DisplayName','phi=90 (Script)','Color','red') 
% imp2l=plot(Apparent_Impedance2(181),'+','DisplayName','phi=180(Script)','Color','#A2142F') 
% imp3l=plot(Apparent_Impedance2(271),'+','DisplayName','phi=270 (Script)','Color','black') 
% legend([true imp0g imp1g imp2g imp3g imp3l impMV1 impMV2 impMV3  impMV1l impMV2l impMV3l ]); 
% R=2
% error1=abs(0.0174+0.0417i-Apparent_Impedance1(1))
% error2=abs(0.0219+0.0444i-Apparent_Impedance1(91)) 
% error3=abs(0.0172+0.047i-Apparent_Impedance1(181)) 
% error4=abs(0.0158+0.0442i-Apparent_Impedance1(271)) 
%R=10
%% Absolute magnitude differences between Simulink and Matlab scripts
% error1=abs(0.0716+0.0391i-Apparent_Impedance1(1))
% error2=abs(0.0984+0.0446i-Apparent_Impedance1(91)) 
% error3=abs(0.0772+0.0682i-Apparent_Impedance1(181)) 
% error4=abs(0.0659+0.0538i-Apparent_Impedance1(271)) 
% imp0l imp1l imp2l imp3l
% impMV1 impMV2 impMV3
legend([ imp1 berlin imp0g imp1g imp2g imp3g])
hold off



