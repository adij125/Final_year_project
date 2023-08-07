
purple=[107 76 154]./255;
%% SImulink outputs 
% v_simu=out.relay_voltage;
% i_simu=out.relay_current_a; 
% seq_simu=out.relay_voltage_a;
% true_imp_ll=(v_simu(1)-v_simu(2))/(seq_simu(1)-seq_simu(2))
% true_imp_lg=v_simu(1)/(i_simu(1)+K0*seq_simu(3))
% fault=i_simu-curr_grid
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
%     R10 = [0.03441    0.29318];
%     L10 = [0.00148    0.00374];
%     C10 = [7.84869e-09 5.67336e-09]; 
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
Lt_pu = 0.04;
Zt = 2*(1i*Lt_pu)*Zb; % 2 x because of two windings

% Mho Setting
reach = 0.8;

% Source Impedance 

Lg=0.05*Zb/(w); 
Zg=(1i*w*Lg);


%% Theory
limit_LG=0; 
limit_LL=0;

V2 = 1;
Rfault = Rf/Zb;
Transformer_imp=Zt/Zb; 
line_impedance1=ZL1/Zb; 
line_impedance0=ZL0/Zb; 
S2=Zs/Zb;
S1=Zg/Zb;




%phase

for phi=1:361 
% phi=131;
current_injection=1*exp(1i*(phi-1)*pi/180);
V1=current_injection*(Transformer_imp + S1+S2 + line_impedance1)+V2; 

%% Magnitude and angle for Simulink
% V1_test=1*exp(1i*(271-1)*pi/180)*(Transformer_imp + S1+S2 + line_impedance1)+V2; 
% V1_simu=abs(V1_test); 
% PhaseA=-30+(angle(V1_test)*180/pi);
% PhaseB=210+(angle(V1_test)*180/pi); 
% PhaseC=90+(angle(V1_test)*180/pi);

V_th=((V1-V2)/(Transformer_imp + S1+S2 + line_impedance1))*((1-m)*line_impedance1+S2)+V2;
 
Z_1_thev=1/(1/(Transformer_imp+S1+m*line_impedance1)+1/(S2+(1-m)*line_impedance1)); 
Z_2_thev=Z_1_thev; 
Z_0_thev=1/(1/(Transformer_imp+m*line_impedance0)+1/(S2+(1-m)*line_impedance0)); 

%% line to ground fault
if Fault_Type==1
I_fault=V_th/(Z_1_thev+Z_2_thev+Z_0_thev+3*Rfault); 


I_r1=I_fault*(((1-m)*line_impedance1+S2)/(S1+S2+Transformer_imp+line_impedance1))+((V1-V2)/(Transformer_imp + S2 +S1+ line_impedance1)); 
I_r2=I_fault*(((1-m)*line_impedance1+S2)/(S1+S2+Transformer_imp+line_impedance1)); 
I_r0=I_fault*(((1-m)*line_impedance0+S2)/(S2+Transformer_imp+line_impedance0));

V_r1=I_r1*(m*line_impedance1) + (V_th-I_fault*(Z_1_thev)); 
V_r2=I_r2*(m*line_impedance1) -I_fault*(Z_2_thev);
V_r0= -I_fault*(Z_0_thev)+I_r0*(m*line_impedance0);


I_AR=I_r1+I_r2+I_r0;
V_AR=V_r1+V_r2+V_r0;

Apparent_Impedance1(phi)=V_AR/(I_AR+K0*I_r0); 

%% if statement to allow plotting of multiple outputs
% if p==1 
%     temp1(phi)=Apparent_Impedance1(phi); 
% elseif p==2 
%     temp2(phi)=Apparent_Impedance1(phi); 
% elseif p==3 
%     temp3(phi)=Apparent_Impedance1(phi); 
% elseif p==4 
%     temp4(phi)=Apparent_Impedance1(phi);
% 
% end 

%% Code that allows detection of fault resistance where relay impedance falls outside mho
% error=(abs(Apparent_Impedance1)-abs(true_imp))*100/(abs(true_imp))
% if angle(reach*line_impedance1-Apparent_Impedance1(phi))-angle(Apparent_Impedance1(phi))>=(pi/2) && limit_LG==0
% limit_LG=Apparent_Impedance1(phi) 
% limit_Rf=Rfault*Zb
% end

end


%% line to line fault 
if Fault_Type==2
I_a1=V_th/(Z_1_thev+Z_2_thev+Rfault);
I_r1=I_a1*(((1-m)*line_impedance1+S2)/(S1+S2+Transformer_imp+line_impedance1))+((V1-V2)/(Transformer_imp + S1+S2 + line_impedance1)); 
I_r2=-I_a1*(((1-m)*line_impedance1+S2)/(S1+S2+Transformer_imp+line_impedance1)); 

I_fault=-i*sqrt(3)*I_a1;
I_AR=I_r1+I_r2; 

V_r1=I_r1*m*line_impedance1+I_a1*(Rfault+Z_2_thev); 

V_r2=I_r2*m*line_impedance1+(I_a1*(Rfault+Z_2_thev)-I_a1*Rfault);

V_AR=V_r1+V_r2;

Apparent_Impedance2(phi)=(V_r1-V_r2)/(I_r1-I_r2); 

% if p==1 
%     temp1l(phi)=Apparent_Impedance2(phi); 
% elseif p==2 
%     temp2l(phi)=Apparent_Impedance2(phi); 
% elseif p==3 
%     temp3l(phi)=Apparent_Impedance2(phi); 
% elseif p==4 
%     temp4l(phi)=Apparent_Impedance2(phi);
% 
% end 
end

% if angle(reach*line_impedance1-Apparent_Impedance2(phi))-angle(Apparent_Impedance2(phi))>=(pi/2) && limit_LL==0
% limit_LL=Apparent_Impedance2(phi) 
% R_limit=Rfault*Zb
% end
% end


end
hold on
true=plot(m*line_impedance1,'x','DisplayName','position of fault along line impedance','Color','black');
    
%% Plotting classical mho
theta=-2*pi:0.00001:2*pi;
x=(reach/2)*(abs(line_impedance1))*(exp(1i*(theta))+exp(1i*angle(line_impedance1))); 
% y=(reach/4)*(abs(line_impedance1))*(exp(1i*(theta))+exp(1i*angle(line_impedance1)));
%% Plotting the apparent impedances where 1 is the line to ground impedance and 2 is the line to line impedance


plot(x,'Color','black') 


%% Plotting quadrilateral mho
plot([0, i*reach*abs(line_impedance1)-0.02,i*reach*abs(line_impedance1)+abs(line_impedance1),abs(line_impedance1)-0.02,0],'black');
     

test=[0,line_impedance1];
plot(test,'--')
     
% imp1=plot(Apparent_Impedance1,'DisplayName','L-G fault relay impedance (Script)','Color','blue'); 

%% Relay impedance values obtained from SImulink for verification purposes
%% R=2 l-g
% berlin=plot(0.0129+0.0426i,'*','DisplayName','Simulink relay Impedances','Color','magenta'); 
% plot(0.0143+0.0435i,'*','Color','red');  
% plot(0.0129+0.0446i,'*','Color','#A2142F');  
% plot(0.0122+0.0436i,'*','Color','black'); 
%% R=10 l-g 
% berlin=plot(0.0506+0.0428i,'*','DisplayName','Simulink relay Impedances','Color','magenta'); 
% plot(0.0598+0.0446i,'*','Color','red');  
% plot(0.0542+0.0536i,'*','Color','#A2142F');  
% plot(0.0485+0.049i,'*','Color','black'); 
%% R=2 l-l
% berlin=plot(0.0127+0.0417i,'*','DisplayName','Simulink relay Impedances','Color','magenta'); 
% plot(0.0155+0.0432i,'*','Color','red');  
% plot(0.0127+0.045i,'*','Color','#A2142F');  
% plot(0.0117+0.0433i,'*','Color','black');  

%R=10 
berlin=plot(0.0495+0.0379i,'*','DisplayName','Simulink relay Impedances','Color','magenta'); 
plot(0.065+0.0415i,'*','Color','red');  
plot(0.0532+0.0544i,'*','Color','#A2142F');  
plot(0.0461+0.0468i,'*','Color','black'); 

%% Code to show relay impedance plots at all injected angles for different network parameters
imp2=plot(Apparent_Impedance2,'--','DisplayName','L-L fault relay impedance (Script)','Color','blue');
% impMV1=plot(temp1,'DisplayName','L-G fault relay impedance:25%','Color','blue'); 
% impMV2=plot(temp2,'DisplayName','L-G fault relay impedance:50%','Color','magenta'); 
% impMV3=plot(temp3,'DisplayName','L-G fault relay impedance:75%','Color','red'); 
% impMV4=plot(temp4,'DisplayName','L-G fault relay impedance:95%','Color','green'); 
% % 
% impMV1l=plot(temp1l,'--','DisplayName','L-L fault relay impedance:25%','Color','blue'); 
% impMV2l=plot(temp2l,'--','DisplayName','L-L fault relay impedance:50%','Color','magenta'); 
% impMV3l=plot(temp3l,'--','DisplayName','L-L fault relay impedance:75%','Color','red'); 
% impMV4l=plot(temp4l,'--','DisplayName','L-L fault relay impedance:95%','Color','green'); 
% plot(true_imp,'x'); 

%%Code to plot relay impedance values at a specific current angle
% imp0g=plot(Apparent_Impedance1(1),'+','DisplayName','phi=0 (Script)','Color','magenta'); 
% imp1g=plot(Apparent_Impedance1(91),'+','DisplayName','phi=90 (Script)','Color','red'); 
% imp2g=plot(Apparent_Impedance1(181),'+','DisplayName','phi=180 (Script)','Color','#A2142F'); 
% imp3g=plot(Apparent_Impedance1(271),'+','DisplayName','phi=270 (Script)','Color','black'); 
% % plot(imp,'x'); 
imp0l=plot(Apparent_Impedance2(1),'+','DisplayName','phi=0 (Script)','Color','magenta'); 
imp1l=plot(Apparent_Impedance2(91),'+','DisplayName','phi=90 (Script)','Color','red'); 
imp2l=plot(Apparent_Impedance2(181),'+','DisplayName','phi=180(Script)','Color','#A2142F'); 
imp3l=plot(Apparent_Impedance2(271),'+','DisplayName','phi=270 (Script)','Color','black'); 
% legend([impMV1 impMV2 impMV3 impMV4 true impMV1l impMV2l impMV3l impMV4l imp0g imp1g imp2g imp3g]) 

%%Code to calculate absolute magnitude differences between SImulink and script relay impedances

%%lg R=2
% error1=abs(0.0129+0.0426i-Apparent_Impedance1(1))
% error2=abs(0.0143+0.0435i-Apparent_Impedance1(91))
% error3=abs(0.0129+0.0446i-Apparent_Impedance1(181))
% error4=abs(0.0122+0.0436i-Apparent_Impedance1(271)) 
%%l R=10
% error1=abs(0.0506+0.0428i-Apparent_Impedance1(1))
% error2=abs(0.0598+0.0446i-Apparent_Impedance1(91))
% error3=abs(0.0542+0.0536i-Apparent_Impedance1(181))
% error4=abs(0.0485+0.049i-Apparent_Impedance1(271))

legend([imp2 berlin imp0l imp1l imp2l imp3l])
hold off

 


