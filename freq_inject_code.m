%% SImulink outputs 
% v_simu=out.relay_voltage;
% seq_simu=out.relay_voltage_a;
% i_simu=out.relay_current_a; 
% i_grid=out.grid_current;
% imp=(v_simu(1)-v_simu(2))/i_simu(1) 
% 
% f_curr=i_simu-i_grid;

%% Base parameters 
Vb = 315e3;
Sb = 300e6;
Zb = Vb^2/Sb;
Ib = Sb/(sqrt(3)*Vb);
n=2;
freq=n*60;
w = 2*pi*freq; 
%% Fault Params
% pos=[0.25  0.5  0.75  0.95]
%   for p=1:4
m=0.75; 
Rf = linspace(0,10,100);
Fault_Type = 1; % LG = 1, LL = 2; 


     % Line Parameters Model
    leng = 50; %km
    R10 = [ 0.02715    0.25155];
    L10 = [0.00100    0.00364];
    C10 = [ 1.16622e-08  6.99588e-09];
%% Line impedance at 120hz
%     ZL1 = (2.8323 + 1i*77.1077)*leng/100; 
%     ZL2 = ZL1;
%     ZL0 = (27.651+288.38i)*leng/100;

%% Line impedance at 240hz
ZL1=(3.251+165.7i)*leng/100; 
ZL2 = ZL1; 
ZL0=(38.95+686.43i)*leng/100;
    

ZL = [ZL1; ZL2; ZL0];
    K0 = (ZL0/ZL1-1);
    

% Grid Impedance
Psc = 5000e6; % Short Circuit level
ZG_abs = Vb^2/Psc;
XR_ratio = 10;
Rs = ZG_abs/sqrt(1+XR_ratio^2);
Ls = (Rs*XR_ratio)/(2*pi*60);
Zs = Rs + 1i*Ls*w; 

% Transformer Winding Impedance
Rt_pu = 0;
Lt_pu = 0.04;

Lt=0.04*Zb/(2*pi*60);
Zt = 2*(Rt_pu + 1i*w*Lt); % 2 x because of two windings

% Mho Setting
reach = 0.8;


%phase
limit_LG=0; 
limit_LL=0;
V2=1;

z=1;
PhaseA=-30 + z-1; 
PhaseB=210 + z-1; 
PhaseC=90 + z-1;

%% Theory
for phi=1:100 
Rfault = Rf(phi)/Zb;
Transformer_imp=Zt/Zb; 
line_impedance1=ZL1/Zb; 
line_impedance0=ZL0/Zb; 
S2=Zs/Zb;
curr = 0.25*exp(1i*0*pi/180);
%% Line to ground fault
if Fault_Type==1
V_th=(curr*((1-m)*line_impedance1+S2)); 
Z_1_thev=(S2+(1-m)*line_impedance1);
Z_0_thev=(1/(1/(S2+(1-m)*line_impedance0)+1/(Transformer_imp+m*line_impedance0))); 
Z_2_thev=Z_1_thev; 

I_a1=V_th/(Z_1_thev+Z_2_thev+Z_0_thev+3*Rfault); 
I_fault=3*I_a1;
I_r1=curr; 
I_r2=0; 
I_r0=I_a1*(S2+(1-m)*line_impedance0)/(Transformer_imp+S2+line_impedance0); 
I_AR=I_r1+I_r2+I_r0;

V_r1=(I_r1*m*line_impedance1)+(V_th-I_a1*Z_1_thev); 

V_r2=-I_a1*Z_2_thev;

V_r0=I_r0*m*line_impedance0-I_a1*Z_0_thev;


V_AR=V_r1+V_r2+V_r0; 

Apparent_Impedance1(phi)=V_AR/(I_AR+K0*I_r0); 
true_imp=v_simu(1)/(i_simu(1)+K0*seq_simu(3));
%% Code to identify fault resistance where classical mho is exceeded
% if angle(reach*line_impedance1-Apparent_Impedance1(phi))-angle(Apparent_Impedance1(phi))>=(pi/2) && limit_LG==0
% limit_LG=Apparent_Impedance1(phi) 
% limit_Rf=Rfault*Zb
% end
%% If statement to allow plotting of multiple relay impedances with varying fault resistances for varying network parameters
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


end
%% Line to line fault 
if Fault_Type==2
Z_1_thev=(S2+(1-m)*line_impedance1);
Z_2_thev=Z_1_thev; 

V_th=(curr*((1-m)*line_impedance1+S2)); 
 
I_a1=V_th/(Z_1_thev+Z_2_thev+Rfault);

I_fault=-i*sqrt(3)*I_a1; 

V_r1=(curr*m*line_impedance1+(V_th-(I_a1*Z_1_thev))); 

V_r2=(I_a1*(Rfault+Z_2_thev)*(Z_2_thev/(Z_2_thev+Rfault))); 

V_AR=V_r1+V_r2; 

I_AR=curr;

Apparent_Impedance2(phi)=(V_r1-V_r2)/curr; 

if angle(reach*line_impedance1-Apparent_Impedance2(phi))-angle(Apparent_Impedance2(phi))>=(pi/2) && limit_LL==0
limit_LL=Apparent_Impedance1(phi) 
limit_Rf_LL=Rfault*Zb
end
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

end 
%% Plotting circular mho 
theta=-2*pi:0.00001:2*pi;
x=(reach/2)*(abs(line_impedance1))*(exp(i*(theta))+exp(i*angle(line_impedance1)));

%% Plotting quadrilateral mho
plot([0, i*reach*abs(line_impedance1)-0.02,i*reach*abs(line_impedance1)+abs(line_impedance1),abs(line_impedance1)-0.02,0],'black'); 
     
hold on
plot(x,'black') 
test=[0,line_impedance1];
plot(test,'--'); 
true=plot(m*line_impedance1,'X','DisplayName','Position of fault along line'); 
%      end
imp1=plot(Apparent_Impedance1,'DisplayName','L-G fault relay impedance'); 

%% Plotting values given by Simulink for verification purposes
%% lg 
% plot(0.00837094111720902 + 0.0852783130100865i,'x'); 
% plot(0.0137260077887392 + 0.0859687544437430i,'*');  
% plot(0.0189673292199148 + 0.0870209207994856i,'.');  
% plot(0.0238595469764596 + 0.0883530830634956i,'+'); 
%% ll 
% plot(0.0069 + 0.0858i,'x'); 
% plot(0.0107635752848137 + 0.0861359682261649i,'*');  
% plot(0.0145480425008113 + 0.0866122688135489i,'.');  
% plot(0.0181326851517763 + 0.0872423464053490i,'+'); 
%% 240hz 
%lg 
% plot(0.00780510096069457 + 0.165433349700472i,'x'); 
% plot(0.0132885301281304 + 0.165768066617194i,'*');  
% plot(0.0187458377393607 + 0.166276774870447i,'.');  
% plot(0.0239523038837554 + 0.166926815124635i,'+');



% impMV1=plot(temp1,'DisplayName','L-G fault relay impedance 25%','Color','blue'); 
% impMV2=plot(temp2,'DisplayName','L-G fault relay impedance 50%','Color','magenta'); 
% impMV3=plot(temp3,'DisplayName','L-G fault relay impedance 75%','Color','Red'); 
% impMV4=plot(temp4,'DisplayName','L-G fault relay impedance 95%','Color','Green'); 
% 



% imp2=plot(Apparent_Impedance2,'--','DisplayName','Line to Line fault relay impedance','Color','BLue'); 
% impMV1l=plot(temp1l,'--','DisplayName','L-L fault relay impedance 25%','Color','blue'); 
% impMV2l=plot(temp2l,'--','DisplayName','L-L fault relay impedance 50%','Color','magenta'); 
% impMV3l=plot(temp3l,'--','DisplayName','L-L fault relay impedance 75%','color','Red'); 
% impMV4l=plot(temp4l,'--','DisplayName','L-L fault relay impedance 95%','Color','Green');
% % plot(imp,'X'); 
impMV1x=plot(Apparent_Impedance1(26),'x','DisplayName','L-L fault relay impedance R=2.5','Color','black'); 
% plot(temp2(20),'x','DisplayName','L-G fault relay impedance R=2','Color','black'); 
impMV1y=plot(Apparent_Impedance1(51),'*','DisplayName','L-L fault relay impedance R=5'); 
impMV1z=plot(Apparent_Impedance1(76),'.','DisplayName','L-L fault relay impedance: R=7.5'); 
impMV1q=plot(Apparent_Impedance1(100),'+','DisplayName','L-L fault relay impedance: R=10'); 



legend([ imp1 true impMV1x impMV1y impMV1z impMV1q]) 
hold off 




