clear all;
clc;
tic
SD0=6*10^(7);  ED0=3*10^(6); ID0=1.5*10^(6); RD0=0; 
SH0=1.29*10^(9); EH0=8*10^(6);IH0=8.9*10^(6); RH0=0;
Ad=3*10^(6); alphad=1; Cd=0.3; md=0.0056; miud=0.3; bhetadd=1.58*10^(-7);
 deltaepsd=3.6; Bh=0.314; bhetadh=2.29*10^(-8);alphah=1; mh=0.0074; miuh=0.7;
 deltahgammah=2.4; deltagammad=2.4; deltahepsh=3.6;
 B1=100; B2=400; B3=100; B4=400;
upsdmax=1; rhodmax=1; upshmax=1; rhohmax=1;
h=0.1;
t=0:h:10;
N=length(t);
SD=zeros(N,1); SDo=zeros(N,1);
ED=zeros(N,1); EDo=zeros(N,1);
ID=zeros(N,1); IDo=zeros(N,1);
RD=zeros(N,1); RDo=zeros(N,1);
SH=zeros(N,1); SHo=zeros(N,1);
EH=zeros(N,1); EHo=zeros(N,1);
IH=zeros(N,1); IHo=zeros(N,1);
RH=zeros(N,1); RHo=zeros(N,1);

d1=zeros(N,1); d1o=zeros(N,1);
d2=zeros(N,1); d2o=zeros(N,1);
d3=zeros(N,1); d3o=zeros(N,1);
d4=zeros(N,1); d4o=zeros(N,1);
d5=zeros(N,1); d5o=zeros(N,1);
d6=zeros(N,1); d6o=zeros(N,1);
d7=zeros(N,1); d7o=zeros(N,1);
d8=zeros(N,1); d8o=zeros(N,1);

upsd=zeros(N,1); upsdo=zeros(N,1);
rhod=zeros(N,1); rhodo=zeros(N,1);
upsh=zeros(N,1); upsho=zeros(N,1);
rhoh=zeros(N,1); rhoho=zeros(N,1);
tes=1;
it=0;
while tes>1e-5
    upsdo=upsd; rhodo=rhod; upsho=upsh; rhoho=rhoh;
    SDo=SD; EDo=ED; IDo=ID; RDo=RD; SHo=SH; EHo=EH; IHo=IH; RHo=RH;
    d1o=d1; d2o=d2; d3o=d3; d4o=d4; d5o=d5; d6o=d6; d7o=d7; d8o=d8;
    SD(1)=SD0; ED(1)=ED0; ID(1)=ID0; RD(1)=RD0;
    SH(1)=SH0; EH(1)=EH0; IH(1)=IH0; RH(1)=RH0;
    J(it+1)=0;
   for i=1:N-1
  y=[SD(i)  ED(i)  ID(i)  RD(i)   SH(i)  EH(i) IH(i) RH(i)];
        k1=h*ST(y,Ad,alphad,Cd,md,miud,upsdo(i),rhodo(i),bhetadd,bhetadh,deltahepsh,deltaepsd,deltagammad,Bh,alphah,mh,miuh,upsho(i),rhoho(i),deltahgammah);
        k2=h*ST(y+0.5*k1,Ad,alphad,Cd,md,miud,upsdo(i),rhodo(i),bhetadd,bhetadh,deltahepsh,deltaepsd,deltagammad,Bh,alphah,mh,miuh,upsho(i),rhoho(i),deltahgammah);
        k3=h*ST(y+0.5*k2,Ad,alphad,Cd,md,miud,upsdo(i),rhodo(i),bhetadd,bhetadh,deltahepsh,deltaepsd,deltagammad,Bh,alphah,mh,miuh,upsho(i),rhoho(i),deltahgammah);
        k4=h*ST(y+k3,Ad,alphad,Cd,md,miud,upsdo(i),rhodo(i),bhetadd,bhetadh,deltahepsh,deltaepsd,deltagammad,Bh,alphah,mh,miuh,upsho(i),rhoho(i),deltahgammah);
        y=y+(1/6)*(k1+2*k2+2*k3+k4);
        
        SD(i+1)=y(1);
        ED(i+1)=y(2);
        ID(i+1)=y(3);
        RD(i+1)=y(4);
        SH(i+1)=y(5);
        EH(i+1)=y(6);
        IH(i+1)=y(7);
        RH(i+1)=y(8);
    J(i+1)=J(it+1)+h*((ED(i))+(EH(i))+(ID(i))+(IH(i))+((B1/2)*(upsdo(i)^2))+((B2/2)*(rhodo(i)^2))+((B3/2)*(upsho(i)^2))+((B4/2)*(rhoho(i)^2)));
   end
  if it==0
        figure(1)
        subplot(2,1,1)
        plot(t,SD,'y','linewidth',1.5);
        hold on;
        subplot(2,1,2)
        plot(t,SH,'m','linewidth',1.5);
        hold on;
        figure(2)
        subplot(2,1,1)
        plot(t,ED,'y','linewidth',1.5);
        hold on;
        subplot(2,1,2)
        plot(t,EH,'m','linewidth',1.5);
        hold on; 
        figure(3)
        subplot(2,1,1)
        plot(t,ID,'y','linewidth',1.5);
        hold on;
        subplot(2,1,2)
        plot(t,IH,'m','linewidth',1.5);
        hold on;
        figure(4)
        subplot(2,1,1)
        plot(t,RD,'y','linewidth',1.5); 
        hold on;
        subplot(2,1,2)
        plot(t,RH,'m','linewidth',1.5); 
        hold on;
  end
   d1(N)=0; d2(N)=0; d3(N)=0; d4(N)=0; d5(N)=0; d6(N)=0; d7(N)=0; d8(N)=0;
   upsd(N)=0;rhod(N)=0;upsh(N)=0; rhoh(N)=0;    
        
       for i=1:N-1
        je=N-i;  
        y=[d1(je+1) d2(je+1) d3(je+1) d4(je+1)  d5(je+1) d6(je+1) d7(je+1) d8(je+1)];
        
        k1=h*CS(y,alphad,Cd,md,miud,upsdo(je+1),rhodo(je+1),bhetadd,bhetadh,deltahepsh,deltaepsd,deltagammad,alphah,mh,miuh,upsho(je+1),rhoho(je+1),deltahgammah,SD(je+1),ID(je+1),SH(je+1));
        k2=h*CS(y+0.5*k1,alphad,Cd,md,miud,upsdo(je+1),rhodo(je+1),bhetadd,bhetadh,deltahepsh,deltaepsd,deltagammad,alphah,mh,miuh,upsho(je+1),rhoho(je+1),deltahgammah,SD(je+1),ID(je+1),SH(je+1));
        k3=h*CS(y+0.5*k2,alphad,Cd,md,miud,upsdo(je+1),rhodo(je+1),bhetadd,bhetadh,deltahepsh,deltaepsd,deltagammad,alphah,mh,miuh,upsho(je+1),rhoho(je+1),deltahgammah,SD(je+1),ID(je+1),SH(je+1));
        k4=h*CS(y+k3,alphad,Cd,md,miud,upsdo(je+1),rhodo(je+1),bhetadd,bhetadh,deltahepsh,deltaepsd,deltagammad,alphah,mh,miuh,upsho(je+1),rhoho(je+1),deltahgammah,SD(je+1),ID(je+1),SH(je+1));
        y=y-(1/6)*(k1+2*k2+2*k3+k4);
        
        d1(je)=y(1); d2(je)=y(2); d3(je)=y(3); d4(je)=y(4); d5(je)=y(5); d6(je)=y(6); d7(je)=y(7); d8(je)=y(8);
        
        Banding1=(((d1(je)-d4(je))*SD(je)+((d2(je)-d1(je))*bhetadd*ID(je)*SD(je)))/B1);
        Banding2=(((d2(je)-d4(je))*ED(je))+((d3(je)-d2(je))*deltagammad*ED(je))/B2);
        Banding3=(((d5(je)-d8(je))*SH(je))+((d6(je)-d5(je))*bhetadh*SH(je)*ID(je))/B3);
        Banding4=(((d6(je)-d8(je))*EH(je)+(d7(je)-d6(je))*deltahgammah*EH(je))/B4);
        temp1=max([0 Banding1]);
        temp2=max([0 Banding2]);
        temp3=max([0 Banding3]);
        temp4=max([0 Banding4]);
        upsd(je)=min([temp1 upsdmax]);
        rhod(je)=min([temp2 rhodmax]);
        upsh(je)=min([temp3 upshmax]);
        rhoh(je)=min([temp4 rhohmax]);
       end
       
    eSD=sum(abs(SD-SDo));
    eED=sum(abs(ED-EDo));
    eID=sum(abs(ID-IDo));
    eRD=sum(abs(RD-RDo));
    eSH=sum(abs(SH-SHo));
    eEH=sum(abs(EH-EHo));
    eIH=sum(abs(IH-IHo));
    eRH=sum(abs(RH-RHo));
    ed1=sum(abs(d1-d1o));
    ed2=sum(abs(d2-d2o));
    ed3=sum(abs(d3-d3o));
    ed4=sum(abs(d4-d4o));
    ed5=sum(abs(d5-d5o));
    ed6=sum(abs(d6-d6o));
    ed7=sum(abs(d7-d7o));
    ed8=sum(abs(d8-d8o));
    eupsd=sum(abs(upsd-upsdo));
    erhod=sum(abs(rhod-rhodo));
    eupsh=sum(abs(upsh-upsho));
    erhoh=sum(abs(rhoh-rhoho));
    tes=eSD+eED+eID+eRD+eSH+eEH+eIH+eRH+ed1+ed2+ed3+ed4+ed5+ed6+ed7+ed8+eupsd+erhod+eupsh+erhoh;
   it=it+1;
   upsd=(0.5*upsd+0.5*upsdo);
   rhod=(0.5*rhod+0.5*rhodo);
   upsh=(0.5*upsh+0.5*upsho);
   rhoh=(0.5*rhoh+0.5*rhoho);
end
%Plot populasi dengan kontrol
figure(1)
subplot(2,1,1)
plot(t,SD,'r--','LineWidth',1.5);
xlabel('Time (years)');
ylabel('SD(t)');
grid on;
legend('SD without control','SD with control');
title('Susceptible dogs');
hold on;
subplot(2,1,2)
plot(t,SH,'b--','LineWidth',1.5);
xlabel('Time (years)');
ylabel('SH(t)');
grid on;
legend('SH without control','SH with control');
title('Susceptible humans');
hold on;
figure(2)
subplot(2,1,1)
plot(t,ED,'r--','LineWidth',1.5);
xlabel('Time (years)');
ylabel('ED(t)');
grid on;
legend('ED  without control','ED with control');
title('Exposed dogs');
hold on;
subplot(2,1,2)
plot(t,EH,'b--','LineWidth',1.5);
xlabel('Time (years)');
ylabel('EH(t)');
grid on;
legend('EH without control','EH with control');
title('Exposed humans');
hold on;
figure(3)
subplot(2,1,1)
plot(t,ID,'r--','LineWidth',1.5);
xlabel('Time (years)');
ylabel('ID(t)');
grid on;
legend('ID without control','ID with control');
title('Infected dogs');
hold on;
subplot(2,1,2)
plot(t,IH,'b--','LineWidth',1.5);
xlabel('Time (years)');
ylabel('IH(t)');
grid on;
legend('IH without control','IH with control');
title('Infected Humans');
hold on;
figure(4)
subplot(2,1,1)
plot(t,RD,'r--','LineWidth',1.5);
xlabel('Time (years)');
ylabel('RD(t)');
grid on;
legend('RD without control','RD with control');
title('Recovered dogs');
hold on;
subplot(2,1,2)
plot(t,RH,'b--','LineWidth',1.5);
xlabel('Time (years)');
ylabel('RH(t)');
grid on;
legend('RH without control','RH with control');
title('Recovered humans');
hold on;
 figure(5)
         plot(t,upsd,'m--','LineWidth',2);
         hold on;
       plot(t,rhod,'r--','LineWidth',2);
       hold on;
        plot(t,upsh,'b--','LineWidth',2);
        hold on;
       plot(t,rhoh,'g--','LineWidth',2);
       hold on;
         xlabel('Time (years)');
       legend('u1','u2','u3','u4');
       grid on;
        figure(6)
        subplot(2,2,1)
        plot(t,upsd,'m--','LineWidth',2);
        xlabel('Waktu (Tahun)');
        legend('u1');
        grid on;
        hold on;
        subplot(2,2,2)
        plot(t,rhod,'r--','LineWidth',2);
        xlabel('Waktu (Tahun)');
        legend('u2');
        grid on;
        hold on;
        subplot(2,2,3)
        plot(t,upsh,'b--','LineWidth',2);
        xlabel('Waktu (Tahun)');
        legend('u3');
        grid on;
        hold on;
        subplot(2,2,4)
        plot(t,rhoh,'g--','LineWidth',2);
        xlabel('Waktu (Tahun)');
        legend('u4');
        grid on;
        hold on;
        figure(7)
        plot(t,J,'m','LineWidth',2);
        xlabel('Waktu (Tahun)');
        legend('J(u1,u3)','J(u2,u4)','J(u1,u2,u3,u4)');
        grid on;
        hold on;    
toc