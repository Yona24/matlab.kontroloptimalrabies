function dd=CS(y,alphad,Cd,md,miud,upsd,rhod,bhetadd,bhetadh,deltahepsh,deltaepsd,deltagammad,alphah,mh,miuh,upsh,rhoh,deltahgammah,SD,ID,SH)
d1=y(1); d2=y(2); d3=y(3); d4=y(4);d5=y(5); d6=y(6); d7=y(7); d8=y(8);
x(1)=(d1*(((1-upsd)*bhetadd*ID)+md+upsd))-(d2*((1-upsd)*bhetadd*ID))-(d4*upsd);
x(2)=d2*((1-rhod)*deltagammad+md+rhod+deltaepsd+Cd)-(d1*deltaepsd)-(d3*(1-rhod)*deltagammad)-(d4*rhod)-1;
x(3)=(d3*(md+miud))+(d1*(1-upsd)*bhetadd*SD)+(d5*(1-upsh)*bhetadh*SH)-(d2*(1-upsd)*bhetadd*SD)-(d6*(1-upsh)*bhetadh*SH)-1;
x(4)=(d4*(md+alphad)-(d1*alphad));
x(5)=(d5*(((1-upsh)*bhetadh*ID)+mh+upsh))-(d6*(1-upsh)*bhetadh*ID)-(d8*upsh);
x(6)=(d6*((1-rhoh)*deltahgammah+mh+rhoh+deltahepsh))-(d5*deltahepsh)-(d7*(1-rhoh)*deltahgammah)-(d8*rhoh)-1;
x(7)=(d7*(mh+miuh))-1;
x(8)=(d8*(mh+alphah))-(d5*alphah);

dd=[x(1) x(2) x(3) x(4) x(5) x(6) x(7) x(8)];


