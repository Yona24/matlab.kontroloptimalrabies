function dk=ST(y,Ad,alphad,Cd,md,miud,upsd,rhod,bhetadd,bhetadh,deltahepsh,deltaepsd,deltagammad,Bh,alphah,mh,miuh,upsh,rhoh,deltahgammah)
x(1)=Ad-((1-upsd)*bhetadd*y(1)*y(3))-((md+upsd)*y(1))+(deltaepsd*y(2)+alphad*y(4));
x(2)=((1-upsd)*bhetadd*y(1)*y(3))-((1-rhod)*deltagammad+md+rhod+deltaepsd+Cd)*y(2);
x(3)=((1-rhod)*deltagammad*y(2))-((md+miud)*y(3));
x(4)=(upsd*y(1))+(rhod*y(2))-((md+alphad)*y(4));
x(5)=Bh-((1-upsh)*bhetadh*y(5)*y(3))-((mh+upsh)*y(5))+(deltahepsh*y(6)+alphah*y(8));
x(6)=((1-upsh)*bhetadh*y(5)*y(3))-((1-rhoh)*deltahgammah+mh+rhoh+deltahepsh)*y(6);
x(7)=((1-rhoh)*deltahgammah*y(6))-(mh+miuh)*y(7);
x(8)=(upsh*y(5))+(rhoh*y(6))-((mh+alphah)*y(8));
dk=[x(1) x(2) x(3) x(4) x(5) x(6) x(7) x(8)];
