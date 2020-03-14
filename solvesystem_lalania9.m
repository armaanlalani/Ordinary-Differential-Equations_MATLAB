function [t,y]= solvesystem_lalania9(f1,f2,t0,tN,x0,h)
 y1=[x0(1)];
 y2=[x0(2)];
 t=[t0];
 j=1;
 for i=t0:h:tN
    m1=f1(i,y1(j),y2(j));
    m3=f2(i,y1(j),y2(j));
    m2=f1(i+h,y1(j)+h*m1,y2(j)+h*m3);
    m4=f2(i+h,y1(j)+h*m1,y2(j)+h*m3);
    y1(j+1)=y1(j) + h*(m1+m2)/2;
    y2(j+1)=y2(j) + h*(m3+m4)/2;
    t(j+1)=t(j)+h;
    j=j+1;
 end
 y=[y1;y2];