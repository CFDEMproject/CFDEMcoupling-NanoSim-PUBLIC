% postproc.m: loads the data necessary and saves the graphs 

clear;
clc;
close all;

x=linspace(1.5,6,10);
val=linspace(1,4.5,10);
sol=linspace(1,4.5,10);


load ../c3po_dataStorage_particles/time_0_0.75_processor0.h5 ;
val(1)= U_Favre_0(1,1);
close all ;

load ../c3po_dataStorage_particles/time_0_1.0_processor0.h5 ;
val(2)= U_Favre_0(1,1);
close all ;

load ../c3po_dataStorage_particles/time_0_1.25_processor0.h5 ;
val(3)= U_Favre_0(1,1);
close all ;

load ../c3po_dataStorage_particles/time_0_1.5_processor0.h5 ;
val(4)= U_Favre_0(1,1);
close all ;

load ../c3po_dataStorage_particles/time_0_1.75_processor0.h5 ;
val(5)= U_Favre_0(1,1);
close all ;

load ../c3po_dataStorage_particles/time_0_2.0_processor0.h5 ;
val(6)= U_Favre_0(1,1);
close all ;

load ../c3po_dataStorage_particles/time_0_2.25_processor0.h5 ;
val(7)= U_Favre_0(1,1);
close all ;

load ../c3po_dataStorage_particles/time_0_2.5_processor0.h5 ;
val(8)= U_Favre_0(1,1);
close all ;

load ../c3po_dataStorage_particles/time_0_2.75_processor0.h5 ;
val(9)= U_Favre_0(1,1);
close all ;

load ../c3po_dataStorage_particles/time_0_3.0_processor0.h5 ;
val(10)= U_Favre_0(1,1);
close all ;



for i=1:10

sol(i)= 1 - 3/2 * (x(i)*x(i)-1)/(x(i)*x(i)*x(i)-1);

i=i+1;
end

plot(x,val,'o',x,sol,'--')

set(gca,'Fontsize',14)
legend('CPPPO','Analytical solution')
xlabel('R_f/R')
ylabel('U_f/U')

print('-depsc','validation.eps')

clear;
clc;
close all;

coarse = 'data/stokes_coarse.dat';
c = dlmread(coarse);
t = c(:,1);
z = c(:,6);

medium = 'data/stokes_medium.dat';
m = dlmread(medium);
u = m(:,1);

r = linspace(0.5,10,100);

stokesSol = linspace(1,4.5,100);

for i=1:100
 stokesSol(i) = 1 - 0.5*0.5*0.5/( 4* r(i)*r(i)*r(i) ) - 3*0.5/(4*r(i));

i=i+1;
end

plot(z,t,'o',z,u,'+',r,stokesSol);

set(gca,'Fontsize',14)
legend('coarse grid','fine grid','Analytical solution')
title('Stokes Flow past a sphere');
xlabel('r');
ylabel('U (theta = p/2)');
grid on
print('-depsc','solution.eps')

