% Simulation of a 1-D Linear Convection model by a time march (Finite
...Difference Method)
% Numerical scheme used is a first order upwind in both space and time
%%
%Specifying Parameters
width = 1;
timerange = 1;
dt=0.008;              %Width of each time step
dx = 0.01;
nt = timerange / dt;
nx = width / dx + 1;
a=1.0;                %Velocity of wave propagation
% dx=1/(nx-1);          %Width of space step
x=-0.5:dx:0.5;             %Range of x (0,2) and specifying the grid points
u=zeros(1,nx);        %Preallocating u
u1=zeros(1,nx);        %Preallocating u1
u2=zeros(1,nx);        %Preallocating u2
un=zeros(1,nx);       %Preallocating un
sigma=abs(a)*dt/dx;   %Courant-Freidrich-Lewy number
epsilon2 = 0.2;
%%
%Initial Conditions: A square wave
u = f2(x, nx);
%%
%Evaluating velocity profile for each time step

for it=0:nt
    un=u;
    u_exact = f2(x-dt*it, nx);
    
    %% plot the wave
    h1 = plot(x,u, 'color', [216, 50, 48]/256, 'DisplayName','Numerical solution');       %plotting the velocity profile
    legend('Numerical solution')
    axis([-1 1 -2 2])
    xlabel('x')
    ylabel('u')
    title({['Third-order Runge-Kutta scheme in time, second-order central difference scheme in space']; ['with {\itc} = ',num2str(a),' and CFL (\sigma) = ',num2str(sigma)];['time(\itt) = ',num2str(dt*it)]})
    hold on
    h2 = plot(x, u_exact, 'color', [76, 198, 252]/256, 'DisplayName','Exact solution');
    legend('show', 'Location','northwest');
    drawnow;
    pause(0.1)
    delete([h1,h2])
%% Three order Runge-Kutta
    for j=2:nx-1
        u1(j) = un(j) - dt/dx*(un(j+1)-un(j-1))/2 + sign(a) * epsilon2*dt / dx * (un(j+1)-2*un(j)+un(j-1));
    end
    u1(1) = un(1) - dt/dx*(un(2)-un(end))/2 + sign(a) * epsilon2 *dt/ dx * (un(2)-2*un(1)+un(end));
    u1(end) = un(end) - dt/dx*(un(1)-un(end-1))/2 + sign(a) * epsilon2 *dt/ dx * (un(1)-2*un(end)+un(end-1));

    for j=2:nx-1
        u2(j) = 3/4*un(j)+ 1/4*u1(j) - 1/4*dt/dx*(u1(j+1)-u1(j-1))/2 + 1/4*sign(a) * epsilon2*dt / dx * (u1(j+1)-2*u1(j)+u1(j-1));
    end
    u2(1) = 3/4*un(1)+ 1/4*u1(1) - dt/4/dx*(u1(2)-u1(end))/2+ 1/4*sign(a) * epsilon2*dt / dx * (u1(2)-2*u1(1)+u1(end));
    u2(end) = 3/4*un(end)+ 1/4*u1(end) - dt/4/dx*(u1(1)-u1(end-1))/2+ 1/4*sign(a) * epsilon2*dt / dx * (u1(1)-2*u1(end)+u1(end-1));
    
    for j=2:nx-1
        u(j) = 1/3*un(j)+2/3*u2(j) - 2/3*dt/dx * (u2(j+1)-u2(j-1))/2 + 2/3*sign(a) * epsilon2*dt / dx * (u2(j+1)-2*u2(j)+u2(j-1));
    end
    u(1) = 1/3*un(1)+2/3*u2(1) - 2/3*dt/dx * (u2(2)-u2(end))/2 +2/3*sign(a) * epsilon2*dt / dx * (u2(2)-2*u2(1)+u2(end));
    u(end) = 1/3*un(end)+2/3*u2(end) - 2/3*dt/dx * (u2(1)-u2(end-1))/2 +2/3*sign(a) * epsilon2*dt / dx * (u2(1)-2*u2(end)+u2(end-1));
%     u=u1;
end
%% function to get exact solution
function u_exact = f(x, nx)
   u_exact=zeros(1,nx);       %Preallocating u_exact
   for j=1:nx
      flag = mod(x(j), 1);
      if flag >0.5
          flag = flag - 1;
      end
    if ((-0.5<=flag)&&(flag<-0.25))
        u_exact(j)=0;
    elseif ((-0.25<=flag)&& (flag<=0.25))
        u_exact(j)=1;
    elseif ((0.25<flag) && (flag<=0.5))
        u_exact(j)=0;
    end
   end
end

%% function to get exact solution
function u_exact = f2(x, nx)
   u_exact=zeros(1,nx);       %Preallocating u_exact
   for j=1:nx
      flag = mod(x(j), 1);
      if flag >0.5
          flag = flag - 1;
      end
    u_exact(j) = sin(4*pi*flag);
   end
end