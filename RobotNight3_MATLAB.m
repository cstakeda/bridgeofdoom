%clear

%These equations were derived before I understood the need to sub u with
%a*t

%Timesteps
t = linspace(0,3.2,200);

%Original function for spiral of doom 
r = [0.396E0.*cos(0.371E1+0.265E1.*t); (-0.99E0).*sin(0.14E1+t); 0*t];

%Tangent vector
That = [(-0.10494E1).*sin(0.371E1+0.265E1.*t).*(0.9801E0.*cos(0.14E1+t).^2+0.110124E1.*sin(0.371E1+0.265E1.*t).^2).^(-1/2);
    (-0.99E0).*cos(0.14E1+t).*(0.9801E0.*cos(0.14E1+t).^2+0.110124E1.*sin(0.371E1+0.265E1.*t).^2).^(-1/2);
    0*t];

%Normal vector
Nhat = [cos(0.14E1+t).*((-0.278091E1).*cos(0.14E1+t).*cos(0.371E1+0.265E1.*t)+(-0.10494E1).*sin(0.14E1+t).*sin(0.371E1+0.265E1.*t)).*(0.9801E0.*cos(0.14E1+t).^2+0.110124E1.*sin(0.371E1+0.265E1.*t).^2).^(-1/2).*(0.1E1.*cos(0.14E1+t).^2+0.11236E1.*sin(0.371E1+0.265E1.*t).^2).^(-1).*((0.1E1.*cos(0.14E1+t).^2+0.11236E1.*sin(0.371E1+0.265E1.*t).^2).^(-3).*(cos(0.14E1+t).^4.*(0.789048E1.*cos(0.371E1+0.265E1.*t).^2+(-0.113994E-15).*sin(0.14E1+t).^2)+0.886574E1.*cos(0.14E1+t).^2.*cos(0.371E1+0.265E1.*t).^2.*sin(0.371E1+0.265E1.*t).^2+0.2809E0.*sin(0.28E1+2.*t).^2.*sin(0.371E1+0.265E1.*t).^2+0.126248E1.*sin(0.14E1+t).^2.*sin(0.371E1+0.265E1.*t).^4+0.372192E0.*csc(0.14E1+t).^2.*sin(0.28E1+2.*t).^3.*sin(0.742E1+0.53E1.*t)+0.167278E1.*sin(0.28E1+2.*t).*sin(0.371E1+0.265E1.*t).^2.*sin(0.742E1+0.53E1.*t))).^(-1/2);
sin(0.371E1+0.265E1.*t).*(0.294776E1.*cos(0.14E1+t).*cos(0.371E1+0.265E1.*t)+0.111236E1.*sin(0.14E1+t).*sin(0.371E1+0.265E1.*t)).*(0.9801E0.*cos(0.14E1+t).^2+0.110124E1.*sin(0.371E1+0.265E1.*t).^2).^(-1/2).*(0.1E1.*cos(0.14E1+t).^2+0.11236E1.*sin(0.371E1+0.265E1.*t).^2).^(-1).*((0.1E1.*cos(0.14E1+t).^2+0.11236E1.*sin(0.371E1+0.265E1.*t).^2).^(-3).*(cos(0.14E1+t).^4.*(0.789048E1.*cos(0.371E1+0.265E1.*t).^2+(-0.113994E-15).*sin(0.14E1+t).^2)+0.886574E1.*cos(0.14E1+t).^2.*cos(0.371E1+0.265E1.*t).^2.*sin(0.371E1+0.265E1.*t).^2+0.2809E0.*sin(0.28E1+2.*t).^2.*sin(0.371E1+0.265E1.*t).^2+0.126248E1.*sin(0.14E1+t).^2.*sin(0.371E1+0.265E1.*t).^4+0.372192E0.*csc(0.14E1+t).^2.*sin(0.28E1+2.*t).^3.*sin(0.742E1+0.53E1.*t)+0.167278E1.*sin(0.28E1+2.*t).*sin(0.371E1+0.265E1.*t).^2.*sin(0.742E1+0.53E1.*t))).^(-1/2);
0*t];

%Linear speed
nu = (0.9801E0.*cos(0.14E1+t).^2+0.110124E1.*sin(0.371E1+0.265E1.*t).^2).^(1/2);

%Angular velocity
w = (0.1E1.*cos(0.14E1+t).^2+0.11236E1.*sin(0.371E1+0.265E1.*t).^2).^(-3).*((-0.2809E1).*cos(0.14E1+t).^5.*cos(0.371E1+0.265E1.*t)+(-0.106E1).*cos(0.14E1+t).^4.*sin(0.14E1+t).*sin(0.371E1+0.265E1.*t)+(-0.631238E1).*cos(0.14E1+t).^3.*cos(0.371E1+0.265E1.*t).*sin(0.371E1+0.265E1.*t).^2+(-0.238203E1).*cos(0.14E1+t).^2.*sin(0.14E1+t).*sin(0.371E1+0.265E1.*t).^3+(-0.35463E1).*cos(0.14E1+t).*cos(0.371E1+0.265E1.*t).*sin(0.371E1+0.265E1.*t).^4+(-0.133823E1).*sin(0.14E1+t).*sin(0.371E1+0.265E1.*t).^5);

%Left and right velocities
VL = (0.9801E0.*cos(0.14E1+t).^2+0.110124E1.*sin(0.371E1+0.265E1.*t).^2).^(1/2)+(0.1E1.*cos(0.14E1+t).^2+0.11236E1.*sin(0.371E1+0.265E1.*t).^2).^(-3).*(0.356743E0.*cos(0.14E1+t).^5.*cos(0.371E1+0.265E1.*t)+0.13462E0.*cos(0.14E1+t).^4.*sin(0.14E1+t).*sin(0.371E1+0.265E1.*t)+0.801673E0.*cos(0.14E1+t).^3.*cos(0.371E1+0.265E1.*t).*sin(0.371E1+0.265E1.*t).^2+0.302518E0.*cos(0.14E1+t).^2.*sin(0.14E1+t).*sin(0.371E1+0.265E1.*t).^3+0.45038E0.*cos(0.14E1+t).*cos(0.371E1+0.265E1.*t).*sin(0.371E1+0.265E1.*t).^4+0.169955E0.*sin(0.14E1+t).*sin(0.371E1+0.265E1.*t).^5);
VR = (0.9801E0.*cos(0.14E1+t).^2+0.110124E1.*sin(0.371E1+0.265E1.*t).^2).^(1/2)+(0.1E1.*cos(0.14E1+t).^2+0.11236E1.*sin(0.371E1+0.265E1.*t).^2).^(-3).*((-0.356743E0).*cos(0.14E1+t).^5.*cos(0.371E1+0.265E1.*t)+(-0.13462E0).*cos(0.14E1+t).^4.*sin(0.14E1+t).*sin(0.371E1+0.265E1.*t)+(-0.801673E0).*cos(0.14E1+t).^3.*cos(0.371E1+0.265E1.*t).*sin(0.371E1+0.265E1.*t).^2+(-0.302518E0).*cos(0.14E1+t).^2.*sin(0.14E1+t).*sin(0.371E1+0.265E1.*t).^3+(-0.45038E0).*cos(0.14E1+t).*cos(0.371E1+0.265E1.*t).*sin(0.371E1+0.265E1.*t).^4+(-0.169955E0).*sin(0.14E1+t).*sin(0.371E1+0.265E1.*t).^5);

%% Deliverable 1 animated

figure(1)
for n = 1:length(t)% loop through each of the points
    plot3(r(1,:),r(2,:),r(3,:)), axis ([-1.4 1.4 -1.4 1.4 -0.1 0.1]),hold on % plot the entire curve
    quiver3(r(1,n),r(2,n),r(3,n),That(1,n),That(2,n),That(3,n),'r') % plot the unit tangent
    quiver3(r(1,n),r(2,n),r(3,n),Nhat(1,n),Nhat(2,n),Nhat(3,n),'b') % plot the unit normal
   hold off
   view(2)
   drawnow
end
%% Deliverable 1 at time step 50,100,150

figure(2)

 plot3(r(1,:),r(2,:),r(3,:),'black'), axis ([-1.4 1.4 -1.4 1.4 -0.1 0.1]),hold on % plot the entire curve
 quiver3(r(1,25),r(2,100),r(3,100),That(1,100),That(2,100),That(3,100),'r') % plot the unit tangent
 quiver3(r(1,25),r(2,100),r(3,100),Nhat(1,100),Nhat(2,100),Nhat(3,100),'b') % plot the unit normal
 
 quiver3(r(1,50),r(2,50),r(3,50),That(1,50),That(2,50),That(3,50),'r') % plot the unit tangent
 quiver3(r(1,50),r(2,50),r(3,50),Nhat(1,50),Nhat(2,50),Nhat(3,50),'b') % plot the unit normal
 
 quiver3(r(1,150),r(2,150),r(3,150),That(1,150),That(2,150),That(3,150),'r') % plot the unit tangent
 quiver3(r(1,150),r(2,150),r(3,150),Nhat(1,150),Nhat(2,150),Nhat(3,150),'b') % plot the unit normal
 
 view(2)
 title('Unit Tangent and Unit Normal Vectors at Varying Points')
 xlabel('X (m)')
 ylabel('Y (m)')


%% Deliverable 2

figure(3)
hold on
plot(t,VL,'b')
plot(t,VR,'r')
legend('VL','VR')
title('Left and Right Velocities as a Function of Time')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')

%% Optimizing alpha for a desired maxiumum wheel velocity

VL = (a.^2.*(0.9801E0.*cos(0.14E1+t.*a).^2+0.110124E1.*sin(0.371E1+0.265E1.*t.*a).^2)).^(1/2)+0.356743E0.*a.*(0.1E1.*cos(0.14E1+t.*a).^2+0.11236E1.*sin(0.371E1+0.265E1.*t.*a).^2).^(-3).*(0.1E1.*cos(0.14E1+t.*a).^5.*cos(0.371E1+0.265E1.*t.*a)+0.377358E0.*cos(0.14E1+t.*a).^4.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a)+0.22472E1.*cos(0.14E1+t.*a).^3.*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^2+0.848E0.*cos(0.14E1+t.*a).^2.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^3+0.126248E1.*cos(0.14E1+t.*a).*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^4+0.476406E0.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^5);
VR = (a.^2.*(0.9801E0.*cos(0.14E1+t.*a).^2+0.110124E1.*sin(0.371E1+0.265E1.*t.*a).^2)).^(1/2)+(-0.356743E0).*a.*(0.1E1.*cos(0.14E1+t.*a).^2+0.11236E1.*sin(0.371E1+0.265E1.*t.*a).^2).^(-3).*(0.1E1.*cos(0.14E1+t.*a).^5.*cos(0.371E1+0.265E1.*t.*a)+0.377358E0.*cos(0.14E1+t.*a).^4.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a)+0.22472E1.*cos(0.14E1+t.*a).^3.*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^2+0.848E0.*cos(0.14E1+t.*a).^2.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^3+0.126248E1.*cos(0.14E1+t.*a).*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^4+0.476406E0.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^5);

a_test = linspace(0.001,.21,500);
curr = 1;
best_a = zeros(1,length(a_test));

for j = 0.001:length(a_test)
    
    a = a_test(1,curr);
    newt = 3.2/a;
    t = linspace(0,newt,200);
    VL = (a.^2.*(0.9801E0.*cos(0.14E1+t.*a).^2+0.110124E1.*sin(0.371E1+0.265E1.*t.*a).^2)).^(1/2)+0.356743E0.*a.*(0.1E1.*cos(0.14E1+t.*a).^2+0.11236E1.*sin(0.371E1+0.265E1.*t.*a).^2).^(-3).*(0.1E1.*cos(0.14E1+t.*a).^5.*cos(0.371E1+0.265E1.*t.*a)+0.377358E0.*cos(0.14E1+t.*a).^4.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a)+0.22472E1.*cos(0.14E1+t.*a).^3.*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^2+0.848E0.*cos(0.14E1+t.*a).^2.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^3+0.126248E1.*cos(0.14E1+t.*a).*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^4+0.476406E0.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^5);
    VR = (a.^2.*(0.9801E0.*cos(0.14E1+t.*a).^2+0.110124E1.*sin(0.371E1+0.265E1.*t.*a).^2)).^(1/2)+(-0.356743E0).*a.*(0.1E1.*cos(0.14E1+t.*a).^2+0.11236E1.*sin(0.371E1+0.265E1.*t.*a).^2).^(-3).*(0.1E1.*cos(0.14E1+t.*a).^5.*cos(0.371E1+0.265E1.*t.*a)+0.377358E0.*cos(0.14E1+t.*a).^4.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a)+0.22472E1.*cos(0.14E1+t.*a).^3.*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^2+0.848E0.*cos(0.14E1+t.*a).^2.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^3+0.126248E1.*cos(0.14E1+t.*a).*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^4+0.476406E0.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^5);
    maxVL = max(VL);
    maxVR = max(VR);
    
    if maxVL < 3
        if maxVR < 3
            best_a(1,curr) = a;
        end
    end
               
    curr = curr+1;
end
    
%% Recalculating values based upon optimized value of alpha

a = 0.204555110220441;
newt = 3.2/a;
t = linspace(0,newt,100);
VL = (a.^2.*(0.9801E0.*cos(0.14E1+t.*a).^2+0.110124E1.*sin(0.371E1+0.265E1.*t.*a).^2)).^(1/2)+0.356743E0.*a.*(0.1E1.*cos(0.14E1+t.*a).^2+0.11236E1.*sin(0.371E1+0.265E1.*t.*a).^2).^(-3).*(0.1E1.*cos(0.14E1+t.*a).^5.*cos(0.371E1+0.265E1.*t.*a)+0.377358E0.*cos(0.14E1+t.*a).^4.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a)+0.22472E1.*cos(0.14E1+t.*a).^3.*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^2+0.848E0.*cos(0.14E1+t.*a).^2.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^3+0.126248E1.*cos(0.14E1+t.*a).*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^4+0.476406E0.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^5);
VR = (a.^2.*(0.9801E0.*cos(0.14E1+t.*a).^2+0.110124E1.*sin(0.371E1+0.265E1.*t.*a).^2)).^(1/2)+(-0.356743E0).*a.*(0.1E1.*cos(0.14E1+t.*a).^2+0.11236E1.*sin(0.371E1+0.265E1.*t.*a).^2).^(-3).*(0.1E1.*cos(0.14E1+t.*a).^5.*cos(0.371E1+0.265E1.*t.*a)+0.377358E0.*cos(0.14E1+t.*a).^4.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a)+0.22472E1.*cos(0.14E1+t.*a).^3.*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^2+0.848E0.*cos(0.14E1+t.*a).^2.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^3+0.126248E1.*cos(0.14E1+t.*a).*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^4+0.476406E0.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^5);
w = a.*(0.1E1.*cos(0.14E1+t.*a).^2+0.11236E1.*sin(0.371E1+0.265E1.*t.*a).^2).^(-3).*((-0.2809E1).*cos(0.14E1+t.*a).^5.*cos(0.371E1+0.265E1.*t.*a)+(-0.106E1).*cos(0.14E1+t.*a).^4.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a)+(-0.631238E1).*cos(0.14E1+t.*a).^3.*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^2+(-0.238203E1).*cos(0.14E1+t.*a).^2.*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^3+(-0.35463E1).*cos(0.14E1+t.*a).*cos(0.371E1+0.265E1.*t.*a).*sin(0.371E1+0.265E1.*t.*a).^4+(-0.133823E1).*sin(0.14E1+t.*a).*sin(0.371E1+0.265E1.*t.*a).^5);
nu = (a.^2.*(0.9801E0.*cos(0.14E1+t.*a).^2+0.110124E1.*sin(0.371E1+0.265E1.*t.*a).^2)).^(1/2);
That = [(-0.10494E1).*a.*sin(0.371E1+0.265E1.*t.*a).*(a.^2.*(0.9801E0.*cos(0.14E1+t.*a).^2+0.110124E1.*sin(0.371E1+0.265E1.*t.*a).^2)).^(-1/2);
    (-0.99E0).*a.*cos(0.14E1+t.*a).*(a.^2.*(0.9801E0.*cos(0.14E1+t.*a).^2+0.110124E1.*sin(0.371E1+0.265E1.*t.*a).^2)).^(-1/2);
    0*t];

%% Deliverable 3

load('wheel_encoder.mat') %Loading data from encoder data

%Finding the change in time and change in position to find the velocity
time = dataset(1:201,1);
timestep = diff(dataset(:,1));
VLexp = diff(dataset(:,2));
VRexp = diff(dataset(:,3));
VLacc = (VLexp(1:201,1))./timestep(1:201,1);
VRacc = (VRexp(1:201,1))./timestep(1:201,1);

Vacc = (VLacc+VRacc)/2;
wacc = (VRacc-VLacc)/0.254;

%Scale factor to properly have timesteps of theoretical and actual line up
scal1 = max(time)/max(t);

theta = zeros(length(time),1);
x1 = zeros(length(time),1);
y1 = zeros(length(time),1);

for i=2:length(time)
   theta(i) = theta(i-1) + wacc(i)*timestep(i);
   x1(i) = x1(i-1) + (Vacc(i)*timestep(i))*cos(theta(i));
   y1(i) = y1(i-1) + (Vacc(i)*timestep(i))*sin(theta(i));

end

racc = [x1, y1];
tanacc = normr(diff(r));

a = 2;
newt = 3.2/a;
t = linspace(0,newt,50);
r = [0.396E0.*cos(0.371E1+0.265E1.*t.*a);(-0.99E0).*sin(0.14E1+t.*a); 0*t];

figure(5)
hold on
plot(r(1,:),r(2,:),'r')
rotate(plot((x1/1.14)+.2,y1-1,'b--'),[0 0 1],-28)
plot(tanacc(:,1),tanacc(:,2))
legend('Theoretical','Actual')
title('Theoretical and Actual Paths Followed by Neato')


%% Plotting various actual and theoretical values for Neato

figure(6)%Linear speed
hold on
plot(t,nu,'black')
plot(time/scal1,Vacc,'r')
title('Linear Speed as a Function of Time')
xlabel('Time (s)')
ylabel('Linear Speed (m/s)')
legend('Theoretical','Actual')

figure(7)%Angular velocity
hold on
plot(t,w,'m')
plot(time/scal1,wacc,'b')
title('Angular Velocity as a Function of Time')
xlabel('Time (s)')
ylabel('Angular Velocity (radians/s)')
legend('Theoretical','Actual')

figure(8)
hold on
plot(t,VL,'r')
plot(time/scal1,VLacc,'m')
hold off

figure(9)
hold on
plot(t,VR,'b')
plot(time/scal1,VRacc,'c')
legend('Left Velocity','Right Velocity')
title('Left and Right Velocities as a Function of Time')
xlabel('Time (sec)')
ylabel('Velocity (m/s)')
hold off

figure(10)
hold on
plot(time/scal1,VLacc,'b')
plot(time/scal1,VRacc,'c')


%% Programming the neato to follow my desired path

function doom(t,VL,VR)
%This function is an absolute piece of dogshit

pubvel = rospublisher('/raw_vel');
message = rosmessage(pubvel); 
time=1;

tic
while time < length(t)
    pause(0.01)
    if toc >= t(1,time+1)
        time = time+1;
        message.Data = [VL(1,time),VR(1,time)];
        send(pubvel,message);
    end
end  

end