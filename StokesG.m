clear all
close all

%2D Stokes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pressure Differential (Right side is 100)
p0 = 0;

%Viscosity
mu = 2;

%grav. constant
g = 9.8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute Stokes2DG with M=50 and use it as the "analytical" solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num = 30;
%[P U V X Y] = Stokes2DG(g, num, p0, mu, 1);
[P U V X Y] = Stokes2D(num, p0, mu, 1);

%Compute the L2E with 20 different grid spacing values
testPoints = 20;

L2EP = zeros(testPoints,1);
L2EU = zeros(testPoints,1);
L2EV = zeros(testPoints,1);
d = zeros(testPoints,1);

% Compute the L2 error
%2D: || u(x,y) || = sqrt( 1/M^2 sum_{j=1}^M sum_{k=1]^M u_{j,k}^2 }
for i = 4:testPoints
    
    %[ p u v x y] = Stokes2DG(g, i, p0, mu, 0);
    [ p u v x y] = Stokes2D(i, p0, mu, 0);
    d(i) = (1-0) / (i - 1);
    
    for j=1:i
        for k = 1:i
            % interp2(X,Y,F,x,y) interpolates F(X,Y) at x,y
            
            % Don't Want:  X is numx1 vector, Y is a numx1, P is a num^2x1 vector
            % Want: X is numxnum array, Y is numxnum array, P is a numxnum array
            L2EP(i) = L2EP(i) + (p(j,k) - interp2(X,Y,P,x(j,k),y(j,k), 'cubic'))^2;
            L2EU(i) = L2EU(i) + (u(j,k) - interp2(X,Y,U,x(j,k),y(j,k), 'cubic'))^2;
            L2EV(i) = L2EV(i) + (v(j,k) - interp2(X,Y,V,x(j,k),y(j,k), 'cubic'))^2;
        end
    end
    L2EP(i) = sqrt(L2EP(i) / (i^2));
    L2EU(i) = sqrt(L2EU(i) / (i^2));
    L2EV(i) = sqrt(L2EV(i) / (i^2));
end

d = d(4:end);
L2EP = L2EP(4:end);
L2EU = L2EU(4:end);
L2EV = L2EV(4:end);

figure(4)
loglog(d,L2EP,'-',d,d.^2,'--');
xlabel('\Delta x');
ylabel('Discrete L2 Error (p)');
%title('L2 Error for P (Pressure)');

figure(5)
loglog(d,L2EU,'-',d,d.^2,'--');
xlabel('\Delta x');
ylabel('Discrete L2 Error (u)');
%title('L2 Error for U (Horizontal Velocity)');

figure(6)
loglog(d,L2EV,'-',d,d.^2,'--');
xlabel('\Delta x');
ylabel('Discrete L2 Error (v)');
%title('L2 Error for V (Vertical Velocity)');

figure(7);
quiver(X(1:4:end,1:4:end),Y(1:4:end,1:4:end),U(1:4:end,1:4:end),V(1:4:end,1:4:end));
streamline(X,Y,U,V,1.0,0.7);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Now compute the L2 error with the actual analytical solution
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% testPoints = 20;
% 
% L2EP = zeros(testPoints,1);
% L2EU = zeros(testPoints,1);
% L2EV = zeros(testPoints,1);
% d = zeros(testPoints,1);
% 
% % Compute the L2 error
% %2D: || u(x,y) || = sqrt( 1/M^2 sum_{j=1}^M sum_{k=1]^M u_{j,k}^2 }
% for i = 4:testPoints
%     
%     [ p u v x y] = Stokes2DG(g, i, p0, mu, 0);
%     d(i) = (1-0) / (i - 1);
%     
%     for j=1:i
%         for k = 1:i
%             L2EP(i) = L2EP(i) + (p(j,k) - (x(j,k) * (100 - p0) + p0))^2;
%             L2EU(i) = L2EU(i) + (u(j,k) - ((1/(2*mu)) * (100 - p0) * y(j,k) * (y(j,k) - 1)))^2;
%             %L2EV(i) = L2EV(i) + (v(j,k) - ?))^2;
%         end
%     end
%     L2EP(i) = sqrt(L2EP(i) / (i^2));
%     L2EU(i) = sqrt(L2EU(i) / (i^2));
%     %L2EV(i) = sqrt(L2EV(i) / (i^2));
% end
% 
% d = d(4:end);
% L2EP = L2EP(4:end);
% L2EU = L2EU(4:end);
% L2EV = L2EV(4:end);
% 
% figure(7)
% loglog(d,L2EP,'-',d,d.^2,'--');
% title('L2 Error for P (Pressure)');
% 
% figure(8)
% loglog(d,L2EU,'-',d,d.^2,'--');
% title('L2 Error for U (Horizontal Velocity)');
% 
% %figure(9)
% %loglog(d,L2EV,'-',d,d.^2,'--');
% %title('L2 Error for V (Vertical Velocity)');





