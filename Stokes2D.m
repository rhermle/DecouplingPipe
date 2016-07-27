function [ P U V X Y] = Stokes2D( M, p0, mu, toGraph )

%2D Laplace
% p(0,y) = p0, p(1,y) = 100, p_y(x,0) = 0, p_y(x,1) = 0

% Ap=b

%Vectors p and A are indexed as:
% 1 
% 2
% 3
% 4 
% 5
% 6
% 7
% 8
% 9
% 10
% 11
% 12
% 13
% 14
% 15
% 16

%p:
% 1,1 
% 2,1
% 3,1
% 4,1
% 1,2 
% 2,2
% 3,2
% 4,2
% 1,3 
% 2,3
% 3,3
% 4,3
% 1,4  
% 2,4
% 3,4
% 4,4

%The following description is for M=4:
%A and p are indexed from 1 to 16 based on the subscripts for p above, in
%order to create 16 different equations.
%The second index for A is another 1 to 16 based on the same subscripts so
%each equation (row) for A matches up with the elements of p.

%p(i-1,j) - 4p_(i,j) + p(i+1,j) + p(i,j+1) + p(i,j-1) = 0


%M=30;

A = zeros(M*M,M*M);
b = zeros(M*M,1);

%interior points
for i=2:M-1
    for j=2:M-1
        
        index = sub2ind([M,M],i,j);
        
        A(index, sub2ind([M,M], i-1,j)) = 1;
        A(index, index) = -4;
        A(index, sub2ind([M,M], i+1,j)) = 1;
        A(index, sub2ind([M,M], i,j+1)) = 1;
        A(index, sub2ind([M,M], i,j-1)) = 1;

        b(index) = 0;        
    end
end

%border conditions

%p(1,j) = p0
%p(M, j) = 100
for j = 1:M
    index1 = sub2ind([M,M],1,j);
    indexM = sub2ind([M,M],M,j);
    
    A(index1, index1) = 1;
    A(indexM, indexM) = 1;
    b(index1) = p0;
    b(indexM) = 100;
end

%O(2) forward difference formula for bottom endpoint
%py(i,1) = 0 = -3p(i,1)+4p(i,2)-1p(i,3)
%O(2) backward difference formula for top endpoint
%py(i,M) = 0 = 3p(i,M)-4p(i,M-1)+1p(i,M-2)
for i = 2:M-1
    index1 = sub2ind([M,M], i,1);
    indexM = sub2ind([M,M], i,M);
    
    A(index1, index1) = -3;
    A(index1, sub2ind([M,M], i,2)) = 4;
    A(index1, sub2ind([M,M], i,3)) = -1;
    b(index1) = 0;
    
    A(indexM, indexM) = 3;
    A(indexM, sub2ind([M,M], i,M-1)) = -4;
    A(indexM, sub2ind([M,M], i,M-2)) = 1;
    b(indexM) = 0;
end


p = A\b;
%rcond(A) OR use cond

P = zeros(M,M);

for i = 1:M
   for j = 1:M
      k = M * (i-1);
      P(i,j) = p(j + k); 
   end
end

if(toGraph)
    %transform p into a 2x2 matrix for graphing

    X = linspace(0,1,M);
    Y = linspace(0,1,M);

    figure(1);
    surf(X,Y,P);
    %title('Pressure');
    zlabel('p');
    xlabel('x');
    ylabel('y');
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Poisson Equation:
% u_xx + u_yy = p_x / mu
% u_x(0,y) = 0
% u_x(1,y) = 0
% u(x,0) = 0
% u(x,1) = 0

%Vectors u and A are indexed as:

% 1 
% 2
% 3
% 4 
% 5
% 6
% 7
% 8
% 9
% 10
% 11
% 12
% 13
% 14
% 15
% 16

%p:
% 1,1 
% 2,1
% 3,1
% 4,1
% 1,2 
% 2,2
% 3,2
% 4,2
% 1,3 
% 2,3
% 3,3
% 4,3
% 1,4  
% 2,4
% 3,4
% 4,4


%The following description is for M=4:
%A is first indexed from 1 to 16 based on the subscripts for u above, in
%order to create 16 different equations.
%The second index for A is another 1 to 16 based on the same subscripts so
%each equation (row) for A matches up with the elements of u.

%u(i-1,j) - 4u_(i,j) + u(i+1,j) + u(i,j+1) + u(i,j-1) = d*(p(i+1,j)-p(i-1,j))/2

%M=30;

d = (1-0)/(M-1);
A = zeros(M*M,M*M);
b = zeros(M*M,1);

%interior points
for i=2:M-1
    for j=2:M-1
        
        index = sub2ind([M,M],i,j);
        
        A(index, sub2ind([M,M], i-1,j)) = 1;
        A(index, sub2ind([M,M], i,j)) = -4;
        A(index, sub2ind([M,M], i+1,j)) = 1;
        A(index, sub2ind([M,M], i,j+1)) = 1;
        A(index, sub2ind([M,M], i,j-1)) = 1;

        b(index) = d * (p(sub2ind([M,M], i+1,j)) - p(sub2ind([M,M], i-1,j)))/(2 * mu);      
    end
end

%border conditions

% u(x,0) = 0
% u(x,1) = 0
for i = 1:M
    index1 = sub2ind([M,M],i,1);
    indexM = sub2ind([M,M],i,M);
    
    A(index1, index1) = 1;
    A(indexM, indexM) = 1;
    b(index1) = 0;
    b(indexM) = 0;
end

% u_x(0,y) = 0
% u_x(1,y) = 0
%O(2) forward difference formula for left endpoint
%ux(1,j) = 0 = -3u(1,j)+4u(2,j)-1u(3,j)
%O(2) backward difference formula for right endpoint
%ux(M,j) = 0 = 3u(M,j)-4u(M-1,j)+1u(M-2,j)
for j = 2:M-1
    index1 = sub2ind([M,M], 1,j);
    indexM = sub2ind([M,M], M,j);
    
    A(index1, index1) = -3;
    A(index1, sub2ind([M,M], 2,j)) = 4;
    A(index1, sub2ind([M,M], 3,j)) = -1;
    b(index1) = 0;
    
    A(indexM, indexM) = 3;
    A(indexM, sub2ind([M,M], M-1,j)) = -4;
    A(indexM, sub2ind([M,M], M-2,j)) = 1;
    b(indexM) = 0;
end

u = A\b;

U = zeros(M,M);

for i = 1:M
   for j = 1:M
      k = M * (i-1);
      U(i,j) = u(j + k); 
   end
end

    X = linspace(0,1,M);
    Y = linspace(0,1,M);

if(toGraph)
    %transform u into a 2x2 matrix for graphing

    figure(2);
    surf(X,Y,U);
    %title('U (Horizontal Velocity)');
    zlabel('u');
    xlabel('x');
    ylabel('y');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D Poisson Equation:
% v_xx + v_yy = p_y / mu
% v(0,y) = 0
% v(1,y) = 0
% u(x,0) = 0
% u(x,1) = 0

%M=30;

d = (1-0)/(M-1);
A = zeros(M*M,M*M);
b = zeros(M*M,1);

%interior points
for i=2:M-1
    for j=2:M-1
        
        index = sub2ind([M,M],i,j);
        
        A(index, sub2ind([M,M], i-1,j)) = 1;
        A(index, sub2ind([M,M], i,j)) = -4;
        A(index, sub2ind([M,M], i+1,j)) = 1;
        A(index, sub2ind([M,M], i,j+1)) = 1;
        A(index, sub2ind([M,M], i,j-1)) = 1;

        b(index) = d * (p(sub2ind([M,M], i,j+1)) - p(sub2ind([M,M], i,j-1)))/(2 * mu);      
    end
end

%border conditions

% v(x,0) = 0
% v(x,1) = 0
for i = 1:M
    index1 = sub2ind([M,M],i,1);
    indexM = sub2ind([M,M],i,M);
    
    A(index1, index1) = 1;
    A(indexM, indexM) = 1;
    b(index1) = 0;
    b(indexM) = 0;
end

% v(0,y) = 0
% v(1,y) = 0
for j = 1:M
    index1 = sub2ind([M,M], 1,j);
    indexM = sub2ind([M,M], M,j);
    
    A(index1, index1) = 1;
    A(indexM, indexM) = 1;
    b(index1) = 0;
    b(indexM) = 0;
end

v = A\b;

X = linspace(0,1,M);
Y = linspace(0,1,M);

[X Y] = meshgrid(X,Y);

V = zeros(M,M);

for i = 1:M
   for j = 1:M
      k = M * (i-1);
      V(i,j) = v(j + k); 
   end
end

if(toGraph)
    %transform u into a 2x2 matrix for graphing

    figure(3);
    surf(X,Y,V);
    %title('V (Vertical Velocity)');
    zlabel('v');
    xlabel('x');
    ylabel('y');
end

if(toGraph)
    figure(7)
    quiver(X,Y,U,V)
end

end

