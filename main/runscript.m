% Cell apostasis

D_GFP = 24;

D_n = @(D_GFP, m_GFP, m_protein) D_GFP * (m_GFP/m_protein)^(1/3);
D = 2;

% Grid
L = 10;
N = L*10;
x = linspace(1,30,N);
y = x;
[X,Y] = meshgrid(x,y);
dx = L/N;
dt = 0.001;
% time
nmax = 1000;

% Initialize the concentration matrix

% random initial conditions
c_n = rand(N, N);

% half alive initial conditions
% c_n = zeros(N,N);
% c_n(N/2+1:end,:) = ones(N/2,N);

figure;
axis square;
for n = 1:nmax
    c_n(:,:,n+1) = rk2(c_n(:,:,n),D,dx,@rhsdiff,dt);
    imagesc(c_n(:,:,n)); colorbar; colormap bone; axis square;
    drawnow;
end

function k = rhsdiff(D,c_n,dx)
N = length(c_n(:,1));
k = zeros(size(c_n));
for i = 2:N-1
    for j=2:N-1
        k(i,j) = D * (1/dx^2) * (c_n(i+1,j) + c_n(i-1,j) + ...
            c_n(i,j+1) + c_n(i,j-1) - 4*c_n(i,j));
    end
end
% zero flux BCs
k(1,:) = k(2,:);
k(end-1,:) = k(end,:);
k(:,1) = k(:,2);
k(:,end-1) = k(:,end);
end

function u = rk2(u,D,dx,dynfunc,dt)
    k1 = dynfunc(D,u,dx);
    k2 = dynfunc(D,u+k1*dt,dx);
    u = u + dt*((k1+k2)/2);
    u(1,:) = u(2,:);
u(end-1,:) = u(end,:);
u(:,1) = u(:,2);
u(:,end-1) = u(:,end);
end