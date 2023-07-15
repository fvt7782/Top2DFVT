% TOPOLOGY OPTIMIZATION OF 2D STRUCTURES BY ZEROTH ORDER FINITE-VOLUME THEORY
function Top2DFVT(L,H,nx,ny,volfrac,penal,frad,ft,varargin)
%______________________________________________________________USER-DEFINED
P     = -1;                                                                % applied concentrated load
E0    =  1;                                                                % Young's modulus of solid material
Emin  = E0*(1e-9);                                                         % soft (void) material stiffness to avoid singularity
nu    =  0.3;                                                              % Poisson ratio
model = 'SIMP';                                                            % penalization method: "SIMP" or "RAMP"
eta   = 1/2.6;                                                             % damping factor
move  = 0.2;                                                               % move-limit
maxit = 100;                                                               % maximum number of iterations
%__________________________________INITIALIZATION OF DESIGN VARIABLE VECTOR
l = L/nx; h = H/ny;                                                        % subvolume dimensions
x = repmat(volfrac,nx,ny);                                                 % regular distribution of density
dvdx  = repmat(l*h,nx,ny);                                                 % volume constraint gradient
fig   = PlotTopology(l,h,x);                                               % optimized topology plot initialization
axis equal off tight; colormap(gray); caxis([0 1]);                        % plot configuration
%_____________________________________PREPARE FINITE-VOLUME THEORY ANALYSIS
[dof,ndof,iK,jK] = DOFassembly(nx,ny);                                     % degrees of freedom
K0 = LocalStiffMatrix(nu,l,h);                                             % local stiffness matrix
%______________________________________________________PENALIZATION METHODS
if (strcmp(model,'SIMP'))
    MatInt = @(p,x) {Emin+(E0-Emin)*x.^p,p*(E0-Emin)*x.^(p-1)};            % SIMP - Solid Isotropic Material with Penalization
elseif (strcmp(model,'RAMP'))
    MatInt = @(p,x) {Emin+(E0-Emin)*x./(1+p*(1-x)),...                     % RAMP - Rational Approximation of Material Properties
        (1+p)*(E0-Emin)./(1+p*(1-x)).^2};
end
%__________________________ASSEMBLY LOAD VECTOR AND GLOBAL STIFFNESS MATRIX
if (isempty(varargin))
    F = sparse(dof(nx*(ny+1)/2,4)',1,P,ndof,1);                            % global force vector (resultant forces on the subvolumes' faces)
    StiffnessAssemblage = @(sK) sparse(iK(:),jK(:),sK(:),ndof,ndof);       % stiffness matrix assemblage
elseif (strcmp(varargin,'fast'))
    if (~exist("fsparse\"))                                                % logical test to check if the fsparse folder is installed
        warning('check if the fsparse is correctly set up');               % display message
        return
    end
    addpath(genpath('./fsparse'));                                         % set fsparse folder
    F = fsparse(dof(nx*(ny+1)/2,4)',1,P,[ndof,1]);                         % global force vector (resultant forces on the subvolumes' faces)
    StiffnessAssemblage = @(sK) fsparse(iK(:),jK(:),sK(:),[ndof,ndof]);    % stiffness matrix assemblage
end
supp = unique(dof(1:nx:end-nx+1,7:8));                                     % fixed degrees of freedom - supporting conditions
fdof = setdiff(dof(:),supp(:));                                            % free degrees of freedom
%_____IMPLICIT FUNCTIONS TO MODIFY THE SUBVOLUME SENSITIVITIES OR DENSITIES
if (ft == 0) %
    sensitivity = @(x,dfdx,dvdx){dfdx,dvdx};
    density = @(x) x;
elseif (ft == 1) % Sensitivity filter
    [H,Hs] = Filtering(l,h,nx,ny,frad);
    sensitivity = @(x,dfdx,dvdx){H*(x(:).*dfdx(:))./Hs./max(1e-3,x(:)),dvdx};
    density = @(x) x;
elseif (ft == 2) % Density filter
    [H,Hs] = Filtering(l,h,nx,ny,frad);
    sensitivity = @(x,dfdx,dvdx){(H*dfdx(:)./Hs),H*(dvdx(:)./Hs)};
    density = @(x) (H*x(:))./Hs;
end
%______________________________________________________OPTIMIZATION PROCESS
U = zeros(ndof,1);
tic
for i = 1:length(penal(:))
    [change,loop,xPhys] = deal(1, 0, x);
    fprintf('\nPenalty factor: %1.2f\n',penal(i));                         % print current penalty factor
    while (change > 0.01 && loop < maxit), loop = loop+1;                  % start optmization process
        %_____________________________________FINITE-VOLUME THEORY ANALYSIS
        Mat = MatInt(penal(i),xPhys);                                      % material interpolation
        E = Mat{1}; dEdx = Mat{2};
        sK = K0(:)*E(:)';                                                  % stiffness interpolation
        K = StiffnessAssemblage(sK); K = (K+K')/2;                         % assemblage of the stiffness matrix
        U(fdof) = K(fdof,fdof)\F(fdof);                                    % compute global displacements
        %_______________________________________COMPLIANCE AND SENSITIVITY
        fe = reshape(sum((U(dof)*K0).*U(dof),2),nx,ny);
        f  = sum(sum(E.*fe));                                              % objective function: compliance
        dfdx = -dEdx.*fe;                                                  % sensitivity
        %_________________________________MODIFICATION OF THE SENSITIVITIES
        sens = sensitivity(x,dfdx,dvdx);
        dfdx(:) = sens{1};
        dvdx(:) = sens{2};
        %_________________UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
        [xUpp, xLow] = deal (x + move, x - move);                          % Upp. & low. limits
        OcC = x.*((-dfdx./dvdx).^eta);                                     % Opt. parameter
        inL = [0, mean(OcC(:))/volfrac];                                      % Lag. Mul. range
        while (inL(2)-inL(1))/(inL(2)+inL(1))> 1e-3
            lmid = 0.5*(inL(2)+ inL(1));
            xnew = max(0,max(xLow,min(1,min(xUpp,OcC/lmid))));
            if mean(xnew(:))>volfrac, inL(1) = lmid; else, inL(2) = lmid; end
        end
        xPhys(:) = density(xnew);
        change = max(abs(xnew(:)-x(:)));
        x = xnew;
        %_____________________________________________________PRINT RESULTS
        fprintf('It: %i\tObjec.: %1.4f\tVol.: %1.3f\tChange: %1.3f\n',...
            loop,f,mean(xPhys(:)),change);
    end
end
%___________________________________________________________PLOT DESIGN
set(fig,'FaceColor','flat','CData',1-xPhys(:)); drawnow
clearvars;
PrintTime(toc);
%_______________________________________________ORDENING DEGREES OF FREEDOM
function [dof,ndof,iK,jK] = DOFassembly(nx,ny)
[i,j] = ndgrid(1:nx,1:ny);
q = (i+(j-1)*nx);
faces = [q(:),q(:)+j(:)+nx*(ny+1),q(:)+nx,q(:)+j(:)-1+nx*(ny+1)];
dof = zeros(nx*ny,8);
dof(:,2:2:8) = 2*faces; dof(:,1:2:7) = 2*faces-1;                          % create matrix of degrees of freedom
ndof = max(dof(:));                                                        % total number of degrees of freedom
iK = reshape(kron(dof,ones(8,1))',64*nx*ny,1);                             % iK - line mapping indice to the global stiffness location
jK = reshape(kron(dof,ones(1,8))',64*nx*ny,1);                             % jK - column mapping indice to the global stiffness location
%______________LOCAL STIFFNESS MATRIX FOR ZEROTH ORDER FINITE VOLUME THEORY
function K0 = LocalStiffMatrix(nu,l,h)
C0  = 1/(1-nu^2)*[1,nu,0;nu,1,0;0,0,(1-nu)/2];
a = [1,0;0,1;1,0;0,1;1,0;0,1;1,0;0,1];
N1 = [0,0,-1;0,-1 0]; N2 = [1,0,0;0,0,1];
N3 = [0,0,1;0,1,0]; N4 = [-1,0,0;0,0,-1];
N = [N1,zeros(2,9);zeros(2,3),N2,zeros(2,6);
    zeros(2,6),N3,zeros(2,3);zeros(2,9),N4];
A1 = [0,0,1/l,0,0,0,-1/l,0;-1/h,0,0,0,1/h,0,0,0;
    0,0,2/l^2,0,0,0,2/l^2,0;2/h^2,0,0,0,2/h^2,0,0,0;
    0,0,0,1/l,0,0,0,-1/l;0,-1/h,0,0,0,1/h,0,0;
    0,0,0,2/l^2,0,0,0,2/l^2;0,2/h^2,0,0,0,2/h^2,0,0];
E = [1,0,0,0,0,0,0,0;0,0,0,0,0,1,0,-3*h/2;0,1,0,-3*h/2,1,0,0,0;
    1,0,3*l/2,0,0,0,0,0;0,0,0,0,0,1,0,0;0,1,0,0,1,0,3*l/2,0;
    1,0,0,0,0,0,0,0;0,0,0,0,0,1,0,3*h/2;0,1,0,3*h/2,1,0,0,0;
    1,0,-3*l/2,0,0,0,0,0;0,0,0,0,0,1,0,0;0,1,0,0,1,0,-3*l/2,0];
B = N*[C0,zeros(3,9);zeros(3,3),C0,zeros(3,6);
    zeros(3,6),C0,zeros(3,3);zeros(3,9),C0]*E;
sumB = (B(1:2,:)+B(5:6,:))*l+(B(3:4,:)+B(7:8,:))*h;
ab = (sumB*A1*a)\(sumB*A1);
Ab = A1*(eye(8)-a*ab);
K0 = B*Ab;
K0 = [K0(1:2,:)*l;K0(3:4,:)*h;K0(5:6,:)*l;K0(7:8,:)*h];
%_____________________________________________________PREPARING THE FILTER
function [H,Hs] = Filtering(l,h,nx,ny,frad)
neig = ceil(frad);
numSubvolumes = nx * ny * (2 * (neig - 1) + 1)^2;
iH = ones(numSubvolumes, 1);
jH = ones(numSubvolumes, 1);
sH = zeros(numSubvolumes, 1);
k = 0;
% Compute filter values
for j1 = 1:ny
    for i1 = 1:nx
        e1 = (j1 - 1) * nx + i1;
        for j2 = max(j1 - neig, 1):min(j1 + neig, ny)
            for i2 = max(i1 - neig, 1):min(i1 + neig, nx)
                e2 = (j2 - 1) * nx + i2;
                k = k + 1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,frad-sqrt(((i1-i2)*l)^2+((j1-j2)*h)^2));
            end
        end
    end
end
% Create sparse matrix and compute row sums
H = sparse(iH, jH, sH);
Hs = sum(sparse(iH, jH, sH), 2);
%_________________________________________PLOTTING THE OPTIMIZED TOPOLOGY
function fig = PlotTopology(l,h,x)
[nx,ny] = size(x);
[k,face,vert] = deal(1,zeros(nx*ny,4),zeros(4*nx*ny,2));
for j = 1:ny
    for i = 1:nx, q = i+(j-1)*nx;
        face(q,:) = k:k+3;
        vert(k:k+3,:) = [(i-1)*l,(j-1)*h;i*l,(j-1)*h;i*l,j*h;...
            (i-1)*l,j*h];
        k = k+4;
    end
end
fig = patch('Faces',face,'Vertices',vert,'FaceVertexCData',1-x(:),...
    'FaceColor','flat','EdgeColor','none');
%_____________________________________________PRINTING THE PROCESSING TIME
function PrintTime(t)
H = fix(t/3600); M = fix((t-3600*H)/60); S = fix((t-3600*H-60*M));
disp(['Timing of iterative process: ' num2str(H) 'h ' num2str(M) 'min '...
    num2str(S) 's']);
