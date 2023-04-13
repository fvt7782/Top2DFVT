% TOPOLOGY OPTIMIZATION OF STRUCTURES BY THE STANDARD FINITE-VOLUME THEORY
function Top2DFVT(L,H,nx,ny,volfrac,penal,neig,ft,varargin)
%______________________________________________________________USER-DEFINED
P     = -1;                                                                % applied concentrated load
E0    =  1;                                                                % Young's modulus of solid material
Emin  = E0*(1e-9);                                                         % soft (void) material stiffness to avoid singularity
nu    =  0.3;                                                              % Poisson ratio
model = 'RAMP';                                                            % penalization method: "SIMP" or "RAMP"
eta   = 1/2;                                                               % damping factor
move  = 0.2;                                                               % move-limit
%__________________________________INITIALIZATION OF DESIGN VARIABLE VECTOR
l = L/nx; h = H/ny;                                                        % subvolume dimensions
x = repmat(volfrac,nx,ny);                                                 % regular distribution of density
dvdx  = l*h*ones(nx,ny);                                                   % volume constraint gradient
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
    F = sparse(dof(nx*ny-nx+1,6)',1,P,ndof,1);                             % global surface-averaged force vector
    StiffnessAssemblage = @(sK) sparse(iK(:),jK(:),sK(:),ndof,ndof);
elseif (strcmp(varargin,'fast'))                                           
    if ~exist("fsparse\")                                                  % logical testo to check if the fsparce folder is installed
        warning('check if the fsparse is correctly set up');               % displyed message
        return
    end
    addpath(genpath('./fsparse'));                                         % set fsparse folder
    F = fsparse(dof(nx*ny-nx+1,6)',1,P,[ndof,1]);                          % global surface-averaged force vector
    StiffnessAssemblage = @(sK) fsparse(iK(:),jK(:),sK(:),[ndof,ndof]);
end
supp = unique([dof(1:nx:end-nx+1,7);dof(nx,2)]);                           % fixed degrees of freedom - supporting conditions
fdof = setdiff(dof(:),supp(:));                                            % free degrees of freedom
%_____IMPLICIT FUNCTIONS TO MODIFY THE SUBVOLUME SENSITIVITIES OR DENSITIES
[H,Hs] = Filtering(l,h,nx,ny,neig);
if (ft == 1)                                                               % sensitivity filter
    sensitivity = @(x,dfdx,dvdx){H*(x(:).*dfdx(:))./Hs./max(1e-3,x(:)),dvdx};
    density = @(x) x;
elseif (ft == 2)                                                           % density filter
    sensitivity = @(x,dfdx,dvdx){(H*dfdx(:))./Hs,(H*dvdx(:))./Hs};
    density = @(x) (H*x(:))./Hs;
end
%______________________________________________________OPTIMIZATION PROCESS
U = zeros(ndof,1);
tic
for i = 1:length(penal(:))
    change = 1; loop = 0;
    xPhys = x;
    fprintf('\nPenalty factor: %1.2f\n',penal(i));                         % print current penalty factor
    while (change > 0.01), loop = loop+1;                                  % start optmization process
        %_____________________________________FINITE-VOLUME THEORY ANAHSIS
        Mat = MatInt(penal(i),xPhys);                                      % material interpolation
        E = Mat{1}; dEdx = Mat{2};
        sK = K0(:)*E(:)';                                                  % stiffness interpolation
        K = StiffnessAssemblage(sK);                                            % assemblage of stiffness matrix
        K = (K+K')/2;                                                      % correct numeric rounding errors to certify the matrix symmetry
        U(fdof) = K(fdof,fdof)\F(fdof);                                    % compute global displacements
        %________________________________________COMPLIENCE AND SENSITIVITY
        fe = reshape(sum((U(dof)*K0).*U(dof),2),nx,ny);
        f  = sum(sum(E.*fe));                                              % objective function: compliance
        dfdx = -dEdx.*fe;                                                  % sensitivity
        %_____________________________________MODIFICATION OF SENSITIVITIES
        sens = sensitivity(x,dfdx,dvdx);
        dfdx(:) = sens{1};
        dvdx(:) = sens{2};
        %_________________UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
        l1 = 0; l2 = 1e5;                                                  % define initial values for bissection method
        while ((l2-l1)/(l1+l2) > 1e-4)
            lmid = (l1+l2)/2;                                              % average lambda
            beta = -dfdx./dvdx/lmid;                                       % density correction factor
            xnew = max(0,max(x-move,min(1,min(x+move,x.*beta.^eta))));     % update of design variables
            xPhys(:) = density(xnew);
            if (mean(xPhys(:)) > volfrac), l1 = lmid; else, l2 = lmid; end 
        end
        xPhys = reshape(xPhys,size(x,1),size(x,2));
        change = max(abs(xnew(:)-x(:)));
        x = xnew;                                                        
        %_____________________________________________________PRINT RESULTS
        fprintf('It: %i\tObjec.: %1.4f\tVol.: %1.3f\tChange: %1.3f\n',...
            loop,f,mean(xPhys(:)),change);
    end
    %_________________________________________PRINT RESULTS AND PLOT DESIGN
    set(fig,'FaceColor','flat','CData',1-xPhys(:)); drawnow
end
t = toc; clearvars -except t;
PrintTime(t);
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
%____________________________________________________________PREPARE FILTER
function [H,Hs] = Filtering(l,h,nx,ny,neig)
neig = round(neig,0);                                                      % Enforce that neig variable to be an integer                                                                               
rmin = neig*sqrt(l^2+h^2);                                                 % radius of the filter                       
iH = ones(nx*ny*(2*(neig-1)+1)^2,1);
jH = ones(size(iH));                    
if (neig ~= 0)
    sH = zeros(size(iH));
    k = 0;
    for j1 = 1:ny
        for i1 = 1:nx
            e1 = (j1-1)*nx+i1;
            for j2 = max(j1-neig,1):min(j1+neig,ny)
                for i2 = max(i1-neig,1):min(i1+neig,nx)
                    e2 = (j2-1)*nx+i2;
                    k = k+1;
                    iH(k) = e1; jH(k) = e2;
                    dx = (i1-i2)*l; dy = (j1-j2)*h;
                    sH(k) = max(0,rmin-sqrt(dx^2+dy^2));
                end
            end
        end
    end
else
    sH = ones(size(iH));
end
H = sparse(iH,jH,sH); Hs = sum(H,2);
%_________________________________________HANDLE FOR PLOTTING DESIGN DOMAIN
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
%________________________________________________________PRINT PROCESS TIME
function PrintTime(t)
H = fix(t/3600); M = fix((t-3600*H)/60); S = fix((t-3600*H-60*M));
disp(['Timing of iterative process: ' num2str(H) 'h ' num2str(M) 'min '...
    num2str(S) 's']);
