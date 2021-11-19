function mout = moist_neq1_solver(params, kplot, qnk)
%
% mout = moist_neq1_solver(params, kplot, qnk)
%
% [Return] mout: a struct containing dispersion and wave field. 
% [Input]  params: model parameters (details below)
%          kplot: zonal wavenumber for display
%          qnk: zonal wavenumbers
%
% params: parameters for the moist shallow water model, default values are:
%     params.alpha = 0.35;        % moisture relaxation time scale, ~ 1 day
%     params.cri = 0.17;          % cri (r) = 0.17, cloud-radiative feedback
%     params.Q0 = 0.22;           % Q0: linear decrease of dQ/dy, Q0 = 0.22 
                                  % gives the MJO, use small values, 0 or -0.05 gives horizontal
                                  % tilt and poleward propagation, hence the BSIBSO
%     params.Lambda = -0.03;      % WISHE + Zonal moisture gradient
%     params.bwishe = 0;          % WISHE for temperature. Its effect is small. Just let it be 0
%     params.Gamma = 0.16;        % Gross moist stability
%     params.bdamping = 0.01;     % Newtonian damping, ~ 30 days
%     params.udamping = 0.03;     % Raleigh friction,  ~ 10 days
% kplot: the k-th zonal wave field to be calculated, default: kplot=1; plot zonal wave 1
% qnk: zonal wavenubmer, default: qnk = -5:0.5:5
%
% Example: Using default values and return wave 1 field: 
% mout = moist_neq1_solver_4github();
% % plot the wave field:
% figure; clf
% quiver(mout.x(1:3:end),mout.y(1:8:end), real(mout.U(1:8:end,1:3:end)), real(mout.V(1:8:end,1:3:end)), 'color', [1 1 1]/2);
% hold on;
% [cc,hh]=contour(mout.x, mout.y, real(mout.phi));
% maxw = max(real(mout.W(:)));
% contour(mout.x, mout.y, real(mout.W), [-maxw, maxw ]*1/4, 'k--');
%
% %  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    July 2021
%
%  Author:
%
%    S. Wang (wangsg at outlook.com)
%
% Comments: 
% Solve the normal modes of the first baroclinic mode moist shallow water 
% model for the MJO/BSISO.
%
%  Reference:
%
% Wang, S., A. H. Sobel. A Unified moisture mode theory of the MJO and
% BSISO. 2021. JC

%% input
if nargin < 1
    params.ptitle='MJO';
    %params.Q = 0;        % Q ~ Qy*(alpha+ep)
    params.cri = 0.17;    % ep = 0.35*0.17; %0.07; %0.06;
    params.alpha = 0.35;
    params.Q0 = 0.22; %*(params.ep+params.alpha);
    params.bwishe = 0;
    params.Lambda = -0.03;
    params.Gamma = 0.16;
    params.bdamping = 0.01;
    params.udamping = 0.03;
end

ptitle = params.ptitle;
Q0 = params.Q0;
cri = params.cri;
alpha = params.alpha;
bwishe = params.bwishe;
Lambda = params.Lambda;
Gamma = params.Gamma;
bdamping = params.bdamping;
udamping = params.udamping;
bwishe = params.bwishe;
 
disp(params)

if nargin < 2
    kplot = 1;
end

if nargin < 3
    qnk = -5:0.5:5;
end

%----------------------------------------------------------
% combine parameters
ep = cri*params.alpha;           % cloud-radiative feedback
f0 = Gamma*alpha + ep*(Gamma-1); % alpha*effective gross moist stability: alpha*[Gamma*(1+r)-r]
Q = Q0*(ep+alpha);               % Q ~ Q0*alpha*(1+r) = Q0*(alpha+ep)
Lambda2 = Lambda*(alpha+ep);     % WISHE parameter*(alpha+ep)

ai=1j;

% model constant for non-dimensionalization
Req = 2*pi*6400e3;       % circumference of earth at equation, with radius Re = 6400e3 km
c = 50;                  % phase speed of first baroclinic mode, 50 m/s
beta = 2.2e-11;          % beta, s^{-1} 

qn = (qnk)*sqrt(c/beta)/Req*(2*pi);

Nk = length(qnk);
sigg = nan(Nk,8); % contain omega, 8 solutions for 8th order polynomial equation
xi1_arr = sigg;

for j = 1:Nk
    
    syms x  % x is omega (frequency/growth rate)
    k = qn(j);
   
    % 2*a1^2  + a1*d0 -d0^2 - 9*d2*a0 = 0 % dispersion equation
    a0 = (f0-ai*x)*(x+ai*udamping);
    a1 = ai*(Lambda2+bwishe*(alpha-ai*x)) + (x+ai*udamping)*Q;
    d0 = ai*(Lambda2+bwishe*(alpha-ai*x)) + (x+ai*udamping)*Q - k*(f0-ai*x) - k^2*(f0-ai*x)*x  - (ai*(x+ai*udamping)*(x+ai*bdamping)*(x+ai*alpha) - (Lambda2+bwishe*(alpha-ai*x))*k*ai)*x;
    d2 = ai*(x+ai*alpha)*(x+ai*bdamping)+k*Q;
    
    verstr = version;
    i0 = findstr(verstr, 'R');
    versionyear = str2num(verstr(i0+1:i0+4));
    if versionyear <= 2017 % earlier version has no str2sym, just convert it a string
        eqn1 = sym([char(2*a1^2) ' - ' char(d0^2) ' - ' char(9*d2*a0) ' + ' char(a1*d0) ]); % dispersion equation
    else
        eqn1 = 2*a1^2 - d0^2 - 9*d2*a0 + a1*d0; % dispersion equation
    end
    
    if 1 == 11 %versionyear <= 2017
        S=solve(eqn1,x);  % solve omega from the dispersion equation, dimensions: number of zonal wavenumbers, number of roots
        aout = vpa(S);
    else
        %tic
        cout = coeffs(eqn1, x, 'All'); % unpack the coefficients of the dispersion equation and use "roots". This is 10 times faster than solve.
        cout1 = formula(cout);  
        %aout = roots(vpa(cout1(end:-1:1), 50)); 
        aout = roots(vpa(cout1, 50)); 
        %toc
    end
    
    sigg(j,1:length(aout)) = aout;  
    
    for n = 1:size(sigg,2)
        xx = sigg(j,n);       
        x = xx;  % copy from above
        % 2*a1^2  + a1*d0 -d0^2 - 9*d2*a0 = 0 % dispersion equation
        a0 = (f0-ai*x)*(x+ai*udamping);
        a1 = ai*(Lambda2+bwishe*(alpha-ai*x)) + (x+ai*udamping)*Q;
        d0 = ai*(Lambda2+bwishe*(alpha-ai*x)) + (x+ai*udamping)*Q  - k*(f0-ai*x) - k^2*(f0-ai*x)*x - (ai*(x+ai*udamping)*(x+ai*bdamping)*(x+ai*alpha) - (Lambda2+bwishe*(alpha-ai*x))*k*ai)*x;
        d2 = ai*(x+ai*alpha)*(x+ai*bdamping)+k*Q;
        res = 2*a1^2  + a1*d0 -d0^2 - 9*d2*a0;
        xi1 = (a1+d0)./a0/6;                    % shape parameter
        
        if real(xi1) <= 0 || real(xx) < 0 
        % make sure RE(xi)>0 and frequency is positive (negative frequencies are redundant).
            sigg(j,n) = nan;
            xi1 = nan;
        end
        if abs(res)>1e-4  
            % sanity check, make sure dispersion equation is closed
            % if its too big, the solver fails
            sigg(j,n) = nan;
            warning('res > 1e-5')
            disp([j,n])
            xi1 = nan;
        end        
        xi1_arr(j,n) = xi1;        
    end
end

sigsort = sigg*complex(nan, nan);
xi1sort = xi1_arr*complex(nan, nan);
for j = 1:Nk
    %[~, ii] = sort(abs(imag(sigg(j,:))),'descend');
    [~, ii] = sort(imag(sigg(j,:)),'descend');
    in = 1;
    for ij = 1:size(sigg,2)
        if ~isnan(sigg(j,ii(ij)))
            sigsort(j,in)= sigg(j,ii(ij));
            xi1sort(j,in)= xi1_arr(j,ii(ij));
            in = in +1;
        end
    end
end
Tscale = 1/sqrt(c*beta)/86400; % days
 
% extract intraseasonal frequency with phase speed between -25 to 25 m/s
ig = nan(Nk,1); omg = ig; cm=ig;
for m =1:Nk
    ss = (real(sigsort(m,:))./qn(m)'*c);        % phase speed
    ii = find((ss>=-25)&(ss<25) | qnk(m) == 0); % phase speed between -25 to 25 m/s and k=0 
    if length(ii) >= 1
        ig(m) = ii(1);
        omg(m) = sigsort(m,ii(1));
        cm(m) = ss(ii(1));
    end
end
omg(isnan(omg)) = complex(nan,nan);

% ------------------------------------------------------------------------------------------------------------------------------------------------------------------------

%%% Now we get the detailed solution, and check the budget for u, v, phi, qv ...
%% Get omega for the most unstable mode and \xi_1. Also set up the grid

[~,j]=min(abs(qnk-kplot)); % find the nearest zonal wavenumber to plot
if abs(qnk(j)-kplot) > 1e-3
    warning(['No zonal wave ' num2str(kplot) '  Use the nearest: ', num2str(qnk(j))]);
end
%i = ai;

instb = 1; % n-th unstable mode
sig_k_sort = sigsort(j, ~isnan(sigsort(j,:)));
xx = sig_k_sort(instb);
omg_save = omg; 

omg(j) = xx;

x = xx;    % copy from above
k = qn(j); % need qn here, not qnk

a0 = (f0-ai*x)*(x+ai*udamping);
a1 = ai*(Lambda2+bwishe*(alpha-ai*x)) + (x+ai*udamping)*Q;
d0 = ai*(Lambda2+bwishe*(alpha-ai*x)) + (x+ai*udamping)*Q - k*(f0-ai*x) - k^2*(f0-ai*x)*x - (ai*(x+ai*udamping)*(x+ai*bdamping)*(x+ai*alpha) - (Lambda2+bwishe*(alpha-ai*x))*k*ai)*x;
d2 = ai*(x+ai*alpha)*(x+ai*bdamping)+k*Q;
res = 2*a1^2  + a1*d0 -d0^2 - 9*d2*a0;

if abs(res) > 1e-5
    warning('Something wrong. Dispersion not closed!')
end

xi1 = (a1+d0)/6/a0; % shape parameter
            
Lx = 2*pi/abs(qn(j));
if abs(qn(j))<1e-10
    Lx = 30;
end

x = linspace(0, 2*Lx, 130); 
y = linspace(-4.5,4.5,180);
Nx = length(x); 
Ny = length(y);

[X,Y] = meshgrid(x,y);

% calculate the variables
V = Y.*exp(-xi1.*Y.^2).*exp(ai.*k*X); % not cos(k*X);
v_y = (1-2*xi1*Y.^2).*exp(-xi1.*Y.^2).*exp(ai.*k*X);
v_yy = (-2*xi1*Y.*(1-2*xi1*Y.^2) - 4*xi1*Y).*exp(-xi1.*Y.^2).*exp(ai.*k*X);

G = f0;
L = Lambda2 +bwishe*(alpha-ai*xx);
A = (xx+ai*udamping)*(xx+ai*bdamping)*(xx+ ai*alpha) - k.^2*xx - k*L -ai*G*k.^2; 
% A is actually the v=0 dispersion expression

ytmp = (ai*k*xx - k*G + ai*L - 2*xi1.*(xx+ai*udamping).*(G-ai*xx));
phi = 1/A* ((xx+ai*udamping)*(G-ai*xx) +ytmp.*Y.^2).*exp(-xi1.*Y.^2).*exp(ai*k*X);
phi = phi + (xx+ai*udamping).*Q*Y.*V/A;
 
phi_y = 1/A* ((xx+ai*udamping)*(G-ai*xx)*(-2*xi1.*Y) + ytmp.*(2*Y - 2*xi1*Y.^3)).*exp(-xi1.*Y.^2).*exp(ai*k*X);
phi_y = phi_y + 1/A *((xx+ai*udamping)*Q.*(V + Y.*v_y)); % dPhi/dy

U = ai/(xx+ai*udamping).*Y.*V + k/(xx+ai*udamping).*phi;
W = L/(G-ai*xx).*U - (xx+ai*alpha)*(xx+ai*bdamping)/(G-ai*xx).*phi;
W = W + Q.*Y.*V./(G-ai*xx);

% checking budget of each equation below. This is important actually. 

% u equation: du/dt - yv = -dphi/dx - udamping*u
ueq1 = -ai*omg(j)*U + udamping*U;
ueq2 = -Y.*V;
ueq3 = ai*qn(j)*phi;
u_res = ueq1+ueq2+ueq3;

% v equation: dv/dt + yv = -dphi/dy
veq1 = -ai*omg(j)*V;
veq2 = Y.*U;
%phi_y2(2:end-1,:) = (phi(3:end,:)-phi(1:end-2,:))./(y(3)-y(1)); % finite difference approx.
veq3 = phi_y;
v_res = veq1 + veq2 + veq3;

% residual for the single V equation
vsingle_res = a0*v_yy + a1.*Y.*v_y + (d0 + d2*Y.^2).*V;

% Qv equation: dqv/dt + (Gamma-1)*W = Lambda*U + v*y*qv
qv = (+ai*omg(j)*phi + W - bwishe*U- bdamping*phi)/(alpha + ep) ; % from the phi equation
qeq1 = -ai*omg(j)*qv;
qeq2 = (Gamma-1)*W;
qeq3 = -Lambda*U;
qeq4 = alpha*qv;
%qeq5 = -Y.*V*Q;
qeq5 = -Y.*V*Q/(alpha+ep);
q_res = qeq1+qeq2+qeq3+qeq4+qeq5;

% phi or b equation: dphi/dt - w = -alpha*(1+r)*q - bdamping*phi + bwishe*U
dbeq1 =  -ai*omg(j)*(-phi) + bdamping*(-phi) -bwishe*U;
dbeq2 = W;
dbeq3 = -(alpha+ep)*qv;
b_res = dbeq1+dbeq2 + dbeq3;

% equation for divergence: dudx + dvdy + w = 0
deq1 = ai*qn(j)*U;
deq2 = v_y;
deq3 = W; 
d_res = deq1 + deq2 + deq3; 

% moist static energy h = qv - phi
heq1 = -ai*omg(j)*(qv-phi);
heq2 = Gamma*W;
heq3 = -Lambda*U - bwishe*U;
heq4 = -ep*qv;
heq5 = -Y.*V*Q/(alpha+ep);
heq6 = -bdamping*phi;
h_res = heq1 + heq2 + heq3 + heq4 + heq5 +heq6;
mse = qv-phi;


if max(abs(q_res(:))) > 1e-5
    warning('Qv budget for wave 1 not closed')
end

if max(abs(u_res(:))) > 1e-5
    warning('U budget for wave 1 not closed')
end

if max(abs(v_res(:))) > 1e-5
    warning('V budget for wave 1 not closed')
end

if max(abs(b_res(:))) > 1e-5
    warning('B budget for wave 1 not closed')
end

if max(abs(d_res(:))) > 1e-5
    warning('DIV budget for wave 1 not closed')
end

if max(abs(h_res(:))) > 1e-5
    warning('MSE budget for wave 1 not closed')
end

if max(abs(vsingle_res(:))) > 1e-5
    vrelerr = max(abs(vsingle_res(:)./V(:)));
    warning(['Single V budget for wave ' num2str(kplot) ' not closed, relative error is: ' num2str(vrelerr)])
end

disp(['Residuals for q,u,v,d, b, mse, v1: ',  ...
     num2str([max(abs(q_res(:))), max(abs(u_res(:))), max(abs(v_res(:))), ...
              max(abs(d_res(:))), max(abs(b_res(:))), max(abs(h_res(:))), max(abs(vsingle_res(:)))])]);
disp('done ... ');
          

%% output struct
mout.params = params;
mout.sigg = sigg;
mout.sigsort = sigsort;
mout.Tscale = Tscale;
mout.omg = omg;
mout.cm = cm;
mout.qnk = qnk;
mout.Req = Req;
mout.c = c;
mout.beta = beta;
mout.x = x;
mout.y = y;
mout.U = U;
mout.V = V;
mout.phi = phi;
mout.qv = qv;
mout.W = W;
mout.mse = mse;
mout.xi1 = xi1;
mout.xi1sort = xi1sort;

%% plot the dispersion relationship below

% default color in python: https://matplotlib.org/stable/users/dflt_style_changes.html
c_py = [31	119	180;
    255	127	14;
    44	160	44;
    214	39	40;
    148	103	189;
    140	86	75;
    227	119	194;
    127	127	127	;
    188	189	34;
    23	190	207	]/256;

figure; clf
subplot(221)
for i = 1:size(sigsort,2)
    for kk = 1:length(qnk)
        if imag(sigsort(kk,i)) > 0  % unstable modes
            
            msize =  abs(imag(sigsort(kk,i)))*1080;
            msize = min(max(msize, 5),10);
            
            plot(qnk(kk), real(sigsort(kk,i))/Tscale/2/pi, 'o', ...
                'MarkerEdgeColor',c_py(1,:),...
                'MarkerFaceColor',c_py(1,:),...
                'MarkerSize',msize);
            hold on;
        else                        % stable modes
            
            msize =  -abs(imag(sigsort(kk,i)))*1080;
            msize = -max(min(-msize, -5),-10);
            
            plot(qnk(kk), real(sigsort(kk,i))/Tscale/2/pi, 'o', ...
                'MarkerEdgeColor',c_py(1,:),...
                'MarkerFaceColor','w',...
                'MarkerSize',msize);
            hold on;
        end
    end
end
grid on;


subplot(222)

omg(isnan(omg)) = complex(nan,nan);
plot(qnk, real(omg)/Tscale/2/pi, 'o-')
grid on;
hold on;
%ylim([0,0.5])
title('Frequency(Omega)')
ylabel('1/days')

subplot(223)
plot(qnk, imag(omg)/Tscale, 'o-')
grid on;
title('growth rate')
ylabel('1/days')


subplot(224)
plot(qnk, cm, 'o-')
grid on;
hold on;
ylabel('m/s')
title('phase speed')
 
for i =1:4
    subplot(2,2,i)
    xlim([qnk(1)-0.5 qnk(end)+0.5])
    xlabel('zonal wavenumber')
    set(gca, 'xtick', qnk(1):qnk(end))
end


figure; clf
quiver(x(1:3:end), y(1:8:end), real(U(1:8:end,1:3:end)), real(V(1:8:end,1:3:end)), 'color', [1 1 1]/2);
hold on;
contour(x, y, real(mout.phi));
maxw = max(abs(real(W(:))));
contour(x, y, real(W), [ maxw,  maxw ]*3/4, 'k-');
contour(x, y, real(W), [-maxw, -maxw ]*3/4, 'k--');
title(['k=' num2str(qnk(j))]);
  
