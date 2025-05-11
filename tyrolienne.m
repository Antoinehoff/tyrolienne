
%% LA TYROLIENNE
% Parametres physiques
m     = 0.1e+0;    % la masse de la tyrolienne [m]
l     = 1.0e+0;    % la longueur de la tyrolienne au CM [m]
g     = 9.8e+0;    % la gravite [m2/s]
alpha = 0*pi/20;    % l'angle du cable [rad]
%nu    = 0.5;      % frottement [N/(m/s)]
% Conditions initiales
th0   = pi/4; % angle initial de la barre [rad]
om0   = 0;   % vitesse initiale de la barre [rad/s]
% Proprietes physiques
Ig    = m*l^2;                       % le moment d'inertie de la barre
eta   = Ig/(m*l^2);                  % ratio inertie barre/pendule
w0    = sqrt(g/l);  % frequence propre OH
T     = 2*pi/w0;                     % Periode OH
OH    = @(t) (th0+alpha).*cos(w0*t) + om0/w0 * sin(w0*t) - alpha;    % Solution OH
% Parametres numeriques
Tmax  = 10*T;              % le temps de simulation max [s]
dt    = 1.0e-3;            % le pas de temps [s]
Nstep = ceil(Tmax/dt);     % nombre de pas de temps (arrondit vers le haut)
% pour les plots
h = 2; cable = @(x) -tan(alpha)*x + h;
%------- Solver
% l'equation a resoudre est de la form A_x y^2 + B_x y' + C_x = 0
% et y = x'
A = @(x) cos(alpha+x)*sin(alpha+x)^2;  % terme devant le omega^2
B = @(x) sin(alpha+x).*(sin(alpha+x)^2+eta); % facteur devant le omega dot 
C = @(x) g/l*cos(alpha)*sin(alpha+x)^2; % facteur non differentiel
% Eq diff in the form thetadot = ..., omegadot = ...
tyrol = @(t,y) [y(2);-(A(y(1))*y(2)^2+C(y(1)))/B(y(1))];

% inconnues
yEE   = zeros(2,Nstep); % 1: angle en fonction du temps [rad]
                        % 2: vitesse angulaire en fonction du temps [rad/s]
yEE(1,1) = th0;
yEE(2,1) = om0;

%euler explicit a la con
t = 0;
for it = 1:Nstep-1
    yEE(:,it+1) = yEE(:,it) + dt * tyrol(t,yEE(:,it));
    t = t + dt;
end
tEE = (0:Nstep-1)*dt;

%ODE 45
[tODE,yODE] = ode45(tyrol,[0 Tmax],[th0,om0]);
%ODE 113
% [tODE,yODE] = ode113(tyrol,[0 Tmax],[th0,om0]);

%----- Post processing (to compute the movement of the translation point)
theta = yODE(:,1)';                  %theta
omega = yODE(:,2)';                  %omega
omdot = diff(omega')./diff(tODE); %omega dot, needed in OA acceleration
omdot = [0; omdot]';

Acc   = @(th,w,wdot) (l*w^2*sin(alpha+th) + g*sin(alpha) - l*wdot*cos(alpha+th));

yA = zeros(2,numel(tODE)); % dynamique de OA selon le cable
yA(1,1) = 0.0; % xA(t=0)
yA(2,1) = 0.0; % vA(t=0)

for it = 1:numel(tODE)-1
    dt = tODE(it+1) - tODE(it);
    yA(2,it+1) = yA(2,it) + dt * Acc(theta(it),omega(it),omdot(it));
    yA(1,it+1) = yA(1,it) + dt * yA(2,it);
end

xA = yA(1,:); vA = yA(2,:);

% recover the trajectory of OA from lab referential
OAx = 0.0 + xA.*cos(alpha); % along ex
OAy = h   - xA.*sin(alpha); % along ey
% and its speed
vOAx = 0.0 + vA.*cos(alpha); % along ex
vOAy = 0.0 - vA.*sin(alpha); % along ey

% Trajectoire de OP (bout de la tige)
APx =  2*l*sin(theta);
APy = -2*l*cos(theta);
OPx = OAx + APx;   % fin de la tige
OPy = OAy + APy;   % fin de la tige
OGx = OAx + APx/2; % centre de masse
OGy = OAy + APy/2; % centre de masse
% Et la vitesse
vOPx = vOAx + 2*l*omega.*cos(theta);
vOPy = vOAy + 2*l*omega.*sin(theta);
vOGx = vOAx +   l*omega.*cos(theta);
vOGy = vOAy +   l*omega.*sin(theta);
% Bilan d'énergie, exprimé selon les coordonnées mobiles
Epot = -m*g*(l*cos(theta) + xA*sin(alpha));
Ekin = 0.5*m*(vOGx.^2 + vOGy.^2);
Erot = 0.5*Ig*omega.^2;
Emec = Epot+Ekin+Erot;
%% ----- Plots
if 1
%% Film/trajectory
% Interpolation de la trajectoire pour les films (équidistant en t)
Ax = griddedInterpolant(tODE,OAx); Ay = griddedInterpolant(tODE,OAy);
Px = griddedInterpolant(tODE,OPx); Py = griddedInterpolant(tODE,OPy);
Gx = griddedInterpolant(tODE,OGx); Gy = griddedInterpolant(tODE,OGy);
% window measure
t  = linspace(0,Tmax,numel(tODE));
min_x = min([OPx OAx])-0.5*l; max_x = max([OPx OAx])+0.5*l;
min_y = min([OPy OAy])-0.5*l; max_y = max([OPy OAy])+0.5*l;
figure
for it = 1:1:numel(t)
    t_ = t(it);
    % Le cable
    plot([min_x,max_x],cable([min_x,max_x]),'--k','color',[0.5 0.5 0.5]); hold on;
    % La barre
    plot([Ax(t_),Px(t_)],[Ay(t_),Py(t_)],'k');
    % La fixation
    plot(Ax(t_),Ay(t_),'sk'); hold on;
    % Le centre de masse
    plot(Gx(t_),Gy(t_),'xr'); hold on;
    % Labels et titre
    xlabel('x[m]'); ylabel('y[m]');
    xlim([min_x,max_x]); 
    ylim([min_y,max_y]);
    pbaspect([(max_x-min_x) (max_y-min_y) 1])
    title(['t=',num2str(t(it)),'[s]']);
    drawnow
    if(it<numel(t))
%         pause((t(it+1)-t(it)))
    clf
    end
end
end
if 1
    %%
figure
%--------------angle time evolution
subplot(311)
% solution Euler explicit
plot(tEE,yEE(1,:),'-r'); hold on;
% solution of matlab solver
plot(tODE,yODE(:,1),'-b'); hold on;
% solution Oscillateur harmonique
plot(tEE,OH(tEE),'-k'); 
% Position d'equilibre
plot([tEE(1) tEE(end)], -alpha*[1 1],'--k');
ylabel('[rad]');
set(0,'defaultTextInterpreter','none');
legend('EE','ODE45',...
    '$(\theta_0+\alpha)cos(\omega t)+\omega_0/\omega sin(\omega t) - \alpha$',...
    '$-\alpha$');

% position time evolution
subplot(312)
plot(tODE,OAx,'b-');  hold on;
plot(tODE,OAy,'b--'); hold on;
plot(tODE,OGx,'r-');  hold on;
plot(tODE,OGy,'r--'); hold on;
plot(tODE,sqrt((OAx-OPx).^2+(OAy-OPy).^2),'g-'); hold on;
ylabel('[m]');
legend('OAx','OAx','OGx','OGy','l(t)');

% Energy time evolution
subplot(325)
plot(tODE,Epot-Epot(1),'b-');  hold on;
plot(tODE,Ekin-Ekin(1),'r-'); hold on;
plot(tODE,Erot-Erot(1),'g-');  hold on;
ylabel('[J]'); xlabel('[s]');
legend('Epot-Epot0','Ekin-Ekin0','Erot-Erot0');
subplot(326)
plot(tODE,abs(Emec-Emec(1))/abs(Emec(1))*100,'k-'); hold on;
ylabel('[\%]'); xlabel('[s]');
legend('|Emec-Emec0|/Emec0');
end
%% Storage
% % frequency spectrum
% Freq = 2*pi/dt * [1:floor(numel(tEE)/2)];
% subplot(222)
% AMP = abs(fft(yEE(1,:)));  AMP = AMP(1:floor(numel(AMP)/2));
% loglog(Freq,AMP,'or'); hold on;
% F = griddedInterpolant(tODE,yODE(:,1));
% AMP = abs(fft(F(tEE))); AMP = AMP(1:floor(numel(AMP)/2));
% loglog(Freq,AMP,'ob'); hold on; 
% AMP = abs(fft(OH(tEE))); AMP = AMP(1:floor(numel(AMP)/2));
% loglog(Freq,AMP,'ok'); hold on;



