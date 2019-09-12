// 1. Declare variables and parameters
// -----------------------------------

var  y ystar a astar r rstar d dstar p pstar h c cstar w wstar hstar cn cnstar gdp gdpstar real realstar;
varexo e estar z zstar;
parameters ybar ybarstar rho1 rho2 abar abarstar beta dbar dbarstar kappa alpha delta hbar eta gamma ;


// 2. Calibrate parameters values
// -----------------------------------



eta        = 1.06;
ybar   	   = 1.00;
ybarstar   = eta*ybar;
rho1       = 0.95;
rho2       = 0.95;
abar	   = 1.00;
beta	   = 0.9926;
dbar       = 2.21;
dbarstar   = 0.84;
kappa      = 0.0000000001;
alpha      = 0.50;
delta	   = 0.50;
hbar       = 1.00;
gamma      = 0.62;
abarstar   = eta*abar;

// 3. Declare the model's dynamics
// -----------------------------------

model;

y=ybar+rho1*(y(-1)-ybar)-e;
ystar=ybarstar+rho1*(ystar(-1)-ybarstar)-estar;
a=abar+rho2*(a(-1)-abar)-z;
astar=abarstar+rho2*(astar(-1)-abarstar)-zstar;

r=(1/beta)*(d/dbar)^kappa;
rstar=(1/beta)*(dstar/dbarstar)^kappa;

p=((1-alpha)*c)/(alpha*(a*((h)^gamma)+((delta*wstar*(hbar-h))/(pstar))));
pstar=((1-alpha)*cstar)/(alpha*(astar*(2*hbar-h)^gamma-(delta*wstar*(hbar-h))/(pstar)));

(c/(a*h^gamma+(delta*wstar*(hbar-h))/pstar))^(alpha-1)=beta*r*((c(1))/((a(1))*(h(1))^gamma+(delta*(wstar(1))*(hbar-(h(1))))/(pstar(1))))^(alpha-1);
(cstar/(astar*(2*hbar-h)^gamma-(delta*wstar*(hbar-h))/pstar))^(alpha-1)=beta*rstar*((cstar(1))/((astar(1))*(2*hbar-h(1))^gamma-(delta*(wstar(1))*(hbar-(h(1))))/(pstar(1))))^(alpha-1);

w/wstar=(delta*p)/pstar+(1-delta);

h=(w/(gamma*a*p))^(1/(gamma-1));
2*hbar-h=(wstar/(gamma*astar*pstar))^(1/(gamma-1));


c+d(-1)=(1-delta)*wstar*(hbar-h)+y+d/r;
cstar+dstar(-1)=(delta-1)*wstar*(hbar-h)+ystar+dstar/r;


hstar=2*hbar-h;
cn=a*h^gamma+(delta*wstar*(hbar-h))/pstar;
cnstar=astar*(2*hbar-h)^gamma-(delta*wstar*(hbar-h))/pstar;

gdp=y+a*h^gamma;
gdpstar=ystar+astar*hstar^gamma;

real=w/p;
realstar=wstar/pstar;

end;

// 4. Define initial (steady-state) values
// ----------------------------------------

initval;

a=abar;
astar=abarstar;
y=ybar;
ystar=ybarstar;
d=dbar;
dstar=dbarstar;

r=(1/beta);
rstar=r;

h=hbar*((((ybar-(1-beta)*dbar)*(2*alpha+gamma*(1-alpha)*(ybarstar-(1-beta)*dbarstar))))/((ybar+ybarstar-(1-beta)*(dbar+dbarstar))*(alpha+gamma*(1-alpha))));

c=((2*alpha+gamma*(1-alpha))*(ybar-(1-beta)*dbar)+gamma*(1-alpha)*(ybarstar-(1-beta)*dbarstar))/(2*(alpha+gamma*(1-alpha)));
cstar=((2*alpha+gamma*(1-alpha))*(ybarstar-(1-beta)*dbarstar)+gamma*(1-alpha)*(ybar-(1-beta)*dbar))/(2*(alpha+gamma*(1-alpha)));

p=((1-alpha)/(2*alpha*abar*hbar^gamma))*((((2*alpha+gamma*(1-alpha))*(ybar-(1-beta)*dbar)+gamma*(1-alpha)*(ybarstar-(1-beta)*dbarstar)))/((alpha+gamma*(1-alpha))*(ybar+ybarstar-(1-beta)*(dbar+dbarstar))^(gamma/(gamma-1))))^(1-gamma);
pstar=((1-alpha)/(2*alpha*abar*hbar^gamma))*((((2*alpha+gamma*(1-alpha))*(ybarstar-(1-beta)*dbarstar)+gamma*(1-alpha)*(ybar-(1-beta)*dbar)))/((alpha+gamma*(1-alpha))*(ybar+ybarstar-(1-beta)*(dbar+dbarstar))^(gamma/(gamma-1))))^(1-gamma);


w=(gamma*(1-alpha)*(ybar+ybarstar-(1-beta)*(dbar+dbarstar)))/(2*alpha*hbar);
wstar=w;

hstar=2*hbar-h; 

cn=abar*h^gamma+(delta*wstar*(hbar-h))/pstar;

cnstar=abarstar*(2*hbar-h)^gamma-(delta*wstar*(hbar-h))/pstar;

gdp=ybar+abar*h^gamma;
gdpstar=ybarstar+abarstar*hstar^gamma;

real=w/p;
realstar=wstar/pstar;

end;

steady;

// 5. Simulate the stochastic model
// ----------------------------------------

shocks;
var e=0.06;

end;


check;
stoch_simul(order=1,irf=100, relative_irf);