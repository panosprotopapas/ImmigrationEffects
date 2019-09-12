// 1. Declare variables and parameters
// -----------------------------------

var  y ystar a astar r rstar d dstar p pstar h c cstar hstar cn cnstar gdp gdpstar;
varexo e estar z zstar;
parameters ybar ybarstar rho1 rho2 abar abarstar beta dbar dbarstar kappa alpha delta hbar eta gamma w wstar theta;


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
delta	   = 1.00;
hbar       = 1.00;
gamma      = 0.62;
abarstar   = eta*abar;
w          = 0.628125;
wstar      = 0.634885;
theta      = 1-0.92348;

// 3. Declare the model's dynamics
// -----------------------------------

model;

y=ybar+rho1*(y(-1)-ybar)-e;
ystar=ybarstar+rho1*(ystar(-1)-ybarstar)-estar;
a=abar+rho2*(a(-1)-abar)-z;
astar=abarstar+rho2*(astar(-1)-abarstar)-zstar;

r=(1/beta)*(d/dbar)^kappa;
rstar=(1/beta)*(dstar/dbarstar)^kappa;

p=((1-alpha)*c)/(alpha*(a*((h)^gamma)+((delta*wstar*theta*hstar)/(pstar))));
pstar*(alpha*(astar*(hstar)^gamma-(delta*wstar*theta*hstar)/(pstar)))=((1-alpha)*cstar);

(c/(a*h^gamma+(delta*wstar*theta*hstar)/pstar))^(alpha-1)=beta*r*((c(1))/((a(1))*(h(1))^gamma+(delta*(wstar(1))*theta*hstar(1))/(pstar(1))))^(alpha-1);

(cstar/(astar*hstar^gamma-(delta*wstar*theta*hstar)/pstar))^(alpha-1)=beta*rstar*((cstar(1))/(((astar(1))*hstar(1)^gamma-(delta*(wstar(1))*theta*hstar(1)))/(pstar(1))))^(alpha-1);

c+d(-1)=(1-delta)*wstar*theta*hstar+y+d/r;
cstar+dstar(-1)=(delta-1)*wstar*theta*hstar+ystar+dstar/r;

h=(w/(gamma*a*p))^(1/(gamma-1));
hstar=(wstar/(gamma*astar*pstar))^(1/(gamma-1));

cn=a*h^gamma+(delta*wstar*theta*hstar)/pstar;
cnstar=astar*(1-theta)*hstar-(delta*wstar*theta*hstar)/pstar;

gdp=y+a*h^gamma;
gdpstar=ystar+astar*hstar^gamma;

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

h=0.92348;
c=0.983646;
cstar=1.05378;
p=0.982917;
pstar=0.993496;
hstar=2-0.92348;

cn=abar*h^gamma+(delta*wstar*theta*hstar)/pstar;
cnstar=abarstar*(1-theta)*hstar-(delta*wstar*theta*hstar)/pstar;

gdp=ybar+abar*h^gamma;
gdpstar=ybarstar+abarstar*hstar^gamma;

end;

steady;

// 5. Simulate the stochastic model
// ----------------------------------------

shocks;
var e=0.058;

end;


check;
stoch_simul(order=1,irf=100, relative_irf);