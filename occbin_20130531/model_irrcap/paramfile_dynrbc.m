BETA=0.96;
ALPHA=0.33;
DELTAK=0.10;
GAMMAC=2;
RHOA = 0.9;
PHII = 0.975;
PHIK = 0.0;
PSI = 0;        % adjustment cost for capital
PSINEG = 0;     % adjustment cost if investment is negative

load PARAM_EXTRA

KSS = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));
ISS = DELTAK*KSS;
