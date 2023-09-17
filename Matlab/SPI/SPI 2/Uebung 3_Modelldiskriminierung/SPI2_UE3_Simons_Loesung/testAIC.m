load iddata1 z1;
np = 2;
sys = tfest(z1,np);
value = aic(sys)