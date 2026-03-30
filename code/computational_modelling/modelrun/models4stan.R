mstr1 = "
data{
int<lower=1>  nTrial;
int<lower=1>  nSub;
int<lower=1>  nStim;
int<lower=1>  nResp;
int<lower=1>  K;
int<lower=1>  s[nSub, nTrial];
int<lower=1>  a[nSub, nTrial];
int<lower=1>  ya[nSub, nTrial];
int<lower=-1> r[nSub, nTrial];
int<lower=1>  rew[nSub, nTrial];
int<lower=1>  yax[nTrial*nSub];
int<lower=1>  nData;

matrix[nResp,nStim]   Qi;
vector[nStim]         Vi;
} 


parameters{
vector[K] X;
vector<lower=0, upper=20>[K] sdX;
vector<lower=-10,  upper=10>[nSub] x1;
vector<lower=-8,  upper=8>[nSub] x2;
} 

transformed parameters{
real rho[nSub];
real<lower=0, upper=1> epsilon[nSub];

matrix[nResp,nStim] Q;
vector[nResp]       q;
simplex[nResp]      p0;
real                er;
vector[nResp] BGx[nTrial*nSub];
vector[nStim]       valenced;

for(iSub in 1:nSub){
rho[iSub]          <- exp(x1[iSub]);
epsilon[iSub]      <- inv_logit(x2[iSub]);
}

for(iSub in 1:nSub){
Q            <- Qi * rho[iSub];
for (iStim in 1:nStim) {
valenced[iStim] <- 0;
}

for(iTrial in 1:nTrial){
for(iResp in 1:nResp){
q[iResp] <- valenced[s[iSub,iTrial]] * Q[iResp,s[iSub,iTrial]]; 
}

p0         <- softmax(q); 

for(iResp in 1:nResp){
BGx[(iSub-1)*nTrial+iTrial,iResp] <- p0[iResp];
} 

er         <- rho[iSub] * r[iSub,iTrial];  
Q[ya[iSub,iTrial],s[iSub,iTrial]] <- Q[ya[iSub,iTrial],s[iSub,iTrial]] + epsilon[iSub] * (er - Q[ya[iSub,iTrial],s[iSub,iTrial]]);


if (r[iSub,iTrial] != 0) 
valenced[s[iSub,iTrial]] <- 1;

} 
} 

} 

model{
vector[nResp] theta;

X[1]  ~ normal(2, 3);
X[2]  ~ normal(0, 2);
sdX   ~ cauchy(0, 2); 


x1 ~ normal(X[1], sdX[1]);
x2 ~ normal(X[2], sdX[2]);


for(iData in 1:nData){
yax[iData] ~ categorical(BGx[iData]);
}
} 

generated quantities {

vector[nData] log_lik; 

for(iData in 1:nData){
log_lik[iData] <- categorical_log(yax[iData], BGx[iData]);
} 

} 
" 

mstr2 = "
data{
int<lower=1>  nTrial;
int<lower=1>  nSub;
int<lower=1>  nStim;
int<lower=1>  nResp;
int<lower=1>  K;
int<lower=1>  s[nSub, nTrial];
int<lower=1>  a[nSub, nTrial];
int<lower=1>  ya[nSub, nTrial];
int<lower=-1> r[nSub, nTrial];
int<lower=1>  rew[nSub, nTrial];
int<lower=1>  yax[nTrial*nSub];
int<lower=1>  nData;

matrix[nResp,nStim]   Qi;
vector[nStim]         Vi;

} 



parameters{


vector[K] X;
vector<lower=0, upper=20>[K] sdX;


vector<lower=-10,  upper=10>[nSub] x1;
vector<lower=-8,  upper=8>[nSub] x2;
vector<lower=-8,  upper=8>[nSub] x3;
} 



transformed parameters{

real rho[nSub];
real<lower=0, upper=1> epsilon[nSub];
real gobias[nSub];

matrix[nResp,nStim] Q;
vector[nResp]       q;
simplex[nResp]      p0;
real                er;
vector[nResp] BGx[nTrial*nSub];
vector[nStim]       valenced;

for(iSub in 1:nSub){
rho[iSub]          <- exp(x1[iSub]);
epsilon[iSub]      <- inv_logit(x2[iSub]);
gobias[iSub]       <- x3[iSub];
}

for(iSub in 1:nSub){
Q            <- Qi * rho[iSub];
for (iStim in 1:nStim) {
valenced[iStim] <- 0;
}

for(iTrial in 1:nTrial){
for(iResp in 1:nResp){
q[iResp] <- valenced[s[iSub,iTrial]] * Q[iResp,s[iSub,iTrial]]; 
}
for(iResp in 1:2){
q[iResp] <- q[iResp] + gobias[iSub]; 
}
p0         <- softmax(q); 

for(iResp in 1:nResp){
BGx[(iSub-1)*nTrial+iTrial,iResp] <- p0[iResp];
} 

er         <- rho[iSub] * r[iSub,iTrial];  
Q[ya[iSub,iTrial],s[iSub,iTrial]] <- Q[ya[iSub,iTrial],s[iSub,iTrial]] + epsilon[iSub] * (er - Q[ya[iSub,iTrial],s[iSub,iTrial]]);


if (r[iSub,iTrial] != 0) 
valenced[s[iSub,iTrial]] <- 1;

} 
} 

} 



model{
vector[nResp] theta;


X[1]  ~ normal(2, 3);
X[2]  ~ normal(0, 2);
X[3]  ~ normal(0, 3);
sdX   ~ cauchy(0, 2); 


x1 ~ normal(X[1], sdX[1]);
x2 ~ normal(X[2], sdX[2]);
x3 ~ normal(X[3], sdX[3]);



for(iData in 1:nData){
yax[iData] ~ categorical(BGx[iData]);
}
} 



generated quantities {

vector[nData] log_lik; 

for(iData in 1:nData){
log_lik[iData] <- categorical_log(yax[iData], BGx[iData]);
} 

} 

"



# ========================================================================================.
# mstr3 corresponds to M3a in the paper. M3a includes a learning rate (epsilon), 
# feedback sensitivity (rho), go bias, and pavlovian bias (pi).
# ========================================================================================.
mstr3 = "

data{
int<lower=1>  nTrial;
int<lower=1>  nSub;
int<lower=1>  nStim;
int<lower=1>  nResp;
int<lower=1>  K;
int<lower=1>  s[nSub, nTrial];
int<lower=1>  a[nSub, nTrial];
int<lower=1>  ya[nSub, nTrial];
int<lower=-1> r[nSub, nTrial];
int<lower=1>  rew[nSub, nTrial];
int<lower=1>  yax[nTrial*nSub];
int<lower=1>  nData;

matrix[nResp,nStim]   Qi;
vector[nStim]         Vi;

} 



parameters{


vector[K] X;
vector<lower=0, upper=20>[K] sdX;


vector<lower=-10,  upper=10>[nSub] x1; 
vector<lower=-8,  upper=8>[nSub] x2; 
vector<lower=-8,  upper=8>[nSub] x3; 
vector<lower=-8,  upper=8>[nSub] x4; 

} 



transformed parameters{

real rho[nSub];
real<lower=0, upper=1> epsilon[nSub];
real gobias[nSub];
real pibias[nSub];

matrix[nResp,nStim] Q;
vector[nResp]       q;
vector[nStim]       V;
simplex[nResp]      p0;
real                er;
vector[nResp]       BGx[nTrial*nSub];
vector[nStim]       valenced;

for(iSub in 1:nSub){

rho[iSub]          <- exp(x1[iSub]);
epsilon[iSub]      <- inv_logit(x2[iSub]);
gobias[iSub]       <- x3[iSub];
pibias[iSub]       <- x4[iSub];
}

for(iSub in 1:nSub){

Q            <- Qi * rho[iSub];
V            <- Vi;
for (iStim in 1:nStim) {
valenced[iStim] <- 0;
}

for(iTrial in 1:nTrial){

for(iResp in 1:nResp){
q[iResp] <- valenced[s[iSub,iTrial]] * Q[iResp,s[iSub,iTrial]]; 
}
for(iResp in 1:2){
q[iResp] <- q[iResp] + gobias[iSub] + pibias[iSub] * valenced[s[iSub,iTrial]] * V[s[iSub,iTrial]]; 
}

p0         <- softmax(q); 
for(iResp in 1:nResp){
BGx[(iSub-1)*nTrial+iTrial,iResp] <- p0[iResp];
} 


er         <- rho[iSub] * r[iSub,iTrial]; 
Q[ya[iSub,iTrial],s[iSub,iTrial]] <- Q[ya[iSub,iTrial],s[iSub,iTrial]] + epsilon[iSub] * (er - Q[ya[iSub,iTrial],s[iSub,iTrial]]);


if (r[iSub,iTrial] != 0) 
valenced[s[iSub,iTrial]] <- 1;

} 
} 

} 



model{
vector[nResp] theta;


X[1]  ~ normal(2, 3); 
X[2]  ~ normal(0, 2); 
X[3]  ~ normal(0, 3); 
X[4]  ~ normal(0, 3);
sdX   ~ cauchy(0, 2); 


x1 ~ normal(X[1], sdX[1]);
x2 ~ normal(X[2], sdX[2]);
x3 ~ normal(X[3], sdX[3]);
x4 ~ normal(X[4], sdX[4]);


for(iData in 1:nData){
yax[iData] ~ categorical(BGx[iData]);
}
} 



generated quantities {

vector[nData] log_lik; 

for(iData in 1:nData){
log_lik[iData] <- categorical_log(yax[iData], BGx[iData]);
} 

} 

" 



# ========================================================================================.
# mstr4 corresponds to M3b in the paper. M3b includes a learning rate (epsilon), 
# feedback sensitivity (rho), go bias, learning bias (kappa), but not a pavlovian bias (pi).
# ========================================================================================.
mstr4 = "

data{
int<lower=1>  nTrial;
int<lower=1>  nSub;
int<lower=1>  nStim;
int<lower=1>  nResp;
int<lower=1>  K;
int<lower=1>  s[nSub, nTrial];
int<lower=1>  a[nSub, nTrial];
int<lower=1>  ya[nSub, nTrial];
int<lower=-1> r[nSub, nTrial];
int<lower=1>  rew[nSub, nTrial];
int<lower=1>  yax[nTrial*nSub];
int<lower=1>  nData;

matrix[nResp,nStim]   Qi;
vector[nStim]         Vi;

} 



parameters{


vector[K] X;
vector<lower=0, upper=20>[K] sdX;


vector<lower=-10,  upper=10>[nSub] x1; 
vector<lower=-8,  upper=8>[nSub] x2; 
vector<lower=-8,  upper=8>[nSub] x3; 
vector<lower=-8,  upper=8>[nSub] x4; 

} 



transformed parameters{

real rho[nSub];
real<lower=0, upper=1> epsilon[nSub];
real gobias[nSub];
real pibias[nSub];
vector[2] biaseps[nSub];

matrix[nResp,nStim] Q;
vector[nResp]       q;
simplex[nResp]      p0;
real                er;
vector[nResp]       BGx[nTrial*nSub];
vector[nStim]       valenced;

for(iSub in 1:nSub){

rho[iSub]          <- exp(x1[iSub]);
gobias[iSub]       <- x3[iSub];
epsilon[iSub]      <- inv_logit(x2[iSub]);
if(epsilon[iSub] < .5){
biaseps[iSub,2]    <- inv_logit(x2[iSub]-x4[iSub]);
biaseps[iSub,1]    <- 2*epsilon[iSub] - biaseps[iSub,2];
} else {
biaseps[iSub,1]    <- inv_logit(x2[iSub]+x4[iSub]);
biaseps[iSub,2]    <- 2*epsilon[iSub] - biaseps[iSub,1];
}
}

for(iSub in 1:nSub){

Q            <- Qi * rho[iSub];
for (iStim in 1:nStim) {
valenced[iStim] <- 0;
}

for(iTrial in 1:nTrial){

for(iResp in 1:nResp){
q[iResp] <- valenced[s[iSub,iTrial]] * Q[iResp,s[iSub,iTrial]]; 
}
for(iResp in 1:2){
q[iResp] <- q[iResp] + gobias[iSub]; 
}

p0         <- softmax(q); 
for(iResp in 1:nResp){
BGx[(iSub-1)*nTrial+iTrial,iResp] <- p0[iResp];
} 


er         <- rho[iSub] * r[iSub,iTrial]; 

if ((ya[iSub,iTrial]==3) && (r[iSub,iTrial]==-1)) 
Q[ya[iSub,iTrial],s[iSub,iTrial]] <- Q[ya[iSub,iTrial],s[iSub,iTrial]] + biaseps[iSub,2] * (er - Q[ya[iSub,iTrial],s[iSub,iTrial]]);
else if ((!ya[iSub,iTrial] == 3) && (r[iSub,iTrial]==1)) 
Q[ya[iSub,iTrial],s[iSub,iTrial]] <- Q[ya[iSub,iTrial],s[iSub,iTrial]] + biaseps[iSub,1] * (er - Q[ya[iSub,iTrial],s[iSub,iTrial]]);
else 
Q[ya[iSub,iTrial],s[iSub,iTrial]] <- Q[ya[iSub,iTrial],s[iSub,iTrial]] + epsilon[iSub] * (er - Q[ya[iSub,iTrial],s[iSub,iTrial]]);


if (r[iSub,iTrial] != 0) 
valenced[s[iSub,iTrial]] <- 1;

}
} 

}



model{
vector[nResp] theta;


X[1]  ~ normal(2, 3); 
X[2]  ~ normal(0, 2); 
X[3]  ~ normal(0, 3); 
X[4]  ~ normal(0, 2); 
sdX   ~ cauchy(0, 2); 


x1 ~ normal(X[1], sdX[1]);
x2 ~ normal(X[2], sdX[2]);
x3 ~ normal(X[3], sdX[3]);
x4 ~ normal(X[4], sdX[4]);


for(iData in 1:nData){
yax[iData] ~ categorical(BGx[iData]);
}
} 



generated quantities {

vector[nData] log_lik; 

for(iData in 1:nData){
log_lik[iData] <- categorical_log(yax[iData], BGx[iData]);
} 

} 

" 




# ========================================================================================.
# mstr5 corresponds to M3c in the paper. M3c includes a learning rate (epsilon), 
# feedback sensitivity (rho), go bias, pavlovian bias (pi), and learning bias (kappa).
# ========================================================================================.
mstr5 = "


data{
int<lower=1>  nTrial;
int<lower=1>  nSub;
int<lower=1>  nStim;
int<lower=1>  nResp;
int<lower=1>  K;
int<lower=1>  s[nSub, nTrial];
int<lower=1>  a[nSub, nTrial];
int<lower=1>  ya[nSub, nTrial];
int<lower=-1> r[nSub, nTrial];
int<lower=1>  rew[nSub, nTrial];
int<lower=1>  yax[nTrial*nSub];
int<lower=1>  nData;

matrix[nResp,nStim]   Qi;
vector[nStim]         Vi;

} 



parameters{


vector[K] X;
vector<lower=0, upper=20>[K] sdX;


vector<lower=-10,  upper=10>[nSub] x1; 
vector<lower=-8,  upper=8>[nSub] x2; 
vector<lower=-8,  upper=8>[nSub] x3; 
vector<lower=-8,  upper=8>[nSub] x4; 
vector<lower=-8,  upper=8>[nSub] x5;  

} 



transformed parameters{

real rho[nSub];
real<lower=0, upper=1> epsilon[nSub];
real gobias[nSub];
real pibias[nSub];
vector[2] biaseps[nSub];

matrix[nResp,nStim] Q;
vector[nResp]       q;
vector[nStim]       V;
simplex[nResp]      p0;
real                er;
simplex[nResp]      BGx[nTrial*nSub];
vector[nStim]       valenced;

for(iSub in 1:nSub){

rho[iSub]          = exp(x1[iSub]);
gobias[iSub]       = x3[iSub];
pibias[iSub]       = x4[iSub];
epsilon[iSub]      = inv_logit(x2[iSub]);
if(epsilon[iSub] < .5){
biaseps[iSub][2]   = inv_logit(x2[iSub]-x5[iSub]);
biaseps[iSub][1]   = 2*epsilon[iSub] - biaseps[iSub][2];
} else {
biaseps[iSub][1]   = inv_logit(x2[iSub]+x5[iSub]);
biaseps[iSub][2]   = 2*epsilon[iSub] - biaseps[iSub][1];
}
}

for(iSub in 1:nSub){

Q            = Qi * rho[iSub];
V            = Vi;
for (iStim in 1:nStim) {
valenced[iStim] = 0;
}

for(iTrial in 1:nTrial){

for(iResp in 1:nResp){
q[iResp] = valenced[s[iSub,iTrial]] * Q[iResp,s[iSub,iTrial]]; 
}
for(iResp in 1:2){
q[iResp] = q[iResp] + gobias[iSub] + pibias[iSub] * valenced[s[iSub,iTrial]] * V[s[iSub,iTrial]]; 
}

p0         = softmax(q);
BGx[(iSub-1)*nTrial+iTrial] = p0; 


er         = rho[iSub] * r[iSub,iTrial]; 

if ((ya[iSub,iTrial]==3) && (r[iSub,iTrial]==-1)) 
Q[ya[iSub,iTrial],s[iSub,iTrial]] = Q[ya[iSub,iTrial],s[iSub,iTrial]] + biaseps[iSub][2] * (er - Q[ya[iSub,iTrial],s[iSub,iTrial]]);
else if ( ((ya[iSub,iTrial] == 1) || (ya[iSub,iTrial] == 2)) && (r[iSub,iTrial] == 1) )   
Q[ya[iSub,iTrial],s[iSub,iTrial]] = Q[ya[iSub,iTrial],s[iSub,iTrial]] + biaseps[iSub][1] * (er - Q[ya[iSub,iTrial],s[iSub,iTrial]]);
    else
Q[ya[iSub,iTrial],s[iSub,iTrial]] = Q[ya[iSub,iTrial],s[iSub,iTrial]] + epsilon[iSub] * (er - Q[ya[iSub,iTrial],s[iSub,iTrial]]);


if (r[iSub,iTrial] != 0) 
valenced[s[iSub,iTrial]] = 1;

} 
} 

} 


model{
vector[nResp] theta;

X[1]  ~ normal(0, 3); 
X[2]  ~ normal(0, 2); 
X[3]  ~ normal(0, 3); 
X[4]  ~ normal(0, 3); 
X[5]  ~ normal(0, 2); 
sdX   ~ cauchy(0, 2); 

x1 ~ normal(X[1], sdX[1]);
x2 ~ normal(X[2], sdX[2]);
x3 ~ normal(X[3], sdX[3]);
x4 ~ normal(X[4], sdX[4]);
x5 ~ normal(X[5], sdX[5]);

for(iData in 1:nData){
yax[iData] ~ categorical(BGx[iData]);
}
} 


generated quantities {

vector[nData] log_lik; 

for(iData in 1:nData){
log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
} 

} 

"




# ========================================================================================
# mstr6 and 7 are control models for M1 in the paper. M1 includes a learning rate
# (epsilon) and feedback sensitivity (rho).
# These control models don't "re-evaluate" the Q-value after 1st reward/punishment.
# mstr6 assumes that Qi = +/- .5 * rho after 1st reward/punishment, while mstr7 assumes Qi = 0.
# ========================================================================================
mstr6 = "
data {
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];        // state (cue)
  int<lower=1>  a[nSub, nTrial];        // required response (unused here)
  int<lower=1>  ya[nSub, nTrial];       // chosen response (1..nResp)
  int<lower=-1> r[nSub, nTrial];        // outcome: -1/0/1
  int<lower=1>  rew[nSub, nTrial];      // reward flag (unused here)
  int<lower=1>  yax[nTrial * nSub];     // flattened choices for likelihood
  int<lower=1>  nData;

  matrix[nResp, nStim] Qi;              // initial Q table
  vector[nStim]        Vi;              // not used here (kept for compatibility)
}

// -------------------------------------------------------------------------
// parameters in 'fitting space'
// -------------------------------------------------------------------------
parameters {
  // Hierarchical (group-level) parameters
  vector[K] X;                          // group means (e.g., log rho, logit epsilon)
  vector<lower=0, upper=20>[K] sdX;     // group SDs (half-Cauchy prior)

  // Subject-level parameters
  vector<lower=-10, upper=10>[nSub] x1; // individual log(rho)
  vector<lower=-8,  upper=8>[nSub] x2;  // individual logit(epsilon)
}

// -------------------------------------------------------------------------
// transformed parameters: subject parameters + trial-by-trial choice probs
// -------------------------------------------------------------------------
transformed parameters {
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp] q;
  simplex[nResp] p0;
  real er;
  vector[nResp] BGx[nTrial * nSub];     // choice probabilities per trial
  vector[nStim] valenced;               // 0/1 flags: whether cue valence is known

  // retrieve sampled subject parameters
  for (iSub in 1:nSub) {
    rho[iSub]     = exp(x1[iSub]);       // feedback sensitivity > 0
    epsilon[iSub] = inv_logit(x2[iSub]); // learning rate in (0,1)
  }

  // compute subject- and trial-specific choice probabilities
  for (iSub in 1:nSub) {
    // initialize Q and valence-known flags for this subject
    Q = Qi * rho[iSub];

    for (iStim in 1:nStim) {
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial) {
      int idx;
      idx = (iSub - 1) * nTrial + iTrial;

      // Q-value only contributes once valence is known; otherwise 0
      for (iResp in 1:nResp) {
        q[iResp] = valenced[s[iSub, iTrial]] * Q[iResp, s[iSub, iTrial]];
      }

      // softmax over action values
      p0 = softmax(q);

      // store choice probabilities
      for (iResp in 1:nResp) {
        BGx[idx, iResp] = p0[iResp];
      }

      // update whether valence is known (first non-zero outcome)
      if (r[iSub, iTrial] != 0)
        valenced[s[iSub, iTrial]] = 1;

      // once valence is known, update Q-values via delta rule
      if (valenced[s[iSub, iTrial]] == 1) {
        er = rho[iSub] * r[iSub, iTrial];  // scaled outcome

        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          epsilon[iSub] *
          (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      }
    } // end iTrial
  } // end iSub
}

// -------------------------------------------------------------------------
// model: priors + likelihood
// -------------------------------------------------------------------------
model {
  // Hierarchical priors
  X[1] ~ normal(2, 3);   // prior for log rho
  X[2] ~ normal(0, 2);   // prior for logit epsilon
  sdX  ~ cauchy(0, 2);   // half-Cauchy on SDs

  // Subject-level priors
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);

  // Likelihood of chosen actions
  for (iData in 1:nData) {
    yax[iData] ~ categorical(BGx[iData]);
  }
}

// -------------------------------------------------------------------------
// generated quantities: pointwise log-likelihood (for WAIC/LOO)
// -------------------------------------------------------------------------
generated quantities {
  vector[nData] log_lik;

  for (iData in 1:nData) {
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}
"






# ========================================================================================.
# mstr6 and 7 are control models for M1 in the paper. M1 includes a learning rate (epsilon), 
# and feedback sensitivity (rho).
# these control model don't "re-evaluate" the Q-value after 1st reward/punishment.
# mstr6 assumes that Qi=+/-.5*rho after 1st rew/pun, while mstr7 assumes Qi=0.
# ========================================================================================.
mstr7 = mstr6





# ========================================================================================
# mstr8 is a control model for M3a in the paper. M3a includes a learning rate (epsilon), 
# feedback sensitivity (rho), Go bias, and Pavlovian bias (pi).
# This control model implements learning of the Pavlovian value V(s).
# ========================================================================================
mstr8 = "
data {
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];        // state (cue)
  int<lower=1>  a[nSub, nTrial];        // required response (unused here)
  int<lower=1>  ya[nSub, nTrial];       // chosen response (1..nResp)
  int<lower=-1> r[nSub, nTrial];        // outcome: -1/0/1
  int<lower=1>  rew[nSub, nTrial];      // reward flag (unused here)
  int<lower=1>  yax[nTrial * nSub];     // flattened choices for likelihood
  int<lower=1>  nData;

  matrix[nResp, nStim] Qi;              // initial Q table
  vector[nStim]        Vi;              // initial Pavlovian values (used & learned)
}

parameters {
  // group-level parameters
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // subject-level parameters
  vector<lower=-10, upper=10>[nSub] x1; // log(rho)
  vector<lower=-8,  upper=8>[nSub] x2;  // logit(epsilon)
  vector<lower=-8,  upper=8>[nSub] x3;  // Go bias
  vector<lower=-8,  upper=8>[nSub] x4;  // Pavlovian bias
}

transformed parameters {
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp] q;
  vector[nStim] V;
  vector[nResp] p0;
  real er;
  vector[nResp] BGx[nTrial * nSub];
  vector[nStim] valenced;

  // subject-level parameters in natural space
  for (iSub in 1:nSub) {
    rho[iSub]     = exp(x1[iSub]);
    epsilon[iSub] = inv_logit(x2[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
  }

  // trial-by-trial choice probabilities
  for (iSub in 1:nSub) {
    // initialize Q and Pavlovian values for this subject
    Q = Qi * rho[iSub];
    V = Vi;

    // valence-known flags start at 0
    for (iStim in 1:nStim) {
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial) {
      int idx;
      idx = (iSub - 1) * nTrial + iTrial;

      // base action values: only contribute once valence is known
      for (iResp in 1:nResp) {
        q[iResp] = valenced[s[iSub, iTrial]] * Q[iResp, s[iSub, iTrial]];
      }

      // add Go bias + Pavlovian bias * V(s) to the two Go responses
      for (iResp in 1:2) {
        q[iResp] = q[iResp]
                   + gobias[iSub]
                   + pibias[iSub] * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
      }

      // softmax to get choice probabilities
      p0 = softmax(q);

      // store for likelihood
      for (iResp in 1:nResp) {
        BGx[idx, iResp] = p0[iResp];
      }

      // scaled outcome
      er = rho[iSub] * r[iSub, iTrial];

      // instrumental Q-learning update (delta rule)
      Q[ya[iSub, iTrial], s[iSub, iTrial]] =
        Q[ya[iSub, iTrial], s[iSub, iTrial]]
        + epsilon[iSub] *
          (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);

      // Pavlovian V-learning update (Rescorla–Wagner style)
      V[s[iSub, iTrial]] =
        V[s[iSub, iTrial]]
        + epsilon[iSub] *
          (r[iSub, iTrial] - V[s[iSub, iTrial]]);

      // mark this stimulus as having known valence once non-zero outcome occurs
      if (r[iSub, iTrial] != 0)
        valenced[s[iSub, iTrial]] = 1;
    }
  }
}

model {
  // group-level priors
  X[1] ~ normal(0, 3); // log rho
  X[2] ~ normal(0, 2); // logit epsilon
  X[3] ~ normal(0, 3); // Go bias
  X[4] ~ normal(0, 3); // Pavlovian bias
  sdX  ~ cauchy(0, 2);

  // subject-level priors
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);

  // likelihood
  for (iData in 1:nData) {
    yax[iData] ~ categorical(BGx[iData]);
  }
}

generated quantities {
  vector[nData] log_lik;

  for (iData in 1:nData) {
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}
"


mstr9 = "

// ========================================================================================.
// mstr9 has exactly the same computational logic as mstr5 (no EEG used).
// Subject exclusion, if needed, should be handled during data preparation.
// This version only modernizes syntax without changing any computation or update rule.
// ========================================================================================.

data{
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];
  int<lower=1>  a[nSub, nTrial];
  int<lower=1>  ya[nSub, nTrial];
  int<lower=-1> r[nSub, nTrial];
  int<lower=1>  rew[nSub, nTrial];
  int<lower=1>  yax[nTrial * nSub];
  int<lower=1>  nData;

  // Retained only to match the data structure (not used in this model)
  int<lower=1> motivconflict[nStim];
  real<lower=-1,upper=1> tPow[nSub, nTrial];
  real<lower=-1,upper=1> pfcICPC[nSub, nTrial];
  real<lower=-1,upper=1> lICPC[nSub, nTrial];
  real<lower=-1,upper=1> rICPC[nSub, nTrial];

  matrix[nResp, nStim] Qi;
  vector[nStim]        Vi;
}

parameters{
  // Group level
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // Subject level
  vector<lower=-10, upper=10>[nSub] x1; // rho (log-space)
  vector<lower= -8, upper=  8>[nSub] x2; // epsilon (logit-space)
  vector<lower= -8, upper=  8>[nSub] x3; // go bias
  vector<lower= -8, upper=  8>[nSub] x4; // pavlovian bias
  vector<lower= -8, upper=  8>[nSub] x5; // learning bias
}

transformed parameters{
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  vector[2] biaseps[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp] q;
  vector[nStim] V;
  simplex[nResp] p0;                  // softmax probabilities
  real er;
  simplex[nResp] BGx[nTrial * nSub];  // choice probabilities for each (sub × trial)
  vector[nStim] valenced;

  // Parameter transforms
  for (iSub in 1:nSub){
    rho[iSub]     = exp(x1[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    epsilon[iSub] = inv_logit(x2[iSub]);

    if (epsilon[iSub] < 0.5){
      biaseps[iSub, 2] = inv_logit(x2[iSub] - x5[iSub]);
      biaseps[iSub, 1] = 2 * epsilon[iSub] - biaseps[iSub, 2];
    } else {
      biaseps[iSub, 1] = inv_logit(x2[iSub] + x5[iSub]);
      biaseps[iSub, 2] = 2 * epsilon[iSub] - biaseps[iSub, 1];
    }
  }

  // Q-learning and choice probability generation
  for (iSub in 1:nSub){
    Q = Qi * rho[iSub];
    V = Vi;
    for (iStim in 1:nStim) valenced[iStim] = 0;

    for (iTrial in 1:nTrial){
      // Action values for the current stimulus
      for (iResp in 1:nResp)
        q[iResp] = valenced[s[iSub, iTrial]] * Q[iResp, s[iSub, iTrial]];

      // Add go/pavlovian biases only to the two go responses
      for (iResp in 1:2)
        q[iResp] = q[iResp]
                   + gobias[iSub]
                   + pibias[iSub] * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];

      // Softmax probabilities
      p0 = softmax(q);
      BGx[(iSub - 1) * nTrial + iTrial] = p0;

      // RPE update
      er = rho[iSub] * r[iSub, iTrial];

      if ( (ya[iSub, iTrial] == 3) && (r[iSub, iTrial] == -1) ) {
        // punished NoGo
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] =
          Q[ ya[iSub, iTrial], s[iSub, iTrial] ]
          + biaseps[iSub, 2] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      } else if ( (ya[iSub, iTrial] != 3) && (r[iSub, iTrial] == 1) ) {
        // rewarded Go
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] =
          Q[ ya[iSub, iTrial], s[iSub, iTrial] ]
          + biaseps[iSub, 1] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      } else {
        // neutral outcome
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] =
          Q[ ya[iSub, iTrial], s[iSub, iTrial] ]
          + epsilon[iSub] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      }

      // Once a stimulus receives non-zero feedback, its valence is treated as known
      if (r[iSub, iTrial] != 0)
        valenced[ s[iSub, iTrial] ] = 1;
    }
  }
}

model{
  // Priors (same as original)
  X[1]  ~ normal(0, 3); // feedback sensitivity
  X[2]  ~ normal(0, 2); // instrumental learning rate
  X[3]  ~ normal(0, 3); // go bias
  X[4]  ~ normal(0, 3); // pavlovian bias
  X[5]  ~ normal(0, 2); // learning bias
  sdX   ~ cauchy(0, 2); // half-Cauchy (lower bound = 0)

  // Hierarchical structure
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);

  // Likelihood
  for (iData in 1:nData)
    yax[iData] ~ categorical( BGx[iData] );
}

generated quantities{
  vector[nData] log_lik;
  for (iData in 1:nData)
    log_lik[iData] = categorical_lpmf( yax[iData] | BGx[iData] );
}

"




mstr10 = "

// ========================================================================================.
// M4a: epsilon, rho, go bias, pavlovian bias (pi), learning bias (kappa),
//      plus modulation of pi by midfrontal theta (only on motivational-conflict trials).
// This version only modernizes syntax; the computational logic is unchanged.
// ========================================================================================.

data{
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];
  int<lower=1>  a[nSub, nTrial];
  int<lower=1>  ya[nSub, nTrial];
  int<lower=-1> r[nSub, nTrial];
  int<lower=1>  rew[nSub, nTrial];
  int<lower=1>  yax[nTrial * nSub];
  int<lower=1>  nData;
  int<lower=1>  motivconflict[nStim];

  // EEG inputs: only tPow is used; others are retained to align with the data structure
  real<lower=-1,upper=1> tPow[nSub, nTrial];
  real<lower=-1,upper=1> pfcICPC[nSub, nTrial];
  real<lower=-1,upper=1> lICPC[nSub, nTrial];
  real<lower=-1,upper=1> rICPC[nSub, nTrial];

  matrix[nResp, nStim] Qi;
  vector[nStim]        Vi;
}

parameters{
  // Group level
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // Subject level
  vector<lower=-10, upper=10>[nSub] x1; // rho (log-space)
  vector<lower= -8, upper=  8>[nSub] x2; // epsilon (logit-space)
  vector<lower= -8, upper=  8>[nSub] x3; // go bias
  vector<lower= -8, upper=  8>[nSub] x4; // pavlovian bias
  vector<lower= -8, upper=  8>[nSub] x5; // learning bias
  vector<lower= -8, upper=  8>[nSub] x6; // betaEEG
}

transformed parameters{
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  vector[2] biaseps[nSub];
  real betaEEG[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp] q;
  vector[nStim] V;
  simplex[nResp] p0;                   // softmax probabilities
  real er;
  simplex[nResp] BGx[nTrial * nSub];   // choice probabilities for each (sub × trial)
  vector[nStim] valenced;

  // Parameter transforms
  for (iSub in 1:nSub){
    rho[iSub]     = exp(x1[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    epsilon[iSub] = inv_logit(x2[iSub]);

    if (epsilon[iSub] < 0.5){
      biaseps[iSub, 2] = inv_logit(x2[iSub] - x5[iSub]);
      biaseps[iSub, 1] = 2 * epsilon[iSub] - biaseps[iSub, 2];
    } else {
      biaseps[iSub, 1] = inv_logit(x2[iSub] + x5[iSub]);
      biaseps[iSub, 2] = 2 * epsilon[iSub] - biaseps[iSub, 1];
    }
    betaEEG[iSub] = x6[iSub];
  }

  // Q-learning and choice probability generation
  for (iSub in 1:nSub){
    Q = Qi * rho[iSub];
    V = Vi;
    for (iStim in 1:nStim) valenced[iStim] = 0;

    for (iTrial in 1:nTrial){
      // Action values for the current stimulus
      for (iResp in 1:nResp)
        q[iResp] = valenced[s[iSub, iTrial]] * Q[iResp, s[iSub, iTrial]];

      // Add go/pavlovian biases only to the two go responses;
      // for motivational-conflict trials (=1), use (pi + betaEEG * tPow),
      // otherwise use pi only (same as original)
      for (iResp in 1:2){
        if (motivconflict[s[iSub, iTrial]] == 1)
          q[iResp] = q[iResp]
                     + gobias[iSub]
                     + (pibias[iSub] + betaEEG[iSub] * tPow[iSub, iTrial])
                       * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
        else
          q[iResp] = q[iResp]
                     + gobias[iSub]
                     + pibias[iSub] * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
      }

      // Softmax probabilities
      p0 = softmax(q);
      BGx[(iSub - 1) * nTrial + iTrial] = p0;

      // RPE update
      er = rho[iSub] * r[iSub, iTrial];

      if ( (ya[iSub, iTrial] == 3) && (r[iSub, iTrial] == -1) ) {
        // punished NoGo
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] =
          Q[ ya[iSub, iTrial], s[iSub, iTrial] ]
          + biaseps[iSub, 2] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      } else if ( (ya[iSub, iTrial] != 3) && (r[iSub, iTrial] == 1) ) {
        // rewarded Go
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] =
          Q[ ya[iSub, iTrial], s[iSub, iTrial] ]
          + biaseps[iSub, 1] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      } else {
        // neutral outcome
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] =
          Q[ ya[iSub, iTrial], s[iSub, iTrial] ]
          + epsilon[iSub] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      }

      // Once a stimulus receives non-zero feedback, its valence is treated as known
      if (r[iSub, iTrial] != 0)
        valenced[ s[iSub, iTrial] ] = 1;
    }
  }
}

model{
  // Group-level priors (same as original)
  X[1]  ~ normal(0, 3); // feedback sensitivity
  X[2]  ~ normal(0, 2); // instrumental learning rate
  X[3]  ~ normal(0, 3); // go bias
  X[4]  ~ normal(0, 3); // pavlovian bias
  X[5]  ~ normal(0, 2); // learning bias
  X[6]  ~ normal(0, 3); // EEG weight
  sdX   ~ cauchy(0, 2); // half-Cauchy (lower bound = 0)

  // Hierarchical structure
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);

  // Likelihood
  for (iData in 1:nData)
    yax[iData] ~ categorical( BGx[iData] );
}

generated quantities{
  vector[nData] log_lik; // for WAIC
  for (iData in 1:nData)
    log_lik[iData] = categorical_lpmf( yax[iData] | BGx[iData] );
}

"


# ========================================================================================.
# mstr11 is a control model for M5a in the paper. In this control model midfrontal theta 
# has an effect on the Pavlovian bias on all trials.
# ========================================================================================.
mstr11 = "

data {
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];        // state (cue)
  int<lower=1>  a[nSub, nTrial];        // required response (unused here)
  int<lower=1>  ya[nSub, nTrial];       // chosen response (1..nResp)
  int<lower=-1> r[nSub, nTrial];        // outcome: -1/0/1
  int<lower=1>  rew[nSub, nTrial];      // reward flag (unused here)
  int<lower=1>  yax[nTrial * nSub];     // flattened choices for likelihood
  int<lower=1>  nData;
  int<lower=1>  motivconflict[nStim];   // not used in this control model

  real<lower=-1, upper=1> tPow[nSub, nTrial];   // midfrontal theta
  real<lower=-1, upper=1> pfcICPC[nSub, nTrial]; // unused here
  real<lower=-1, upper=1> lICPC[nSub, nTrial];   // unused here
  real<lower=-1, upper=1> rICPC[nSub, nTrial];   // unused here

  matrix[nResp, nStim] Qi;              // initial Q table
  vector[nStim]        Vi;              // initial Pavlovian values
}

// -------------------------------------------------------------------------
// parameters in 'fitting space'
// -------------------------------------------------------------------------
parameters {
  // Hierarchical (group-level) parameters
  vector[K] X;                          // group means: rho, epsilon, go, pi, eps-bias, betaEEG
  vector<lower=0, upper=20>[K] sdX;     // group sds (half-Cauchy prior)

  // Subject-level parameters
  vector<lower=-10, upper=10>[nSub] x1; // log(rho)
  vector<lower=-8,  upper=8>[nSub] x2;  // logit(epsilon)
  vector<lower=-8,  upper=8>[nSub] x3;  // go bias
  vector<lower=-8,  upper=8>[nSub] x4;  // pavlovian bias
  vector<lower=-8,  upper=8>[nSub] x5;  // learning bias on epsilon
  vector<lower=-8,  upper=8>[nSub] x6;  // betaEEG (theta weight)
}

// -------------------------------------------------------------------------
// transformed parameters: subject parameters + trial-by-trial choice probs
// -------------------------------------------------------------------------
transformed parameters {
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  vector[2] biaseps[nSub];              // learning rates for rewarded Go / punished NoGo
  real betaEEG[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp] q;
  vector[nStim] V;
  simplex[nResp] p0;
  real er;
  vector[nResp] BGx[nTrial * nSub];     // choice probabilities per trial
  vector[nStim] valenced;              // whether cue valence is known

  // subject-level transforms
  for (iSub in 1:nSub) {
    rho[iSub]     = exp(x1[iSub]);       // feedback sensitivity > 0
    epsilon[iSub] = inv_logit(x2[iSub]); // base learning rate in (0,1)
    gobias[iSub]  = x3[iSub];            // Go bias
    pibias[iSub]  = x4[iSub];            // Pavlovian bias
    betaEEG[iSub] = x6[iSub];            // EEG weight

    // construct valence-specific learning rates biaseps[1] (rewarded Go) and biaseps[2] (punished NoGo)
    if (epsilon[iSub] < 0.5) {
      biaseps[iSub][2] = inv_logit(x2[iSub] - x5[iSub]);           // punished NoGo
      biaseps[iSub][1] = 2 * epsilon[iSub] - biaseps[iSub][2];     // rewarded Go
    } else {
      biaseps[iSub][1] = inv_logit(x2[iSub] + x5[iSub]);           // rewarded Go
      biaseps[iSub][2] = 2 * epsilon[iSub] - biaseps[iSub][1];     // punished NoGo
    }
  }

  // compute trial-by-trial choice probabilities
  for (iSub in 1:nSub) {
    // initialize Q and Pavlovian values for this subject
    Q = Qi * rho[iSub];
    V = Vi;

    // initially, valence of each cue is unknown
    for (iStim in 1:nStim) {
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial) {
      int idx;
      idx = (iSub - 1) * nTrial + iTrial;

      // 1) compute action values q
      // instrumental Q part: only when valence is known
      for (iResp in 1:nResp) {
        q[iResp] = valenced[s[iSub, iTrial]] * Q[iResp, s[iSub, iTrial]];
      }

      // add Go bias and Pavlovian bias (with EEG modulation) to Go responses (1 and 2)
      for (iResp in 1:2) {
        q[iResp] = q[iResp]
                   + gobias[iSub]
                   + (pibias[iSub] + betaEEG[iSub] * tPow[iSub, iTrial])
                     * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
      }

      // 2) softmax to get choice probabilities
      p0 = softmax(q);

      for (iResp in 1:nResp) {
        BGx[idx, iResp] = p0[iResp];
      }

      // 3) update Q-values with valence-specific learning bias
      er = rho[iSub] * r[iSub, iTrial];  // scaled outcome

      if ((ya[iSub, iTrial] == 3) && (r[iSub, iTrial] == -1)) {
        // punished NoGo
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          biaseps[iSub][2] * (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      } else if ((ya[iSub, iTrial] != 3) && (r[iSub, iTrial] == 1)) {
        // rewarded Go
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          biaseps[iSub][1] * (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      } else {
        // neutral or other outcomes use base epsilon
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          epsilon[iSub] * (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      }

      // 4) mark that cue valence is now known if outcome is non-zero
      if (r[iSub, iTrial] != 0) {
        valenced[s[iSub, iTrial]] = 1;
      }
    } // end iTrial
  }   // end iSub
}

// -------------------------------------------------------------------------
// model: priors + likelihood
// -------------------------------------------------------------------------
model {
  // Hierarchical priors
  X[1] ~ normal(0, 3);   // feedback sensitivity (log rho)
  X[2] ~ normal(0, 2);   // instrumental learning rate (logit epsilon)
  X[3] ~ normal(0, 3);   // Go bias
  X[4] ~ normal(0, 3);   // Pavlovian bias
  X[5] ~ normal(0, 2);   // learning bias
  X[6] ~ normal(0, 3);   // EEG weight
  sdX  ~ cauchy(0, 2);   // half-Cauchy on SDs

  // Subject-level priors
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);

  // Likelihood of chosen actions
  for (iData in 1:nData) {
    yax[iData] ~ categorical(BGx[iData]);
  }
}

// -------------------------------------------------------------------------
// generated quantities: pointwise log-likelihood (for WAIC/LOO)
// -------------------------------------------------------------------------
generated quantities {
  vector[nData] log_lik;

  for (iData in 1:nData) {
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}

"  # end mstr11





mstr12 = "

// ========================================================================================.
// M4b (control model: theta affects the instrumental Q contribution only on conflict trials).
// Syntax modernized; computational logic unchanged.
// ========================================================================================.

data{
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];
  int<lower=1>  a[nSub, nTrial];
  int<lower=1>  ya[nSub, nTrial];
  int<lower=-1> r[nSub, nTrial];
  int<lower=1>  rew[nSub, nTrial];
  int<lower=1>  yax[nTrial * nSub];
  int<lower=1>  nData;

  int<lower=1>  motivconflict[nStim];   // 1 = conflict, 2 = non-conflict

  // EEG inputs (only tPow is used in this model; others are retained for alignment)
  real<lower=-1,upper=1> tPow[nSub, nTrial];
  real<lower=-1,upper=1> pfcICPC[nSub, nTrial];
  real<lower=-1,upper=1> lICPC[nSub, nTrial];
  real<lower=-1,upper=1> rICPC[nSub, nTrial];

  matrix[nResp, nStim] Qi;
  vector[nStim]        Vi;
}

parameters{
  // Group level
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // Subject level
  vector<lower=-10, upper=10>[nSub] x1; // rho (log-space)
  vector<lower= -8, upper=  8>[nSub] x2; // epsilon (logit-space)
  vector<lower= -8, upper=  8>[nSub] x3; // go bias
  vector<lower= -8, upper=  8>[nSub] x4; // pavlovian bias
  vector<lower= -8, upper=  8>[nSub] x5; // learning bias
  vector<lower= -8, upper=  8>[nSub] x6; // betaEEG
}

transformed parameters{
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  vector[2] biaseps[nSub];
  real betaEEG[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp] q;
  vector[nStim] V;
  simplex[nResp] p0;                   // softmax probabilities
  real er;
  simplex[nResp] BGx[nTrial * nSub];   // choice probabilities for each (sub × trial)
  vector[nStim] valenced;

  // Parameter transforms
  for (iSub in 1:nSub){
    rho[iSub]     = exp(x1[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    epsilon[iSub] = inv_logit(x2[iSub]);

    if (epsilon[iSub] < 0.5){
      biaseps[iSub, 2] = inv_logit(x2[iSub] - x5[iSub]);
      biaseps[iSub, 1] = 2 * epsilon[iSub] - biaseps[iSub, 2];
    } else {
      biaseps[iSub, 1] = inv_logit(x2[iSub] + x5[iSub]);
      biaseps[iSub, 2] = 2 * epsilon[iSub] - biaseps[iSub, 1];
    }
    betaEEG[iSub] = x6[iSub];
  }

  // Q-learning and choice probability generation
  for (iSub in 1:nSub){
    Q = Qi * rho[iSub];
    V = Vi;
    for (iStim in 1:nStim) valenced[iStim] = 0;

    for (iTrial in 1:nTrial){
      // Action values for the current stimulus
      for (iResp in 1:nResp){
        q[iResp] = valenced[s[iSub, iTrial]] * Q[iResp, s[iSub, iTrial]];
        // Scale the Q contribution only on conflict trials (affects all responses)
        if (motivconflict[s[iSub, iTrial]] == 1)
          q[iResp] = (1 - betaEEG[iSub] * tPow[iSub, iTrial]) * q[iResp];
      }

      // Add go/pavlovian contributions to the Go channels
      // (same as original: not scaled by theta)
      for (iResp in 1:2){
        q[iResp] = q[iResp]
                   + gobias[iSub]
                   + pibias[iSub] * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
      }

      // Softmax probabilities
      p0 = softmax(q);
      BGx[(iSub - 1) * nTrial + iTrial] = p0;

      // RPE update
      er = rho[iSub] * r[iSub, iTrial];

      if ( (ya[iSub, iTrial] == 3) && (r[iSub, iTrial] == -1) ) {
        // punished NoGo
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] =
          Q[ ya[iSub, iTrial], s[iSub, iTrial] ]
          + biaseps[iSub, 2] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      } else if ( (ya[iSub, iTrial] != 3) && (r[iSub, iTrial] == 1) ) {
        // rewarded Go
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] =
          Q[ ya[iSub, iTrial], s[iSub, iTrial] ]
          + biaseps[iSub, 1] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      } else {
        // neutral outcome
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] =
          Q[ ya[iSub, iTrial], s[iSub, iTrial] ]
          + epsilon[iSub] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      }

      // Once a stimulus receives non-zero feedback, its valence is treated as known
      if (r[iSub, iTrial] != 0)
        valenced[ s[iSub, iTrial] ] = 1;
    }
  }
}

model{
  // Group-level priors (same as original)
  X[1]  ~ normal(0, 3); // feedback sensitivity
  X[2]  ~ normal(0, 2); // instrumental learning rate
  X[3]  ~ normal(0, 3); // go bias
  X[4]  ~ normal(0, 3); // pavlovian bias
  X[5]  ~ normal(0, 2); // learning bias
  X[6]  ~ normal(0, 3); // EEG weight
  sdX   ~ cauchy(0, 2); // half-Cauchy (lower bound = 0)

  // Hierarchical structure
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);

  // Likelihood
  for (iData in 1:nData)
    yax[iData] ~ categorical( BGx[iData] );
}

generated quantities{
  vector[nData] log_lik; // for WAIC
  for (iData in 1:nData)
    log_lik[iData] = categorical_lpmf( yax[iData] | BGx[iData] );
}

"













mstr13 = "

// ========================================================================================.
// M4c: theta affects the trade-off between instrumental (Q) and Pavlovian (V)
//      components only on conflict trials; syntax modernized, logic unchanged.
// ========================================================================================.

data{
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];
  int<lower=1>  a[nSub, nTrial];
  int<lower=1>  ya[nSub, nTrial];
  int<lower=-1> r[nSub, nTrial];
  int<lower=1>  rew[nSub, nTrial];
  int<lower=1>  yax[nTrial * nSub];
  int<lower=1>  nData;

  int<lower=1>  motivconflict[nStim];   // 1 = conflict, 2 = non-conflict

  // EEG inputs (retained for structural alignment; only tPow is used here)
  real<lower=-1,upper=1> tPow[nSub, nTrial];
  real<lower=-1,upper=1> pfcICPC[nSub, nTrial];
  real<lower=-1,upper=1> lICPC[nSub, nTrial];
  real<lower=-1,upper=1> rICPC[nSub, nTrial];

  matrix[nResp, nStim] Qi;
  vector[nStim]        Vi;
}

parameters{
  // Group level
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // Subject level
  vector<lower=-10, upper=10>[nSub] x1; // rho (log-space)
  vector<lower= -8, upper=  8>[nSub] x2; // epsilon (logit-space)
  vector<lower= -8, upper=  8>[nSub] x3; // go bias
  vector<lower= -8, upper=  8>[nSub] x4; // wtradeoff
  vector<lower= -8, upper=  8>[nSub] x5; // learning bias
  vector<lower= -8, upper=  8>[nSub] x6; // betaEEG
}

transformed parameters{
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real wtradeoff[nSub];
  vector[2] biaseps[nSub];
  real betaEEG[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp] q;
  vector[nStim] V;
  simplex[nResp] p0;                    // softmax probabilities
  real er;
  simplex[nResp] BGx[nTrial * nSub];    // choice probabilities for each (sub × trial)
  vector[nStim] valenced;

  // Parameter transforms
  for (iSub in 1:nSub){
    rho[iSub]      = exp(x1[iSub]);
    gobias[iSub]   = x3[iSub];
    wtradeoff[iSub]= x4[iSub];
    epsilon[iSub]  = inv_logit(x2[iSub]);

    if (epsilon[iSub] < 0.5){
      biaseps[iSub, 2] = inv_logit(x2[iSub] - x5[iSub]);
      biaseps[iSub, 1] = 2 * epsilon[iSub] - biaseps[iSub, 2];
    } else {
      biaseps[iSub, 1] = inv_logit(x2[iSub] + x5[iSub]);
      biaseps[iSub, 2] = 2 * epsilon[iSub] - biaseps[iSub, 1];
    }
    betaEEG[iSub] = x6[iSub];
  }

  // Q-learning and choice probability generation
  for (iSub in 1:nSub){
    Q = Qi * rho[iSub];
    V = Vi;
    for (iStim in 1:nStim) valenced[iStim] = 0;

    for (iTrial in 1:nTrial){
      // Instrumental component: start with Q, then apply trade-off scaling
      for (iResp in 1:nResp){
        q[iResp] = valenced[s[iSub, iTrial]] * Q[iResp, s[iSub, iTrial]];
        if (motivconflict[s[iSub, iTrial]] == 1){
          // Conflict: Q × (1 - (w + β·theta))
          q[iResp] = (1 - (wtradeoff[iSub] + betaEEG[iSub] * tPow[iSub, iTrial])) * q[iResp];
        } else {
          // Non-conflict: Q × (1 - w)
          q[iResp] = (1 - wtradeoff[iSub]) * q[iResp];
        }
      }

      // Pavlovian component is added only to the Go channels;
      // on conflict trials the weight is (w + β·theta), otherwise w
      for (iResp in 1:2){
        if (motivconflict[s[iSub, iTrial]] == 1){
          q[iResp] = q[iResp]
                     + gobias[iSub]
                     + (wtradeoff[iSub] + betaEEG[iSub] * tPow[iSub, iTrial])
                       * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
        } else {
          q[iResp] = q[iResp]
                     + gobias[iSub]
                     + wtradeoff[iSub] * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
        }
      }

      // Softmax probabilities
      p0 = softmax(q);
      BGx[(iSub - 1) * nTrial + iTrial] = p0;

      // RPE update
      er = rho[iSub] * r[iSub, iTrial];

      if ( (ya[iSub, iTrial] == 3) && (r[iSub, iTrial] == -1) ) {
        // punished NoGo
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] =
          Q[ ya[iSub, iTrial], s[iSub, iTrial] ]
          + biaseps[iSub, 2] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      } else if ( (ya[iSub, iTrial] != 3) && (r[iSub, iTrial] == 1) ) {
        // rewarded Go
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] =
          Q[ ya[iSub, iTrial], s[iSub, iTrial] ]
          + biaseps[iSub, 1] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      } else {
        // neutral outcome
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] =
          Q[ ya[iSub, iTrial], s[iSub, iTrial] ]
          + epsilon[iSub] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      }

      // Once a stimulus receives non-zero feedback, its valence is treated as known
      if (r[iSub, iTrial] != 0)
        valenced[ s[iSub, iTrial] ] = 1;
    }
  }
}

model{
  // Group-level priors (keep original settings)
  X[1]  ~ normal(0, 3); // feedback sensitivity
  X[2]  ~ normal(0, 2); // instrumental learning rate
  X[3]  ~ normal(0, 3); // go bias
  X[4]  ~ normal(0, 3); // trade-off weight
  X[5]  ~ normal(0, 2); // learning bias
  X[6]  ~ normal(0, 3); // EEG weight
  sdX   ~ cauchy(0, 2); // half-Cauchy

  // Hierarchical structure
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);

  // Likelihood
  for (iData in 1:nData)
    yax[iData] ~ categorical( BGx[iData] );
}

generated quantities{
  vector[nData] log_lik; // for WAIC
  for (iData in 1:nData)
    log_lik[iData] = categorical_lpmf( yax[iData] | BGx[iData] );
}

"




mstr14 = "

// ========================================================================================.
// M4d: theta affects learning bias only on conflict trials;
//      betaEEG ∈ [0,1]; syntax modernized, logic unchanged.
// ========================================================================================.

data{
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];
  int<lower=1>  a[nSub, nTrial];
  int<lower=1>  ya[nSub, nTrial];
  int<lower=-1> r[nSub, nTrial];
  int<lower=1>  rew[nSub, nTrial];
  int<lower=1>  yax[nTrial * nSub];
  int<lower=1>  nData;

  int<lower=1>  motivconflict[nStim];         // 1 = conflict, 2 = non-conflict

  real<lower=-1,upper=1> tPow[nSub, nTrial];   // only tPow is used here; others are retained for alignment
  real<lower=-1,upper=1> pfcICPC[nSub, nTrial];
  real<lower=-1,upper=1> lICPC[nSub, nTrial];
  real<lower=-1,upper=1> rICPC[nSub, nTrial];

  matrix[nResp, nStim] Qi;
  vector[nStim]        Vi;
}

parameters{
  // Group level
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // Subject level
  vector<lower=-10, upper=10>[nSub] x1; // rho (log-space)
  vector<lower= -8, upper=  8>[nSub] x2; // epsilon (logit-space)
  vector<lower= -8, upper=  8>[nSub] x3; // go bias
  vector<lower= -8, upper=  8>[nSub] x4; // pavlovian bias (x4 in the original model)
  vector<lower= -8, upper=  8>[nSub] x5; // learning bias
  vector<lower= -8, upper=  8>[nSub] x6; // betaEEG (mapped to [0,1] by inv_logit)
}

transformed parameters{
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  vector[2] biaseps[nSub];
  real<lower=0, upper=1> betaEEG[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp] q;
  vector[nStim] V;
  simplex[nResp] p0;                       // softmax probabilities
  real er;
  real lr;
  simplex[nResp] BGx[nTrial * nSub];       // choice probabilities for each (sub × trial)
  vector[nStim] valenced;

  // Parameter transforms (same as original)
  for (iSub in 1:nSub){
    rho[iSub]     = exp(x1[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    epsilon[iSub] = inv_logit(x2[iSub]);

    if (epsilon[iSub] < 0.5){
      biaseps[iSub, 2] = inv_logit(x2[iSub] - x5[iSub]);           // NoGo punishment learning rate
      biaseps[iSub, 1] = 2 * epsilon[iSub] - biaseps[iSub, 2];     // Go reward learning rate
    } else {
      biaseps[iSub, 1] = inv_logit(x2[iSub] + x5[iSub]);           // Go reward learning rate
      biaseps[iSub, 2] = 2 * epsilon[iSub] - biaseps[iSub, 1];     // NoGo punishment learning rate
    }
    betaEEG[iSub] = inv_logit(x6[iSub]); // ∈ [0,1]
  }

  // Generate choice probabilities and update Q online
  for (iSub in 1:nSub){
    Q = Qi * rho[iSub];
    V = Vi;
    for (iStim in 1:nStim) valenced[iStim] = 0;

    for (iTrial in 1:nTrial){
      // Base Q -> q
      for (iResp in 1:nResp){
        q[iResp] = valenced[s[iSub, iTrial]] * Q[iResp, s[iSub, iTrial]];
      }
      // Add Pavlovian contributions on the Go channels
      // (no theta modulation of Pavlovian term in this model)
      for (iResp in 1:2){
        q[iResp] = q[iResp]
                   + gobias[iSub]
                   + pibias[iSub] * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
      }

      // Softmax probabilities
      p0 = softmax(q);
      BGx[(iSub - 1) * nTrial + iTrial] = p0;

      // Prediction error
      er = rho[iSub] * r[iSub, iTrial];

      // Theta × conflict modulation of learning rate (conflict trials only)
      if ( (ya[iSub, iTrial] == 3) && (r[iSub, iTrial] == -1) ){
        // punished NoGo
        if (motivconflict[s[iSub, iTrial]] == 1){
          lr = (1 - betaEEG[iSub] * tPow[iSub, iTrial]) * biaseps[iSub, 2]
               + (betaEEG[iSub] * tPow[iSub, iTrial]) * epsilon[iSub];
        } else {
          lr = biaseps[iSub, 2];
        }
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] += lr * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);

      } else if ( (ya[iSub, iTrial] != 3) && (r[iSub, iTrial] == 1) ){
        // rewarded Go
        if (motivconflict[s[iSub, iTrial]] == 1){
          lr = (1 - betaEEG[iSub] * tPow[iSub, iTrial]) * biaseps[iSub, 1]
               + (betaEEG[iSub] * tPow[iSub, iTrial]) * epsilon[iSub];
        } else {
          lr = biaseps[iSub, 1];
        }
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] += lr * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);

      } else {
        // neutral outcome: use epsilon
        Q[ ya[iSub, iTrial], s[iSub, iTrial] ] += epsilon[iSub] * (er - Q[ ya[iSub, iTrial], s[iSub, iTrial] ]);
      }

      // Once a stimulus receives non-zero feedback, its valence is treated as known
      if (r[iSub, iTrial] != 0)
        valenced[ s[iSub, iTrial] ] = 1;
    }
  }
}

model{
  // Group-level priors (keep original settings)
  X[1]  ~ normal(0, 3); // feedback sensitivity
  X[2]  ~ normal(0, 2); // instrumental learning rate
  X[3]  ~ normal(0, 3); // go bias
  X[4]  ~ normal(0, 3); // pavlovian bias
  X[5]  ~ normal(0, 2); // learning bias
  X[6]  ~ normal(0, 3); // EEG weight
  sdX   ~ cauchy(0, 2); // half-Cauchy

  // Hierarchical structure
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);

  // Likelihood
  for (iData in 1:nData)
    yax[iData] ~ categorical( BGx[iData] );
}

generated quantities{
  vector[nData] log_lik; // for WAIC
  for (iData in 1:nData)
    log_lik[iData] = categorical_lpmf( yax[iData] | BGx[iData] );
}

"



mstr15 = "

// ========================================================================================.
// mstr15 corresponds to M5a in the paper. M5a includes a learning rate (epsilon), 
// feedback sensitivity (rho), go bias, pavlovian bias (pi), learning bias (kappa), and
// dlPFC connectivity effect on the Pavlovian bias - ONLY on conflict trials.
// ========================================================================================.

data{
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];
  int<lower=1>  a[nSub, nTrial];
  int<lower=1>  ya[nSub, nTrial];
  int<lower=-1> r[nSub, nTrial];
  int<lower=1>  rew[nSub, nTrial];
  int<lower=1>  yax[nTrial*nSub];
  int<lower=1>  nData;
  int<lower=1>  motivconflict[nStim];
  real<lower=-1,upper=1> tPow[nSub, nTrial];
  real<lower=-1,upper=1> pfcICPC[nSub, nTrial];
  real<lower=-1,upper=1> lICPC[nSub, nTrial];
  real<lower=-1,upper=1> rICPC[nSub, nTrial];

  matrix[nResp,nStim]   Qi;
  vector[nStim]         Vi;
}

// -------------------------------------------------------------------------.
// parameters (hierarchical)
// -------------------------------------------------------------------------.
parameters{
  // Hierarchical parameters.
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // Subject level parameters.
  vector<lower=-10,  upper=10>[nSub] x1; // rho (in log space before exp).
  vector<lower=-8,   upper=8>[nSub] x2;  // epsilon (in logit space before inv_logit).
  vector<lower=-8,   upper=8>[nSub] x3;  // go bias.
  vector<lower=-8,   upper=8>[nSub] x4;  // pavlovian bias.
  vector<lower=-8,   upper=8>[nSub] x5;  // epsilon bias (kappa).
  vector<lower=-8,   upper=8>[nSub] x6;  // betaEEG (dlPFC connectivity weight).
}

// -------------------------------------------------------------------------.
// transformed parameters: compute subject-level parameters & BGx
// -------------------------------------------------------------------------.
transformed parameters{
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  vector[2] biaseps[nSub];
  real betaEEG[nSub];

  matrix[nResp,nStim] Q;
  vector[nResp]       q;
  vector[nStim]       V;
  simplex[nResp]      p0;
  real                er;
  simplex[nResp]      BGx[nTrial*nSub];
  vector[nStim]       valenced;

  // Retrieve sampled subject parameters
  for (iSub in 1:nSub){
    rho[iSub]     = exp(x1[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    epsilon[iSub] = inv_logit(x2[iSub]);
    if (epsilon[iSub] < 0.5){
      biaseps[iSub,2] = inv_logit(x2[iSub] - x5[iSub]);
      biaseps[iSub,1] = 2*epsilon[iSub] - biaseps[iSub,2];
    } else {
      biaseps[iSub,1] = inv_logit(x2[iSub] + x5[iSub]);
      biaseps[iSub,2] = 2*epsilon[iSub] - biaseps[iSub,1];
    }
    betaEEG[iSub] = x6[iSub];
  }

  // Build BGx by iterating over subjects and trials
  for (iSub in 1:nSub){
    // initial values
    Q = Qi * rho[iSub];
    V = Vi;
    for (iStim in 1:nStim) {
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial){
      // retrieve Q-values for current stimulus
      for (iResp in 1:nResp){
        q[iResp] = valenced[s[iSub,iTrial]] * Q[iResp, s[iSub,iTrial]];
      }

      // Pavlovian term with/without dlPFC connectivity on conflict trials
      if (motivconflict[s[iSub,iTrial]] == 1) {
        for (iResp in 1:2){
          q[iResp] = q[iResp]
                     + gobias[iSub]
                     + (pibias[iSub] + betaEEG[iSub] * pfcICPC[iSub,iTrial])
                       * valenced[s[iSub,iTrial]] * V[s[iSub,iTrial]];
        }
      } else {
        for (iResp in 1:2){
          q[iResp] = q[iResp]
                     + gobias[iSub]
                     + pibias[iSub] * valenced[s[iSub,iTrial]] * V[s[iSub,iTrial]];
        }
      }

      // softmax
      p0 = softmax(q);

      // store into BGx (flattened index over sub × trial)
      {
        int idx;
        idx = (iSub - 1) * nTrial + iTrial;
        BGx[idx] = p0;
      }

      // update Q- and V-values given actual response/outcome
      er = rho[iSub] * r[iSub,iTrial];

      if ( (ya[iSub,iTrial] == 3) && (r[iSub,iTrial] == -1) ) {
        // punished nogo
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + biaseps[iSub,2] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      } else if ( (ya[iSub,iTrial] != 3) && (r[iSub,iTrial] == 1) ) {
        // rewarded go
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + biaseps[iSub,1] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      } else {
        // neutral outcome or other cases
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + epsilon[iSub] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      }

      // update whether valence is known
      if (r[iSub,iTrial] != 0) {
        valenced[s[iSub,iTrial]] = 1;
      }
    } // end iTrial
  } // end iSub
}

// -------------------------------------------------------------------------.
// model block: priors and likelihood
// -------------------------------------------------------------------------.
model{
  // Hierarchical priors
  X[1]  ~ normal(0, 3); // feedback sensitivity.
  X[2]  ~ normal(0, 2); // instrumental learning rate.
  X[3]  ~ normal(0, 3); // go bias.
  X[4]  ~ normal(0, 3); // pavlovian bias.
  X[5]  ~ normal(0, 2); // learning bias.
  X[6]  ~ normal(0, 3); // EEG weight.
  sdX   ~ cauchy(0, 2); // half-Cauchy.

  // Subject-level priors
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);

  // Likelihood
  for (iData in 1:nData){
    yax[iData] ~ categorical(BGx[iData]);
  }
}

// -------------------------------------------------------------------------.
// generated quantities: pointwise log-likelihood
// -------------------------------------------------------------------------.
generated quantities {
  vector[nData] log_lik;

  for (iData in 1:nData){
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}

"



mstr16 = "

// ========================================================================================.
// mstr16 corresponds to M5b in the paper. M5b includes a learning rate (epsilon), 
// feedback sensitivity (rho), go bias, pavlovian bias (pi), learning bias (kappa), and
// motor connectivity effect on the Pavlovian bias - ONLY on conflict trials.
// ========================================================================================.

data{
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];
  int<lower=1>  a[nSub, nTrial];
  int<lower=1>  ya[nSub, nTrial];
  int<lower=-1> r[nSub, nTrial];
  int<lower=1>  rew[nSub, nTrial];
  int<lower=1>  yax[nTrial*nSub];
  int<lower=1>  nData;
  int<lower=1>  motivconflict[nStim];
  real<lower=-1,upper=1> tPow[nSub, nTrial];
  real<lower=-1,upper=1> pfcICPC[nSub, nTrial];
  real<lower=-1,upper=1> lICPC[nSub, nTrial];
  real<lower=-1,upper=1> rICPC[nSub, nTrial];

  matrix[nResp,nStim]   Qi;
  vector[nStim]         Vi;
}

// -------------------------------------------------------------------------.
// parameters (hierarchical)
// -------------------------------------------------------------------------.
parameters{
  // Hierarchical parameters.
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // Subject level parameters.
  vector<lower=-10, upper=10>[nSub] x1; // rho (log space before exp)
  vector<lower=-8,  upper=8>[nSub] x2;  // epsilon (logit space before inv_logit)
  vector<lower=-8,  upper=8>[nSub] x3;  // go bias
  vector<lower=-8,  upper=8>[nSub] x4;  // pavlovian bias
  vector<lower=-8,  upper=8>[nSub] x5;  // epsilon bias (kappa)
  vector<lower=-8,  upper=8>[nSub] x6;  // betaEEG
}

// -------------------------------------------------------------------------.
// transformed parameters: compute subject-level parameters & BGx
// -------------------------------------------------------------------------.
transformed parameters{
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  vector[2] biaseps[nSub];
  real betaEEG[nSub];

  matrix[nResp,nStim] Q;
  vector[nResp]       q;
  vector[nStim]       V;
  simplex[nResp]      p0;
  real                er;
  simplex[nResp]      BGx[nTrial*nSub];
  vector[nStim]       valenced;

  // Retrieve sampled subject parameters
  for (iSub in 1:nSub){
    rho[iSub]     = exp(x1[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    epsilon[iSub] = inv_logit(x2[iSub]);
    if (epsilon[iSub] < 0.5){
      biaseps[iSub,2] = inv_logit(x2[iSub] - x5[iSub]);
      biaseps[iSub,1] = 2*epsilon[iSub] - biaseps[iSub,2];
    } else {
      biaseps[iSub,1] = inv_logit(x2[iSub] + x5[iSub]);
      biaseps[iSub,2] = 2*epsilon[iSub] - biaseps[iSub,1];
    }
    betaEEG[iSub] = x6[iSub];
  }

  // Build BGx by iterating over subjects and trials
  for (iSub in 1:nSub){
    // initial values
    Q = Qi * rho[iSub];
    V = Vi;
    for (iStim in 1:nStim) {
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial){
      // retrieve Q-values for current stimulus
      for (iResp in 1:nResp){
        q[iResp] = valenced[s[iSub,iTrial]] * Q[iResp, s[iSub,iTrial]];
      }

      // Pavlovian term with motor connectivity only on conflict trials
      if (motivconflict[s[iSub,iTrial]] == 1){
        // left response: contralateral right motor cortex (rICPC)
        q[1] = q[1] + gobias[iSub]
                    + (pibias[iSub] + betaEEG[iSub] * rICPC[iSub,iTrial])
                      * valenced[s[iSub,iTrial]] * V[s[iSub,iTrial]];
        // right response: contralateral left motor cortex (lICPC)
        q[2] = q[2] + gobias[iSub]
                    + (pibias[iSub] + betaEEG[iSub] * lICPC[iSub,iTrial])
                      * valenced[s[iSub,iTrial]] * V[s[iSub,iTrial]];
      } else {
        for (iResp in 1:2){
          q[iResp] = q[iResp]
                     + gobias[iSub]
                     + pibias[iSub] * valenced[s[iSub,iTrial]] * V[s[iSub,iTrial]];
        }
      }

      // softmax
      p0 = softmax(q);

      // store into BGx (flattened index over sub × trial)
      {
        int idx;
        idx = (iSub - 1) * nTrial + iTrial;
        BGx[idx] = p0;
      }

      // update Q- and V-values given actual response/outcome
      er = rho[iSub] * r[iSub,iTrial];

      if ( (ya[iSub,iTrial] == 3) && (r[iSub,iTrial] == -1) ) {
        // punished nogo
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + biaseps[iSub,2] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      } else if ( (ya[iSub,iTrial] != 3) && (r[iSub,iTrial] == 1) ) {
        // rewarded go
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + biaseps[iSub,1] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      } else {
        // neutral outcome or other cases
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + epsilon[iSub] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      }

      // update whether valence is known
      if (r[iSub,iTrial] != 0){
        valenced[s[iSub,iTrial]] = 1;
      }
    } // end iTrial
  } // end iSub
}

// -------------------------------------------------------------------------.
// model block: priors and likelihood
// -------------------------------------------------------------------------.
model{
  // Hierarchical priors
  X[1]  ~ normal(0, 3); // feedback sensitivity
  X[2]  ~ normal(0, 2); // instrumental learning rate
  X[3]  ~ normal(0, 3); // go bias
  X[4]  ~ normal(0, 3); // pavlovian bias
  X[5]  ~ normal(0, 2); // learning bias
  X[6]  ~ normal(0, 3); // EEG weight
  sdX   ~ cauchy(0, 2); // half-Cauchy

  // Subject-level priors
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);

  // Likelihood
  for (iData in 1:nData){
    yax[iData] ~ categorical(BGx[iData]);
  }
}

// -------------------------------------------------------------------------.
// generated quantities: pointwise log-likelihood
// -------------------------------------------------------------------------.
generated quantities {
  vector[nData] log_lik;

  for (iData in 1:nData){
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}
"



# ========================================================================================
# mstr17: effect of BOTH theta power (tPow) and midfrontal–dlPFC phase synchrony (pfcICPC)
# on the Pavlovian bias, ONLY on motivational conflict trials.
# Per-subject parameters:
#   x1 = log(rho)
#   x2 = logit(epsilon)
#   x3 = go bias
#   x4 = Pavlovian bias pi
#   x5 = learning bias kappa  → biaseps[ ,1:2]
#   x6 = betaPow   (theta power weight on π, conflict trials)
#   x7 = betaICPC  (PFC ICPC weight on π, conflict trials)
# ========================================================================================
mstr17 = "
data {
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];        // state (cue)
  int<lower=1>  a[nSub, nTrial];        // required response (unused here)
  int<lower=1>  ya[nSub, nTrial];       // chosen response (1..nResp)
  int<lower=-1> r[nSub, nTrial];        // outcome: -1/0/1
  int<lower=1>  rew[nSub, nTrial];      // reward flag (unused here)
  int<lower=1>  yax[nTrial * nSub];     // flattened choices for likelihood
  int<lower=1>  nData;

  // 1 = motivational conflict cue, 2 = non-conflict cue
  int<lower=1, upper=2> motivconflict[nStim];

  // EEG regressors
  real<lower=-1, upper=1> tPow[nSub, nTrial];     // midfrontal theta power
  real<lower=-1, upper=1> pfcICPC[nSub, nTrial];  // midfrontal–dlPFC ICPC
  real<lower=-1, upper=1> lICPC[nSub, nTrial];    // unused in this model
  real<lower=-1, upper=1> rICPC[nSub, nTrial];    // unused in this model

  // initial Q-table and Pavlovian values
  matrix[nResp, nStim] Qi;
  vector[nStim]        Vi;
}

// -------------------------------------------------------------------------
// parameters in 'fitting space'
// -------------------------------------------------------------------------
parameters {
  // Hierarchical (group-level) parameters
  vector[K] X;                          // group means
  vector<lower=0, upper=20>[K] sdX;     // group sds (half-Cauchy prior)

  // Subject-level parameters
  vector<lower=-10, upper=10>[nSub] x1; // log(rho)
  vector<lower=-8,  upper=8>[nSub] x2;  // logit(epsilon)
  vector<lower=-8,  upper=8>[nSub] x3;  // go bias
  vector<lower=-8,  upper=8>[nSub] x4;  // Pavlovian bias pi
  vector<lower=-8,  upper=8>[nSub] x5;  // learning bias kappa
  vector<lower=-8,  upper=8>[nSub] x6;  // betaEEG-power
  vector<lower=-8,  upper=8>[nSub] x7;  // betaEEG-ICPC
}

// -------------------------------------------------------------------------
// transformed parameters: subject params + trial-by-trial choice probs
// -------------------------------------------------------------------------
transformed parameters {
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  matrix[nSub, 2] biaseps;      // [sub,1]=Go-reward rate; [sub,2]=NoGo-punish rate
  real betaPow[nSub];
  real betaICPC[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp] q;
  vector[nStim] V;
  simplex[nResp] p0;
  real er;
  vector[nResp] BGx[nTrial * nSub];
  vector[nStim] valenced;

  // 1) individual-level parameters
  for (iSub in 1:nSub) {
    rho[iSub]     = exp(x1[iSub]);       // feedback sensitivity > 0
    epsilon[iSub] = inv_logit(x2[iSub]); // learning rate in (0,1)
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];

    // learning bias: split epsilon into two biased learning rates
    if (epsilon[iSub] < 0.5) {
      biaseps[iSub, 2] = inv_logit(x2[iSub] - x5[iSub]);       // NoGo-punish
      biaseps[iSub, 1] = 2 * epsilon[iSub] - biaseps[iSub, 2]; // Go-reward
    } else {
      biaseps[iSub, 1] = inv_logit(x2[iSub] + x5[iSub]);       // Go-reward
      biaseps[iSub, 2] = 2 * epsilon[iSub] - biaseps[iSub, 1]; // NoGo-punish
    }

    betaPow[iSub]  = x6[iSub];
    betaICPC[iSub] = x7[iSub];
  }

  // 2) trial-by-trial choice probabilities
  for (iSub in 1:nSub) {
    // initialize Q and Pavlovian values for this subject
    Q = Qi * rho[iSub];
    V = Vi;

    // which stimuli already have known valence?
    for (iStim in 1:nStim) {
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial) {
      int stim;
      stim = s[iSub, iTrial];

      // instrumental values: only if valence is known
      for (iResp in 1:nResp) {
        q[iResp] = valenced[stim] * Q[iResp, stim];
      }

      // add Go-level biases (Go = responses 1 and 2 here)
      for (iResp in 1:2) {
        if (motivconflict[stim] == 1) {
          // conflict trials: Pavlovian bias modulated by both theta power and PFC ICPC
          q[iResp] = q[iResp]
                     + gobias[iSub]
                     + (pibias[iSub]
                        + betaICPC[iSub] * pfcICPC[iSub, iTrial]
                        + betaPow[iSub]  * tPow[iSub, iTrial])
                       * valenced[stim] * V[stim];
        } else {
          // non-conflict trials: baseline Pavlovian bias only
          q[iResp] = q[iResp]
                     + gobias[iSub]
                     + pibias[iSub] * valenced[stim] * V[stim];
        }
      }

      // softmax to obtain choice probabilities
      p0 = softmax(q);

      // store BGx (row-major: subject 1..nSub, each with nTrial rows)
      BGx[(iSub - 1) * nTrial + iTrial] = p0;

      // outcome prediction error
      er = rho[iSub] * r[iSub, iTrial];

      // update Q with biased learning rates
      if ((ya[iSub, iTrial] == 3) && (r[iSub, iTrial] == -1)) {
        // punished NoGo
        Q[ya[iSub, iTrial], stim] =
          Q[ya[iSub, iTrial], stim] +
          biaseps[iSub, 2] * (er - Q[ya[iSub, iTrial], stim]);
      } else if ((ya[iSub, iTrial] != 3) && (r[iSub, iTrial] == 1)) {
        // rewarded Go
        Q[ya[iSub, iTrial], stim] =
          Q[ya[iSub, iTrial], stim] +
          biaseps[iSub, 1] * (er - Q[ya[iSub, iTrial], stim]);
      } else {
        // neutral or other outcomes: use baseline epsilon
        Q[ya[iSub, iTrial], stim] =
          Q[ya[iSub, iTrial], stim] +
          epsilon[iSub] * (er - Q[ya[iSub, iTrial], stim]);
      }

      // update Pavlovian value V (simple TD learning)
      V[stim] =
        V[stim] + epsilon[iSub] * (r[iSub, iTrial] - V[stim]);

      // mark valence as known once a non-zero outcome is observed
      if (r[iSub, iTrial] != 0) {
        valenced[stim] = 1;
      }
    } // end iTrial
  }   // end iSub
}

// -------------------------------------------------------------------------
// model: priors + likelihood
// -------------------------------------------------------------------------
model {
  // Hierarchical priors (group-level)
  X[1] ~ normal(0, 3); // feedback sensitivity (log rho)
  X[2] ~ normal(0, 2); // instrumental learning rate (logit epsilon)
  X[3] ~ normal(0, 3); // go bias
  X[4] ~ normal(0, 3); // Pavlovian bias
  X[5] ~ normal(0, 2); // learning bias (kappa)
  X[6] ~ normal(0, 3); // EEG weight - theta power
  X[7] ~ normal(0, 3); // EEG weight - dlPFC ICPC
  sdX  ~ cauchy(0, 2); // half-Cauchy on SDs

  // Subject-level priors
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);
  x7 ~ normal(X[7], sdX[7]);

  // Likelihood
  for (iData in 1:nData) {
    yax[iData] ~ categorical(BGx[iData]);
  }
}

// -------------------------------------------------------------------------
// generated quantities: pointwise log-likelihood (for WAIC/LOO)
// -------------------------------------------------------------------------
generated quantities {
  vector[nData] log_lik;

  for (iData in 1:nData) {
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}
"




# ========================================================================================
# mstr18: effect of BOTH theta power (tPow) and midfrontal–motor phase synchrony (lICPC/rICPC)
# on the Pavlovian bias, ONLY on motivational conflict trials.
# Parameters per subject:
#   rho (x1), epsilon (x2), go bias (x3), Pavlovian bias pi (x4),
#   learning bias kappa (x5 → biaseps), betaPow (x6), betaICPC (x7)
# ========================================================================================
mstr18 = "
data {
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];        // state (cue)
  int<lower=1>  a[nSub, nTrial];        // required response (unused here)
  int<lower=1>  ya[nSub, nTrial];       // chosen response (1..nResp)
  int<lower=-1> r[nSub, nTrial];        // outcome: -1/0/1
  int<lower=1>  rew[nSub, nTrial];      // reward flag (unused here)
  int<lower=1>  yax[nTrial * nSub];     // flattened choices for likelihood
  int<lower=1>  nData;
  int<lower=1>  motivconflict[nStim];   // 1 = conflict cue (else = non-conflict)

  real<lower=-1, upper=1> tPow[nSub, nTrial];    // midfrontal theta power
  real<lower=-1, upper=1> pfcICPC[nSub, nTrial]; // not used in this model
  real<lower=-1, upper=1> lICPC[nSub, nTrial];   // midfrontal–LEFT motor ICPC
  real<lower=-1, upper=1> rICPC[nSub, nTrial];   // midfrontal–RIGHT motor ICPC

  matrix[nResp, nStim] Qi;   // initial Q table
  vector[nStim]        Vi;   // initial Pavlovian values
}

// -------------------------------------------------------------------------
// parameters in 'fitting space'
// -------------------------------------------------------------------------
parameters {
  // Hierarchical parameters
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // Subject level parameters
  vector<lower=-10, upper=10>[nSub] x1; // rho
  vector<lower=-8,  upper=8>[nSub] x2;  // epsilon
  vector<lower=-8,  upper=8>[nSub] x3;  // go bias
  vector<lower=-8,  upper=8>[nSub] x4;  // Pavlovian bias
  vector<lower=-8,  upper=8>[nSub] x5;  // learning bias (kappa)
  vector<lower=-8,  upper=8>[nSub] x6;  // betaEEG-power
  vector<lower=-8,  upper=8>[nSub] x7;  // betaEEG-ICPC (motor)
}

// -------------------------------------------------------------------------
// transformed parameters: subject params + trial-by-trial choice probs
// -------------------------------------------------------------------------
transformed parameters {
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  matrix[nSub, 2] biaseps;          // [sub,1]=Go-reward; [sub,2]=NoGo-punish
  real betaPow[nSub];
  real betaICPC[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp] q;
  vector[nStim] V;
  simplex[nResp] p0;
  real er;
  vector[nResp] BGx[nTrial * nSub];
  vector[nStim] valenced;

  // 1) individual parameters
  for (iSub in 1:nSub) {
    // retrieve sampled subject parameters and initial values
    rho[iSub]     = exp(x1[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    epsilon[iSub] = inv_logit(x2[iSub]);

    // split epsilon into two biased learning rates
    if (epsilon[iSub] < 0.5) {
      biaseps[iSub, 2] = inv_logit(x2[iSub] - x5[iSub]);       // NoGo-punish
      biaseps[iSub, 1] = 2 * epsilon[iSub] - biaseps[iSub, 2]; // Go-reward
    } else {
      biaseps[iSub, 1] = inv_logit(x2[iSub] + x5[iSub]);       // Go-reward
      biaseps[iSub, 2] = 2 * epsilon[iSub] - biaseps[iSub, 1]; // NoGo-punish
    }

    betaPow[iSub]  = x6[iSub];
    betaICPC[iSub] = x7[iSub];
  }

  // 2) trial-by-trial choice probabilities
  for (iSub in 1:nSub) {
    // retrieve and initialize values
    Q = Qi * rho[iSub];
    V = Vi;

    for (iStim in 1:nStim) {
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial) {
      int idx;
      idx = (iSub - 1) * nTrial + iTrial;

      // instrumental Q-values (only if valence is known)
      for (iResp in 1:nResp) {
        q[iResp] = valenced[s[iSub, iTrial]] * Q[iResp, s[iSub, iTrial]];
      }

      // add Go-level Pavlovian + EEG-modulated biases
      if (motivconflict[s[iSub, iTrial]] == 1) {
        // motivational conflict:
        //   left response (1) ← midfrontal–RIGHT motor ICPC
        //   right response (2) ← midfrontal–LEFT motor ICPC
        q[1] = q[1]
               + gobias[iSub]
               + (pibias[iSub]
                  + betaICPC[iSub] * rICPC[iSub, iTrial]
                  + betaPow[iSub]  * tPow[iSub, iTrial])
                 * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];

        q[2] = q[2]
               + gobias[iSub]
               + (pibias[iSub]
                  + betaICPC[iSub] * lICPC[iSub, iTrial]
                  + betaPow[iSub]  * tPow[iSub, iTrial])
                 * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
      } else {
        // no conflict: standard Pavlovian bias only
        for (iResp in 1:2) {
          q[iResp] = q[iResp]
                     + gobias[iSub]
                     + pibias[iSub] * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
        }
      }

      // retrieve choice probabilities with softmax of Q-values
      p0 = softmax(q);   // probability of each action with softmax

      for (iResp in 1:nResp) {
        BGx[idx, iResp] = p0[iResp];
      }

      // update Q-values (and conceptually V-values)
      er = rho[iSub] * r[iSub, iTrial];  // outcome with feedback sensitivity

      // punished NoGo
      if ((ya[iSub, iTrial] == 3) && (r[iSub, iTrial] == -1)) {
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          biaseps[iSub, 2] * (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      } else if ((ya[iSub, iTrial] != 3) && (r[iSub, iTrial] == 1)) {
        // rewarded Go
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          biaseps[iSub, 1] * (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      } else {
        // neutral / other outcomes
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          epsilon[iSub] * (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      }

      // mark stimulus as \"valence known\" once a non-zero outcome is observed
      if (r[iSub, iTrial] != 0) {
        valenced[s[iSub, iTrial]] = 1;
      }
    } // end iTrial
  }   // end iSub
}

// -------------------------------------------------------------------------
// model: priors + likelihood
// -------------------------------------------------------------------------
model {
  vector[nResp] theta;  // unused, kept for compatibility

  // Hierarchical priors
  X[1] ~ normal(0, 3); // feedback sensitivity
  X[2] ~ normal(0, 2); // instrumental learning rate
  X[3] ~ normal(0, 3); // go bias
  X[4] ~ normal(0, 3); // Pavlovian bias
  X[5] ~ normal(0, 2); // learning bias
  X[6] ~ normal(0, 3); // EEG weight - power
  X[7] ~ normal(0, 3); // EEG weight - motor ICPC
  sdX  ~ cauchy(0, 2); // half-Cauchy on SDs

  // Subject level priors
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);
  x7 ~ normal(X[7], sdX[7]);

  // Likelihood: chosen action given choice probabilities
  for (iData in 1:nData) {
    yax[iData] ~ categorical(BGx[iData]);
  }
}

// -------------------------------------------------------------------------
// generated quantities: pointwise log-likelihood (for WAIC/LOO)
// -------------------------------------------------------------------------
generated quantities {
  vector[nData] log_lik; // To calculate WAIC

  for (iData in 1:nData) {
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}
"




mstr19 = "

// ========================================================================================.
// mstr19 models an alternative effect for dlPFC connectivity, reported in supplementary 
// materials. mstr19 includes a learning rate (epsilon), feedback sensitivity (rho), 
// go bias, pavlovian bias (pi), learning bias (kappa), and dlPFC connectivity effect on 
// the instrumental contribution (Q-values) - ONLY on conflict trials.
// ========================================================================================.

data{
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];
  int<lower=1>  a[nSub, nTrial];
  int<lower=1>  ya[nSub, nTrial];
  int<lower=-1> r[nSub, nTrial];
  int<lower=1>  rew[nSub, nTrial];
  int<lower=1>  yax[nTrial*nSub];
  int<lower=1>  nData;
  int<lower=1>  motivconflict[nStim];
  real<lower=-1,upper=1> tPow[nSub, nTrial];
  real<lower=-1,upper=1> pfcICPC[nSub, nTrial];
  real<lower=-1,upper=1> lICPC[nSub, nTrial];
  real<lower=-1,upper=1> rICPC[nSub, nTrial];

  matrix[nResp,nStim]   Qi;
  vector[nStim]         Vi;
}

// -------------------------------------------------------------------------.
// parameters
// -------------------------------------------------------------------------.
parameters{
  // Hierarchical parameters.
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // Subject level parameters.
  vector<lower=-10, upper=10>[nSub] x1; // rho.
  vector<lower=-8,  upper=8>[nSub] x2;  // epsilon.
  vector<lower=-8,  upper=8>[nSub] x3;  // go bias.
  vector<lower=-8,  upper=8>[nSub] x4;  // pavlovian bias.
  vector<lower=-8,  upper=8>[nSub] x5;  // epsilon bias.
  vector<lower=-8,  upper=8>[nSub] x6;  // betaEEG.
}

// -------------------------------------------------------------------------.
// transformed parameters
// -------------------------------------------------------------------------.
transformed parameters{
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  vector[2] biaseps[nSub];
  real betaEEG[nSub];

  matrix[nResp,nStim] Q;
  vector[nResp]       q;
  vector[nStim]       V;
  simplex[nResp]      p0;
  real                er;
  simplex[nResp]      BGx[nTrial*nSub];
  vector[nStim]       valenced;

  // subject-level parameter transforms
  for (iSub in 1:nSub){
    rho[iSub]     = exp(x1[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    epsilon[iSub] = inv_logit(x2[iSub]);

    if (epsilon[iSub] < 0.5){
      biaseps[iSub,2] = inv_logit(x2[iSub] - x5[iSub]);
      biaseps[iSub,1] = 2*epsilon[iSub] - biaseps[iSub,2];
    } else {
      biaseps[iSub,1] = inv_logit(x2[iSub] + x5[iSub]);
      biaseps[iSub,2] = 2*epsilon[iSub] - biaseps[iSub,1];
    }
    betaEEG[iSub] = x6[iSub];
  }

  // build BGx across subjects and trials
  for (iSub in 1:nSub){
    Q = Qi * rho[iSub];
    V = Vi;
    for (iStim in 1:nStim) {
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial){
      // instrumental contribution from Q
      for (iResp in 1:nResp){
        q[iResp] = valenced[s[iSub,iTrial]] * Q[iResp, s[iSub,iTrial]];
        // dlPFC connectivity scales instrumental contribution ONLY on conflict trials
        if (motivconflict[s[iSub,iTrial]] == 1){
          q[iResp] = (1 - betaEEG[iSub] * pfcICPC[iSub,iTrial]) * q[iResp];
        }
      }

      // add go/pavlovian terms (always added; scaling above is only for conflict)
      for (iResp in 1:2){
        q[iResp] = q[iResp] 
                   + gobias[iSub] 
                   + pibias[iSub] * valenced[s[iSub,iTrial]] * V[s[iSub,iTrial]];
      }

      // softmax
      p0 = softmax(q);

      // store into BGx (flattened index over sub × trial)
      {
        int idx;
        idx = (iSub - 1) * nTrial + iTrial;
        BGx[idx] = p0;
      }

      // update Q given actual response/outcome
      er = rho[iSub] * r[iSub,iTrial];

      if ( (ya[iSub,iTrial] == 3) && (r[iSub,iTrial] == -1) ) {
        // punished NoGo
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + biaseps[iSub,2] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      } else if ( (ya[iSub,iTrial] != 3) && (r[iSub,iTrial] == 1) ) {
        // rewarded Go
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + biaseps[iSub,1] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      } else {
        // neutral outcome / other
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + epsilon[iSub] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      }

      // update whether valence is known
      if (r[iSub,iTrial] != 0){
        valenced[s[iSub,iTrial]] = 1;
      }
    } // iTrial
  } // iSub
}

// -------------------------------------------------------------------------.
// model
// -------------------------------------------------------------------------.
model{
  // hierarchical priors
  X[1]  ~ normal(0, 3);
  X[2]  ~ normal(0, 2);
  X[3]  ~ normal(0, 3);
  X[4]  ~ normal(0, 3);
  X[5]  ~ normal(0, 2);
  X[6]  ~ normal(0, 3);
  sdX   ~ cauchy(0, 2); // half-Cauchy on sd

  // subject-level priors
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);

  // likelihood
  for (iData in 1:nData){
    yax[iData] ~ categorical(BGx[iData]);
  }
}

// -------------------------------------------------------------------------.
// generated quantities
// -------------------------------------------------------------------------.
generated quantities {
  vector[nData] log_lik;
  for (iData in 1:nData){
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}

"



# ========================================================================================
# mstr20: alternative motor-connectivity model (supplementary)
# - Learning rate (epsilon), feedback sensitivity (rho), go bias,
#   Pavlovian bias (pi), learning bias (kappa)
# - Motor connectivity (ICPC) effect on *action weights*,
#   ONLY on motivational conflict trials.
# ========================================================================================
mstr20 = "
data {
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];        // state (cue)
  int<lower=1>  a[nSub, nTrial];        // required response (unused here)
  int<lower=1>  ya[nSub, nTrial];       // chosen response (1..nResp)
  int<lower=-1> r[nSub, nTrial];        // outcome: -1/0/1
  int<lower=1>  rew[nSub, nTrial];      // reward flag (unused here)
  int<lower=1>  yax[nTrial * nSub];     // flattened choices for likelihood
  int<lower=1>  nData;
  int<lower=1>  motivconflict[nStim];   // 1 = conflict cue, else = non-conflict

  real<lower=-1, upper=1> tPow[nSub, nTrial];    // not used here
  real<lower=-1, upper=1> pfcICPC[nSub, nTrial]; // not used here
  real<lower=-1, upper=1> lICPC[nSub, nTrial];   // midfrontal–LEFT motor ICPC
  real<lower=-1, upper=1> rICPC[nSub, nTrial];   // midfrontal–RIGHT motor ICPC

  matrix[nResp, nStim] Qi; // initial Q-table
  vector[nStim]        Vi; // initial Pavlovian values
}

// -------------------------------------------------------------------------
// parameters in 'fitting space'
// -------------------------------------------------------------------------
parameters {
  // Hierarchical parameters
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // Subject level parameters
  vector<lower=-10, upper=10>[nSub] x1; // rho
  vector<lower=-8,  upper=8>[nSub] x2;  // epsilon
  vector<lower=-8,  upper=8>[nSub] x3;  // go bias
  vector<lower=-8,  upper=8>[nSub] x4;  // Pavlovian bias
  vector<lower=-8,  upper=8>[nSub] x5;  // learning bias (kappa → biaseps)
  vector<lower=-8,  upper=8>[nSub] x6;  // betaEEG (motor ICPC effect)
}

// -------------------------------------------------------------------------
// transformed parameters: subject parameters + trial-by-trial choice probs
// -------------------------------------------------------------------------
transformed parameters {
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  matrix[nSub, 2] biaseps;          // [sub,1]=Go-reward; [sub,2]=NoGo-punish
  real betaEEG[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp] q;
  vector[nStim] V;
  simplex[nResp] p0;
  real er;
  vector[nResp] BGx[nTrial * nSub];
  vector[nStim] valenced;

  // 1) individual parameters
  for (iSub in 1:nSub) {
    rho[iSub]     = exp(x1[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    epsilon[iSub] = inv_logit(x2[iSub]);

    // split epsilon into two biased learning rates
    if (epsilon[iSub] < 0.5) {
      biaseps[iSub, 2] = inv_logit(x2[iSub] - x5[iSub]);        // NoGo-punish
      biaseps[iSub, 1] = 2 * epsilon[iSub] - biaseps[iSub, 2];  // Go-reward
    } else {
      biaseps[iSub, 1] = inv_logit(x2[iSub] + x5[iSub]);        // Go-reward
      biaseps[iSub, 2] = 2 * epsilon[iSub] - biaseps[iSub, 1];  // NoGo-punish
    }

    betaEEG[iSub] = x6[iSub];
  }

  // 2) trial-by-trial choice probabilities
  for (iSub in 1:nSub) {
    // initialise Q and Pavlovian values
    Q = Qi * rho[iSub];
    V = Vi;

    for (iStim in 1:nStim) {
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial) {
      int idx;
      idx = (iSub - 1) * nTrial + iTrial;

      // instrumental Q-values (only once valence is known)
      for (iResp in 1:nResp) {
        q[iResp] = valenced[s[iSub, iTrial]] * Q[iResp, s[iSub, iTrial]];
      }

      // add Go-level biases
      if (motivconflict[s[iSub, iTrial]] == 1) {
        // motivational conflict: ICPC adds directly to action weights
        // left response (1) ← midfrontal–RIGHT motor ICPC
        // right response (2) ← midfrontal–LEFT motor ICPC
        q[1] = q[1]
               + gobias[iSub]
               + pibias[iSub] * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]]
               + betaEEG[iSub] * rICPC[iSub, iTrial];

        q[2] = q[2]
               + gobias[iSub]
               + pibias[iSub] * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]]
               + betaEEG[iSub] * lICPC[iSub, iTrial];
      } else {
        // no conflict: standard Pavlovian Go-bias (no ICPC term)
        for (iResp in 1:2) {
          q[iResp] = q[iResp]
                     + gobias[iSub]
                     + pibias[iSub] * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
        }
      }

      // softmax → choice probabilities
      p0 = softmax(q);

      for (iResp in 1:nResp) {
        BGx[idx, iResp] = p0[iResp];
      }

      // update Q-values
      er = rho[iSub] * r[iSub, iTrial];

      // punished NoGo
      if (ya[iSub, iTrial] == 3 && r[iSub, iTrial] == -1) {
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          biaseps[iSub, 2] * (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      } else if (ya[iSub, iTrial] != 3 && r[iSub, iTrial] == 1) {
        // rewarded Go
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          biaseps[iSub, 1] * (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      } else {
        // neutral / other outcomes
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          epsilon[iSub] * (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      }

      // once valence is known (non-zero outcome), mark this cue
      if (r[iSub, iTrial] != 0) {
        valenced[s[iSub, iTrial]] = 1;
      }
    } // end iTrial
  }   // end iSub
}

// -------------------------------------------------------------------------
// model: priors + likelihood
// -------------------------------------------------------------------------
model {
  vector[nResp] theta; // unused, kept for compatibility

  // Hierarchical priors
  X[1] ~ normal(0, 3); // feedback sensitivity
  X[2] ~ normal(0, 2); // instrumental learning rate
  X[3] ~ normal(0, 3); // go bias
  X[4] ~ normal(0, 3); // Pavlovian bias
  X[5] ~ normal(0, 2); // learning bias
  X[6] ~ normal(0, 3); // EEG weight (motor ICPC)
  sdX  ~ cauchy(0, 2); // half-Cauchy on SDs

  // Subject level priors
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);

  // Likelihood
  for (iData in 1:nData) {
    yax[iData] ~ categorical(BGx[iData]);
  }
}

// -------------------------------------------------------------------------
// generated quantities: pointwise log-likelihood (for WAIC/LOO)
// -------------------------------------------------------------------------
generated quantities {
  vector[nData] log_lik;

  for (iData in 1:nData) {
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}
"





mstr21 = "

// ========================================================================================.
// mstr21: Motor connectivity (lICPC/rICPC) scales the instrumental Q contribution on
// conflict trials (side-specific), then adds Go/Pav biases as usual (like mstr19).
// Combines the 'motor' source of mstr16 with the 'scale-Q-on-conflict' mechanism of mstr19.
// ========================================================================================.

data{
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];
  int<lower=1>  a[nSub, nTrial];
  int<lower=1>  ya[nSub, nTrial];
  int<lower=-1> r[nSub, nTrial];
  int<lower=1>  rew[nSub, nTrial];
  int<lower=1>  yax[nTrial*nSub];
  int<lower=1>  nData;
  int<lower=1>  motivconflict[nStim];                // 1=conflict, else=non-conflict
  real<lower=-1,upper=1> tPow[nSub, nTrial];         // kept for interface parity
  real<lower=-1,upper=1> pfcICPC[nSub, nTrial];      // kept for interface parity (unused)
  real<lower=-1,upper=1> lICPC[nSub, nTrial];
  real<lower=-1,upper=1> rICPC[nSub, nTrial];

  matrix[nResp,nStim]   Qi;
  vector[nStim]         Vi;
}

// -------------------------------------------------------------------------.
// parameters (hierarchical)
// -------------------------------------------------------------------------.
parameters{
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  vector<lower=-10, upper=10>[nSub] x1; // rho (log space)
  vector<lower=-8,  upper=8>[nSub] x2;  // epsilon (logit space)
  vector<lower=-8,  upper=8>[nSub] x3;  // go bias
  vector<lower=-8,  upper=8>[nSub] x4;  // pavlovian bias
  vector<lower=-8,  upper=8>[nSub] x5;  // epsilon bias (kappa)
  vector<lower=-8,  upper=8>[nSub] x6;  // betaEEG (weight for connectivity)
}

// -------------------------------------------------------------------------.
// transformed parameters: subject-level params & BGx
// -------------------------------------------------------------------------.
transformed parameters{
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  vector[2] biaseps[nSub];
  real betaEEG[nSub];

  matrix[nResp,nStim] Q;
  vector[nResp]       q;
  vector[nStim]       V;
  simplex[nResp]      p0;
  real                er;
  simplex[nResp]      BGx[nTrial*nSub];
  vector[nStim]       valenced;

  // transforms
  for (iSub in 1:nSub){
    rho[iSub]     = exp(x1[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    epsilon[iSub] = inv_logit(x2[iSub]);

    if (epsilon[iSub] < 0.5){
      biaseps[iSub,2] = inv_logit(x2[iSub] - x5[iSub]);
      biaseps[iSub,1] = 2*epsilon[iSub] - biaseps[iSub,2];
    } else {
      biaseps[iSub,1] = inv_logit(x2[iSub] + x5[iSub]);
      biaseps[iSub,2] = 2*epsilon[iSub] - biaseps[iSub,1];
    }
    betaEEG[iSub] = x6[iSub];
  }

  // forward loop for BGx
  for (iSub in 1:nSub){
    Q = Qi * rho[iSub];
    V = Vi;
    for (iStim in 1:nStim) {
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial){
      // base instrumental contribution (pre-bias)
      for (iResp in 1:nResp){
        q[iResp] = valenced[s[iSub,iTrial]] * Q[iResp, s[iSub,iTrial]];
      }

      // ONLY on conflict trials: side-specific motor connectivity scales Q
      if (motivconflict[s[iSub,iTrial]] == 1){
        // left response uses contralateral RIGHT motor connectivity (rICPC)
        q[1] = (1 - betaEEG[iSub] * rICPC[iSub,iTrial]) * q[1];
        // right response uses contralateral LEFT motor connectivity (lICPC)
        q[2] = (1 - betaEEG[iSub] * lICPC[iSub,iTrial]) * q[2];
        // NoGo (q[3]) is left unchanged (matches mstr16's Go-only modulation idea)
      }

      // add Go + Pavlovian terms (as in mstr19/mstr16)
      for (iResp in 1:2){
        q[iResp] = q[iResp]
                   + gobias[iSub]
                   + pibias[iSub] * valenced[s[iSub,iTrial]] * V[s[iSub,iTrial]];
      }

      // softmax & store
      p0 = softmax(q);
      {
        int idx;
        idx = (iSub - 1) * nTrial + iTrial;
        BGx[idx] = p0;
      }

      // update Q given actual response/outcome (same learning logic)
      er = rho[iSub] * r[iSub,iTrial];

      if ( (ya[iSub,iTrial] == 3) && (r[iSub,iTrial] == -1) ) {
        // punished NoGo
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + biaseps[iSub,2] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      } else if ( (ya[iSub,iTrial] != 3) && (r[iSub,iTrial] == 1) ) {
        // rewarded Go
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + biaseps[iSub,1] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      } else {
        // neutral / other
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + epsilon[iSub] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      }

      if (r[iSub,iTrial] != 0){
        valenced[s[iSub,iTrial]] = 1;
      }
    } // iTrial
  } // iSub
}

// -------------------------------------------------------------------------.
// model: priors and likelihood
// -------------------------------------------------------------------------.
model{
  // hierarchical priors
  X[1]  ~ normal(0, 3);
  X[2]  ~ normal(0, 2);
  X[3]  ~ normal(0, 3);
  X[4]  ~ normal(0, 3);
  X[5]  ~ normal(0, 2);
  X[6]  ~ normal(0, 3);
  sdX   ~ cauchy(0, 2);  // half-Cauchy on sd

  // subject-level priors
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);

  // likelihood
  for (iData in 1:nData){
    yax[iData] ~ categorical(BGx[iData]);
  }
}

// -------------------------------------------------------------------------.
// generated quantities
// -------------------------------------------------------------------------.
generated quantities {
  vector[nData] log_lik;
  for (iData in 1:nData){
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}

"





# ========================================================================================
# mstr22: EEG control model for M4b
# - Same as M4b (mstr12): includes epsilon, rho, go bias, pavlovian bias,
#   learning bias (kappa), and modulation of instrumental Q values by
#   midfrontal theta (betaEEG)
# - Difference from mstr12: the scaling effect of betaEEG * tPow on Q is applied
#   on all trials, rather than only on conflict cues with motivconflict == 1.
# ========================================================================================
mstr22 = "
data {
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];        // cue/state
  int<lower=1>  a[nSub, nTrial];        // required response (unused here)
  int<lower=1>  ya[nSub, nTrial];       // chosen response (1..nResp)
  int<lower=-1> r[nSub, nTrial];        // outcome: -1/0/1
  int<lower=1>  rew[nSub, nTrial];      // reward flag (unused here)
  int<lower=1>  yax[nTrial * nSub];     // flattened choices for likelihood
  int<lower=1>  nData;

  int<lower=1>  motivconflict[nStim];   // retained for compatibility with other models; unused here
  real<lower=-1, upper=1> tPow[nSub, nTrial];
  real<lower=-1, upper=1> pfcICPC[nSub, nTrial]; // unused
  real<lower=-1, upper=1> lICPC[nSub, nTrial];   // unused
  real<lower=-1, upper=1> rICPC[nSub, nTrial];   // unused

  matrix[nResp, nStim] Qi;              // initial Q-table
  vector[nStim]        Vi;              // initial Pavlovian values
}

parameters {
  // Hierarchical (group-level) parameters
  vector[K] X;                          // [rho, epsilon, gobias, pibias, kappa, betaEEG]
  vector<lower=0, upper=20>[K] sdX;     // group SDs

  // Subject-level parameters
  vector<lower=-10, upper=10>[nSub] x1; // log(rho)
  vector<lower=-8,  upper=8>[nSub] x2;  // logit(epsilon)
  vector<lower=-8,  upper=8>[nSub] x3;  // go bias
  vector<lower=-8,  upper=8>[nSub] x4;  // pavlovian bias
  vector<lower=-8,  upper=8>[nSub] x5;  // learning bias (kappa)
  vector<lower=-8,  upper=8>[nSub] x6;  // betaEEG
}

transformed parameters {
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  matrix[nSub, 2] biaseps;              // [ ,1]=Go-reward; [ ,2]=NoGo-punish
  real betaEEG[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp] q;
  vector[nStim] V;
  simplex[nResp] p0;
  real er;
  vector[nResp] BGx[nTrial * nSub];
  vector[nStim] valenced;

  // ----- individual parameters -----
  for (iSub in 1:nSub) {
    rho[iSub]     = exp(x1[iSub]);       // > 0
    epsilon[iSub] = inv_logit(x2[iSub]); // (0,1)
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    betaEEG[iSub] = x6[iSub];

    // Split learning rate into Go-reward vs NoGo-punish
    // (kept consistent with the M5a/M3c family)
    if (epsilon[iSub] < 0.5) {
      biaseps[iSub, 2] = inv_logit(x2[iSub] - x5[iSub]);           // punished NoGo
      biaseps[iSub, 1] = 2 * epsilon[iSub] - biaseps[iSub, 2];     // rewarded Go
    } else {
      biaseps[iSub, 1] = inv_logit(x2[iSub] + x5[iSub]);           // rewarded Go
      biaseps[iSub, 2] = 2 * epsilon[iSub] - biaseps[iSub, 1];     // punished NoGo
    }
  }

  // ----- trial-by-trial choice probabilities -----
  for (iSub in 1:nSub) {
    Q = Qi * rho[iSub];
    V = Vi;

    for (iStim in 1:nStim) {
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial) {
      int idx;
      idx = (iSub - 1) * nTrial + iTrial;

      // 1) Start with instrumental Q values
      //    (only when cue valence is known)
      for (iResp in 1:nResp) {
        q[iResp] = valenced[s[iSub, iTrial]] * Q[iResp, s[iSub, iTrial]];
      }

      // 2) Scale the Q contribution by midfrontal theta on all trials
      //    In original M4b/mstr12, scaling applies only when motivconflict == 1;
      //    here it is applied to every trial using (1 - betaEEG * tPow).
      for (iResp in 1:nResp) {
        q[iResp] = (1 - betaEEG[iSub] * tPow[iSub, iTrial]) * q[iResp];
      }

      // 3) Add Go bias and Pavlovian term
      //    (theta no longer directly modulates these terms)
      for (iResp in 1:2) {  // only Go options
        q[iResp] = q[iResp]
                   + gobias[iSub]
                   + pibias[iSub] * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
      }

      // 4) Softmax choice probabilities
      p0      = softmax(q);
      BGx[idx] = p0;

      // 5) Update Q and V
      er = rho[iSub] * r[iSub, iTrial];  // scaled outcome

      if ((ya[iSub, iTrial] == 3) && (r[iSub, iTrial] == -1)) {
        // punished NoGo
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          biaseps[iSub, 2] *
          (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      } else if ((ya[iSub, iTrial] != 3) && (r[iSub, iTrial] == 1)) {
        // rewarded Go
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          biaseps[iSub, 1] *
          (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      } else {
        // neutral or other outcomes
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]] +
          epsilon[iSub] *
          (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      }

      // Pavlovian value learning
      V[s[iSub, iTrial]] =
        V[s[iSub, iTrial]] +
        epsilon[iSub] * (r[iSub, iTrial] - V[s[iSub, iTrial]]);

      // Once a non-zero outcome occurs, mark the cue as valenced
      if (r[iSub, iTrial] != 0)
        valenced[s[iSub, iTrial]] = 1;
    }
  }
}

model {
  // Group-level priors
  X[1] ~ normal(0, 3);  // rho
  X[2] ~ normal(0, 2);  // epsilon
  X[3] ~ normal(0, 3);  // go bias
  X[4] ~ normal(0, 3);  // pavlovian bias
  X[5] ~ normal(0, 2);  // learning bias
  X[6] ~ normal(0, 3);  // EEG weight
  sdX  ~ cauchy(0, 2);  // half-Cauchy on SDs

  // Subject-level priors
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);

  // Likelihood
  for (iData in 1:nData) {
    yax[iData] ~ categorical(BGx[iData]);
  }
}

generated quantities {
  vector[nData] log_lik;

  for (iData in 1:nData) {
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}
"





# ========================================================================================.
# mstr23: theta power (tPow) + dlPFC ICPC (pfcICPC) jointly modulate Q-values
# ONLY on motivational conflict trials.
# 
# Per-subject parameters:
#   x1: rho          (feedback sensitivity, >0)
#   x2: epsilon      (learning rate, 0..1)
#   x3: go bias
#   x4: pavlovian bias pi
#   x5: learning bias kappa (→ biaseps)
#   x6: betaTheta    (theta power effect on Q)
#   x7: betaPFC      (dlPFC ICPC effect on Q)
#   - mstr23: tPow + pfcICPC on instrumental Q
# ========================================================================================.
mstr23 = "

data{
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];
  int<lower=1>  a[nSub, nTrial];
  int<lower=1>  ya[nSub, nTrial];
  int<lower=-1> r[nSub, nTrial];
  int<lower=1>  rew[nSub, nTrial];
  int<lower=1>  yax[nTrial*nSub];
  int<lower=1>  nData;
  int<lower=1>  motivconflict[nStim];          // 1 = conflict, else = non-conflict

  real<lower=-1,upper=1> tPow[nSub, nTrial];   // midfrontal theta power
  real<lower=-1,upper=1> pfcICPC[nSub, nTrial];// dlPFC ICPC
  real<lower=-1,upper=1> lICPC[nSub, nTrial];  // unused here
  real<lower=-1,upper=1> rICPC[nSub, nTrial];  // unused here

  matrix[nResp,nStim]   Qi;
  vector[nStim]         Vi;
}

parameters{

  // Hierarchical parameters.
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // Subject level parameters.
  vector<lower=-10, upper=10>[nSub] x1; // rho (log space).
  vector<lower=-8,  upper=8>[nSub] x2;  // epsilon (logit space).
  vector<lower=-8,  upper=8>[nSub] x3;  // go bias.
  vector<lower=-8,  upper=8>[nSub] x4;  // pavlovian bias pi.
  vector<lower=-8,  upper=8>[nSub] x5;  // learning bias kappa.
  vector<lower=-8,  upper=8>[nSub] x6;  // betaTheta (theta effect on Q).
  vector<lower=-8,  upper=8>[nSub] x7;  // betaPFC   (PFC ICPC effect on Q).
}

transformed parameters{

  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  vector[2] biaseps[nSub];
  real betaTheta[nSub];
  real betaPFC[nSub];

  matrix[nResp,nStim] Q;
  vector[nResp]       q;
  vector[nStim]       V;
  simplex[nResp]      p0;
  real                er;
  vector[nResp]       BGx[nTrial*nSub];
  vector[nStim]       valenced;

  // --------- 1) subject-level parameters ---------
  for (iSub in 1:nSub){
    rho[iSub]     = exp(x1[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    epsilon[iSub] = inv_logit(x2[iSub]);

    // learning bias: split epsilon into two biased learning rates
    if (epsilon[iSub] < 0.5){
      biaseps[iSub,2] = inv_logit(x2[iSub] - x5[iSub]);       // NoGo-punish
      biaseps[iSub,1] = 2*epsilon[iSub] - biaseps[iSub,2];    // Go-reward
    } else {
      biaseps[iSub,1] = inv_logit(x2[iSub] + x5[iSub]);       // Go-reward
      biaseps[iSub,2] = 2*epsilon[iSub] - biaseps[iSub,1];    // NoGo-punish
    }

    betaTheta[iSub] = x6[iSub];
    betaPFC[iSub]   = x7[iSub];
  }

  // --------- 2) build BGx over subjects × trials ---------
  for (iSub in 1:nSub){

    // initial values
    Q = Qi * rho[iSub];
    V = Vi;
    for (iStim in 1:nStim){
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial){

      // 2.1 instrumental Q-values for current stimulus
      for (iResp in 1:nResp){
        q[iResp] = valenced[s[iSub,iTrial]] * Q[iResp, s[iSub,iTrial]];

        // ONLY on conflict trials: joint theta + PFC modulation of Q
        if (motivconflict[s[iSub,iTrial]] == 1){
          q[iResp] = ( 1
                       - betaTheta[iSub] * tPow[iSub,iTrial]
                       - betaPFC[iSub]   * pfcICPC[iSub,iTrial]
                     ) * q[iResp];
        }
      }

      // 2.2 add Go-level bias + Pavlovian term (same for conflict / non-conflict)
      for (iResp in 1:2){
        q[iResp] = q[iResp]
                   + gobias[iSub]
                   + pibias[iSub] * valenced[s[iSub,iTrial]] * V[s[iSub,iTrial]];
      }

      // 2.3 softmax → choice probabilities
      p0 = softmax(q);

      {
        int idx;
        idx = (iSub - 1) * nTrial + iTrial;
        BGx[idx] = p0;
      }

      // 2.4 update Q and V using actual choice & outcome
      er = rho[iSub] * r[iSub,iTrial];

      if ( (ya[iSub,iTrial] == 3) && (r[iSub,iTrial] == -1) ){
        // punished NoGo
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + biaseps[iSub,2] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      } else if ( (ya[iSub,iTrial] != 3) && (r[iSub,iTrial] == 1) ){
        // rewarded Go
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + biaseps[iSub,1] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      } else {
        // neutral / other outcomes
        Q[ya[iSub,iTrial], s[iSub,iTrial]] =
          Q[ya[iSub,iTrial], s[iSub,iTrial]]
          + epsilon[iSub] * (er - Q[ya[iSub,iTrial], s[iSub,iTrial]]);
      }

      // 2.5 mark cue as valence-known once non-zero outcome observed
      if (r[iSub,iTrial] != 0){
        valenced[s[iSub,iTrial]] = 1;
      }

    } // end iTrial
  }   // end iSub
}

model{
  // Hierarchical priors.
  X[1] ~ normal(0, 3); // log rho
  X[2] ~ normal(0, 2); // logit epsilon
  X[3] ~ normal(0, 3); // go bias
  X[4] ~ normal(0, 3); // pavlovian bias pi
  X[5] ~ normal(0, 2); // learning bias kappa
  X[6] ~ normal(0, 3); // betaTheta
  X[7] ~ normal(0, 3); // betaPFC
  sdX  ~ cauchy(0, 2); // half-Cauchy

  // Subject-level priors.
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);
  x7 ~ normal(X[7], sdX[7]);

  // Likelihood: all trials contribute, but conflict trials have EEG-modulated Q.
  for (iData in 1:nData){
    yax[iData] ~ categorical(BGx[iData]);
  }
}

generated quantities{
  vector[nData] log_lik;

  for (iData in 1:nData){
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}

"  # end mstr23







# ========================================================================================.
# mstr24: theta power (tPow) + motor ICPC (lICPC/rICPC) jointly modulate instrumental
# Q-values on motivational conflict trials.
#
# Per-subject parameters:
#   x1: rho          (feedback sensitivity, >0)
#   x2: epsilon      (learning rate, 0..1)
#   x3: go bias
#   x4: pavlovian bias pi
#   x5: learning bias kappa (→ biaseps)
#   x6: betaTheta    (theta power effect on Q)
#   x7: betaMotor    (motor ICPC effect on Q, side-specific)
# ========================================================================================.
mstr24 = "
data {
  int<lower=1>  nTrial;
  int<lower=1>  nSub;
  int<lower=1>  nStim;
  int<lower=1>  nResp;
  int<lower=1>  K;
  int<lower=1>  s[nSub, nTrial];
  int<lower=1>  a[nSub, nTrial];
  int<lower=1>  ya[nSub, nTrial];
  int<lower=-1> r[nSub, nTrial];
  int<lower=1>  rew[nSub, nTrial];
  int<lower=1>  yax[nTrial * nSub];
  int<lower=1>  nData;
  int<lower=1>  motivconflict[nStim];          // 1 = conflict, else = non-conflict

  real<lower=-1, upper=1> tPow[nSub, nTrial];  // midfrontal theta power
  real<lower=-1, upper=1> pfcICPC[nSub, nTrial]; // unused, kept for interface parity
  real<lower=-1, upper=1> lICPC[nSub, nTrial];  // left motor ICPC
  real<lower=-1, upper=1> rICPC[nSub, nTrial];  // right motor ICPC

  matrix[nResp, nStim]   Qi;
  vector[nStim]          Vi;
}

parameters {
  // Hierarchical parameters.
  vector[K] X;
  vector<lower=0, upper=20>[K] sdX;

  // Subject level parameters.
  vector<lower=-10, upper=10>[nSub] x1; // rho (log space).
  vector<lower=-8,  upper=8>[nSub] x2;  // epsilon (logit space).
  vector<lower=-8,  upper=8>[nSub] x3;  // go bias.
  vector<lower=-8,  upper=8>[nSub] x4;  // pavlovian bias pi.
  vector<lower=-8,  upper=8>[nSub] x5;  // epsilon bias kappa.
  vector<lower=-8,  upper=8>[nSub] x6;  // betaTheta (theta effect on Q).
  vector<lower=-8,  upper=8>[nSub] x7;  // betaMotor (motor ICPC effect on Q).
}

transformed parameters {
  real rho[nSub];
  real<lower=0, upper=1> epsilon[nSub];
  real gobias[nSub];
  real pibias[nSub];
  vector[2] biaseps[nSub];
  real betaTheta[nSub];
  real betaMotor[nSub];

  matrix[nResp, nStim] Q;
  vector[nResp]        q;
  vector[nStim]        V;
  simplex[nResp]       p0;
  real                 er;
  vector[nResp]        BGx[nTrial * nSub];
  vector[nStim]        valenced;

  // 1) subject-level transforms
  for (iSub in 1:nSub) {
    rho[iSub]     = exp(x1[iSub]);
    gobias[iSub]  = x3[iSub];
    pibias[iSub]  = x4[iSub];
    epsilon[iSub] = inv_logit(x2[iSub]);

    if (epsilon[iSub] < 0.5) {
      biaseps[iSub, 2] = inv_logit(x2[iSub] - x5[iSub]);      // NoGo-punish
      biaseps[iSub, 1] = 2 * epsilon[iSub] - biaseps[iSub, 2]; // Go-reward
    } else {
      biaseps[iSub, 1] = inv_logit(x2[iSub] + x5[iSub]);      // Go-reward
      biaseps[iSub, 2] = 2 * epsilon[iSub] - biaseps[iSub, 1]; // NoGo-punish
    }

    betaTheta[iSub] = x6[iSub];
    betaMotor[iSub] = x7[iSub];
  }

  // 2) forward loop: build BGx
  for (iSub in 1:nSub) {

    // initial values
    Q = Qi * rho[iSub];
    V = Vi;
    for (iStim in 1:nStim) {
      valenced[iStim] = 0;
    }

    for (iTrial in 1:nTrial) {

      // 2.1 instrumental Q-values (before EEG modulation)
      for (iResp in 1:nResp) {
        q[iResp] = valenced[s[iSub, iTrial]] * Q[iResp, s[iSub, iTrial]];
      }

      // 2.2 ONLY on conflict trials: theta + motor ICPC jointly scale Q for Go responses
      if (motivconflict[s[iSub, iTrial]] == 1) {

        // left Go (response 1) uses contralateral RIGHT motor ICPC (rICPC)
        q[1] = (1
                - betaTheta[iSub] * tPow[iSub, iTrial]
                - betaMotor[iSub] * rICPC[iSub, iTrial])
               * q[1];

        // right Go (response 2) uses contralateral LEFT motor ICPC (lICPC)
        q[2] = (1
                - betaTheta[iSub] * tPow[iSub, iTrial]
                - betaMotor[iSub] * lICPC[iSub, iTrial])
               * q[2];

        // NoGo (response 3) left unchanged
      }

      // 2.3 add go-bias + Pavlovian bias (same as mstr12/mstr21 family)
      for (iResp in 1:2) {
        q[iResp] = q[iResp]
                   + gobias[iSub]
                   + pibias[iSub] * valenced[s[iSub, iTrial]] * V[s[iSub, iTrial]];
      }

      // 2.4 softmax to get choice probabilities
      p0 = softmax(q);

      {
        int idx;
        idx = (iSub - 1) * nTrial + iTrial;
        BGx[idx] = p0;
      }

      // 2.5 update Q- and V-values (same learning rule as mstr12/18)
      er = rho[iSub] * r[iSub, iTrial];

      if ((ya[iSub, iTrial] == 3) && (r[iSub, iTrial] == -1)) {
        // punished NoGo
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]]
          + biaseps[iSub, 2] * (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      } else if ((ya[iSub, iTrial] != 3) && (r[iSub, iTrial] == 1)) {
        // rewarded Go
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]]
          + biaseps[iSub, 1] * (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      } else {
        // neutral / other outcomes
        Q[ya[iSub, iTrial], s[iSub, iTrial]] =
          Q[ya[iSub, iTrial], s[iSub, iTrial]]
          + epsilon[iSub] * (er - Q[ya[iSub, iTrial], s[iSub, iTrial]]);
      }

      // 2.6 mark cue as valence-known once non-zero outcome observed
      if (r[iSub, iTrial] != 0) {
        valenced[s[iSub, iTrial]] = 1;
      }

    } // end iTrial
  }   // end iSub
}

model {
  // Hierarchical priors.
  X[1] ~ normal(0, 3); // log rho
  X[2] ~ normal(0, 2); // logit epsilon
  X[3] ~ normal(0, 3); // go bias
  X[4] ~ normal(0, 3); // pavlovian bias pi
  X[5] ~ normal(0, 2); // learning bias kappa
  X[6] ~ normal(0, 3); // betaTheta
  X[7] ~ normal(0, 3); // betaMotor
  sdX  ~ cauchy(0, 2); // half-Cauchy

  // Subject-level priors.
  x1 ~ normal(X[1], sdX[1]);
  x2 ~ normal(X[2], sdX[2]);
  x3 ~ normal(X[3], sdX[3]);
  x4 ~ normal(X[4], sdX[4]);
  x5 ~ normal(X[5], sdX[5]);
  x6 ~ normal(X[6], sdX[6]);
  x7 ~ normal(X[7], sdX[7]);

  // Likelihood.
  for (iData in 1:nData) {
    yax[iData] ~ categorical(BGx[iData]);
  }
}

generated quantities {
  vector[nData] log_lik;

  for (iData in 1:nData) {
    log_lik[iData] = categorical_lpmf(yax[iData] | BGx[iData]);
  }
}
"
