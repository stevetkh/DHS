data{
  int<lower=0> k; //marginal basis dimension
  int<lower=0> tips;
  int<lower=0> n_m; //number of data
  int<lower=0> n_f; //number of data
  int<lower=0> dx_m[n_m];
  vector<lower=0>[n_m] Ex_m;
  int<lower=0> dx_f[n_f];
  vector<lower=0>[n_f] Ex_f;
  int<lower=0> tp_m[n_m];
  int<lower=0> tp_f[n_f];
  matrix[n_m,k*k] DX_m; //full design matrix
  matrix[n_f,k*k] DX_f;
  matrix[k*k,k*k] penal_age;
  matrix[k*k,k*k] penal_time;
  matrix[k*k,k*k] null_penal;
  matrix[k*k,k*k] full_mat_m;
  matrix[k*k,k*k] full_mat_f;
  vector[k*k] XWz_m;
  vector[k*k] XWz_f;
}


parameters{
  vector[k*k] work_params_m;
  vector[k*k] work_params_f;
  vector<lower=-30,upper=30>[tips] tips_params;
  real<lower=-30,upper=20> log_lambda_age_m;
  real<lower=-30,upper=20> log_lambda_time_m;
  real<lower=-30,upper=20> log_lambda_age_f;
  real<lower=-30,upper=20> log_lambda_time_f;
}

transformed parameters{
vector[k*k] spline_params_m;
vector[k*k] spline_params_f;
real log_jac_m;
real log_jac_f;
{
matrix[k*k,k*k] full_mat_L_m;
matrix[k*k,k*k] full_mat_L_f;
vector[n_m] tp_ind_m;
vector[n_f] tp_ind_f;

for(i in 1:n_m){
tp_ind_m[i] = tp_m[i]>2 ? tips_params[tp_m[i]] : tp_m[i]==2 ? 0 : tips_params[tp_m[i]+1];
//tp_ind_m[i] = tp_m[i]>0 ? tips_params[tp_m[i]] : 0;
}

for(i in 1:n_f){
tp_ind_f[i] = tp_f[i]>2 ? tips_params[tp_f[i]] : tp_f[i]==2 ? 0 : tips_params[tp_f[i]+1];
//tp_ind_f[i] = tp_f[i]>0 ? tips_params[tp_f[i]] : 0;
}

full_mat_L_m = cholesky_decompose(full_mat_m+exp(log_lambda_age_m)*penal_age+exp(log_lambda_time_m)*penal_time+null_penal);
spline_params_m = mdivide_right_tri_low(work_params_m',full_mat_L_m)'+ mdivide_right_tri_low(mdivide_left_tri_low(full_mat_L_m,XWz_m - DX_m' *((to_vector(dx_m)+0.1) .* tp_ind_m))',full_mat_L_m)';
full_mat_L_f = cholesky_decompose(full_mat_f+exp(log_lambda_age_f)*penal_age+exp(log_lambda_time_f)*penal_time+null_penal);
spline_params_f = mdivide_right_tri_low(work_params_f',full_mat_L_f)'+ mdivide_right_tri_low(mdivide_left_tri_low(full_mat_L_f,XWz_f - DX_f' *((to_vector(dx_f)+0.1) .* tp_ind_f))',full_mat_L_f)';

log_jac_m = sum(log(diagonal(full_mat_L_m)));
log_jac_f = sum(log(diagonal(full_mat_L_f)));
}
}

model{
vector[n_m] tp_ind_m;
vector[n_f] tp_ind_f;

spline_params_m ~ multi_normal_prec(rep_vector(0,k*k),exp(log_lambda_age_m)*penal_age+exp(log_lambda_time_m)*penal_time+null_penal);
spline_params_f ~ multi_normal_prec(rep_vector(0,k*k),exp(log_lambda_age_f)*penal_age+exp(log_lambda_time_f)*penal_time+null_penal);

target += -log_jac_m;
target += -log_jac_f;

for(i in 1:n_m){
tp_ind_m[i] = tp_m[i]>2 ? tips_params[tp_m[i]] : tp_m[i]==2 ? 0 : tips_params[tp_m[i]+1];
//tp_ind_m[i] = tp_m[i]>0 ? tips_params[tp_m[i]] : 0;
}

for(i in 1:n_f){
tp_ind_f[i] = tp_f[i]>2 ? tips_params[tp_f[i]] : tp_f[i]==2 ? 0 : tips_params[tp_f[i]+1];
//tp_ind_f[i] = tp_f[i]>0 ? tips_params[tp_f[i]] : 0;
}

dx_m ~ poisson(Ex_m .* exp(DX_m*spline_params_m+tp_ind_m));
dx_f ~ poisson(Ex_f .* exp(DX_f*spline_params_f+tp_ind_f));
}