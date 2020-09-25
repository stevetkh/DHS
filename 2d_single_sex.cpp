#include <TMB.hpp>                                

template<class Type>
Type objective_function<Type>::operator() ()
{
using namespace Eigen;
using namespace density;

DATA_VECTOR(dm);
DATA_VECTOR(Em);
DATA_VECTOR(df);
DATA_VECTOR(Ef);
DATA_MATRIX(Dm);
DATA_MATRIX(Df);
DATA_IVECTOR(tpm);
DATA_IVECTOR(tpf);
DATA_MATRIX(penal_age);
DATA_MATRIX(penal_time);
DATA_MATRIX(null_penal);
PARAMETER(log_lambda_age_m);
PARAMETER(log_lambda_time_m);
PARAMETER(log_lambda_age_f);
PARAMETER(log_lambda_time_f);
PARAMETER_VECTOR(tips_params);
PARAMETER_VECTOR(spline_params_m);
PARAMETER_VECTOR(spline_params_f);


matrix<Type> QQm = exp(log_lambda_age_m)*penal_age+exp(log_lambda_time_m)*penal_time+null_penal;
SparseMatrix<Type> full_penal_m = asSparseMatrix(QQm);

matrix<Type> QQf = exp(log_lambda_age_f)*penal_age+exp(log_lambda_time_f)*penal_time+null_penal;
SparseMatrix<Type> full_penal_f = asSparseMatrix(QQf);

vector<Type> tipsm(tpm.size());
vector<Type> tipsf(tpf.size());

for(int i=0; i < tpm.size() ; i++){
tipsm(i) = tpm(i)>2 ? tips_params(tpm(i)-1) : tpm(i)==2 ? 0 : tips_params(tpm(i));
//tipsm(i) = tpm(i)==0 ? 0 : tips_params(tpm(i)-1);
}

for(int i=0; i < tpf.size() ; i++){
tipsf(i) = tpf(i)>2 ? tips_params(tpf(i)-1) : tpf(i)==2 ? 0 : tips_params(tpf(i));
//tipsf(i) = tpf(i)==0 ? 0 : tips_params(tpf(i)-1);
}

vector<Type> mum;
mum = exp(Dm*spline_params_m+tipsm)*Em;
vector<Type> muf;
muf = exp(Df*spline_params_f+tipsf)*Ef;

Type nll;
nll = GMRF(full_penal_m)(spline_params_m)+GMRF(full_penal_f)(spline_params_f)-sum(dpois(dm,mum,1))-sum(dpois(df,muf,1));

return nll;
}
