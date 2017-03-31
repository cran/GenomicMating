#include <RcppArmadillo.h>
#include <cmath>


using namespace arma; 
using namespace Rcpp;




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::mat Kmatfunc(arma::mat Markers) {
  arma::vec pks(Markers.n_cols);
  double c=0;
  for (int iter = 0; iter < Markers.n_cols; ++iter){
    pks(iter)=sum(Markers.col(iter)+ones(Markers.n_rows))/(2*Markers.n_rows);
    c=c+2*pks(iter)*(1-pks(iter));
    }
  arma::mat W=Markers+1-2*ones(Markers.n_rows,1)*pks.t();
  arma::mat Amat=(1/c)*W*W.t();
  return Amat;

}


// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec mapfunct(arma::vec x){
  arma::vec output(2);

  if (((x(0)==1) & (x(1)==1))){

    output(0)=2;
    output(1)=0;
  }
  if (((x(0)==1)  &(x(1)==0))){

    output(0)=1.5;
    output(1)=.25;
    
  }
  if (((x(0)==0)  &(x(1)==1))){

    output(0)=1.5;
    output(1)=.25;
  }
  if (((x(0)==1)  &(x(1)==-1))){

    output(0)=1;
    output(1)=0;
  }
  if (((x(0)==-1)  &(x(1)==1))){

    output(0)=1;
    output(1)=0;
  }
  if (((x(0)==0)  &(x(1)==0))){

    output(0)=1;
    output(1)=.5;
  }
  if (((x(0)==0)  &(x(1)==-1))){

    output(0)=.5;
    output(1)=.25;
  }
  if (((x(0)==-1)  & (x(1)==0))){

    output(0)=.5;
    output(1)=.25;
  }
  if (((x(0)==-1)  &(x(1)==-1))){

    output(0)=0;
    output(1)=0;
  }
  return output;
}





// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
double calculatecrossvalue(arma::vec parent1, arma::vec parent2,arma::vec markereffects, double impvar){
  arma::uvec tempbool=find(markereffects<0);
  
  parent1(tempbool)=-parent1(tempbool);
  parent2(tempbool)=-parent2(tempbool);
  arma::vec absmarkereffects=abs(markereffects);
 
  arma::mat p1p2 = join_rows(parent1, parent2);
  p1p2=p1p2.t();
  

  arma::mat scores(2,p1p2.n_cols);
  for (int iter = 0; iter < p1p2.n_cols; ++iter)
  {
    
    scores.col(iter)=mapfunct(p1p2.col(iter));
  }
 
  arma::mat scoresv=scores.row(0)+(impvar)*sqrt(scores.row(1))-1;

  double output=as_scalar(scoresv*absmarkereffects);
  return output;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec getstats(arma::mat Markers, arma::mat K, arma::vec markereffects, arma::mat P, double impvar) {
  arma::vec output(3);
  arma::vec EBV=Markers*markereffects;
  output(0)=as_scalar(ones(1,P.n_rows)*P*EBV);

  arma::uvec tempbool;
  arma::vec inbreedvec=(K.diag()-ones(K.n_cols));
  arma::vec f1f2(2);
  double psi;
  arma::vec Psi(P.n_rows);

  for (int iter = 0; iter < P.n_rows; ++iter)
  {

tempbool=find(P.row(iter)>0);
if (tempbool.n_elem==1){

  f1f2(0)=as_scalar(inbreedvec.elem(tempbool));
  f1f2(1)=as_scalar(inbreedvec.elem(tempbool));
}else {f1f2=inbreedvec.elem(tempbool);}

psi=.5-(1/4)*(f1f2(0)+f1f2(1));

Psi(iter)=psi;
}  

  output(1)=as_scalar(ones(1,P.n_rows)*(P*K*P.t()+diagmat(Psi))*ones(P.n_rows,1));
  arma::vec crossvalues(P.n_rows);
  arma::mat Parents(2,Markers.n_cols);
  for (int iter = 0; iter < P.n_rows; ++iter)
  {
  tempbool=find(P.row(iter)>0);
  if (tempbool.n_elem==1){
    Parents.row(0)=Markers.rows(tempbool);
    Parents.row(1)=Markers.rows(tempbool);
    
  } else {Parents=Markers.rows(tempbool);}
  
  crossvalues(iter)=calculatecrossvalue(Parents.row(0).t(), Parents.row(1).t(),markereffects,impvar);
  }
  output(2)=sum(crossvalues);

  return output;
}




