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
arma::vec mapfunctM1(arma::vec x){
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
double calculatecrossvalueM1(arma::vec parent1, arma::vec parent2,arma::vec markereffects){
  arma::uvec tempbool=find(markereffects<0);
  
//parent1(tempbool)=parent1(tempbool);
 // parent2(tempbool)=parent2(tempbool);
  
  arma::vec absmarkereffects=abs(markereffects);
 
  arma::mat p1p2 = join_rows(parent1, parent2);
  p1p2=p1p2.t();
  

  arma::mat scores(2,p1p2.n_cols);
  for (int iter = 0; iter < p1p2.n_cols; ++iter)
  {
    scores.col(iter)=mapfunctM1(p1p2.col(iter));
  }
 
  arma::mat scoresv=scores.row(1);

  double output=as_scalar(scoresv*pow(2*absmarkereffects,2));
  return output;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec getstatsM1(arma::mat Markers, arma::mat K, arma::vec markereffects, arma::mat P) {
  arma::vec output(3);
  arma::vec EBV=Markers*markereffects;
  output(1)=as_scalar(ones(1,P.n_rows)*P*EBV);

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

  output(0)=as_scalar(ones(1,P.n_rows)*(P*K*P.t()+diagmat(Psi))*ones(P.n_rows,1));
  arma::vec crossvalues(P.n_rows);
  arma::mat Parents(2,Markers.n_cols);
  for (int iter = 0; iter < P.n_rows; ++iter)
  {
  tempbool=find(P.row(iter)>0);
  if (tempbool.n_elem==1){
    Parents.row(0)=Markers.rows(tempbool);
    Parents.row(1)=Markers.rows(tempbool);
    
  } else {Parents=Markers.rows(tempbool);}
  
  crossvalues(iter)=calculatecrossvalueM1(Parents.row(0).t(), Parents.row(1).t(),markereffects);
  }
  output(2)=sum(crossvalues);

  return output;
}


//////////////////////////////////////////////////////////


// [[Rcpp::depends("RcppArmadillo")]]
double computeDsM2(arma::vec L1, arma::vec L2) {
    //# D is zero
    if ((L1(0)==L1(1))||(L2(0)==L2(1))) {
        return 0;
    }
    //# D is non-zero
    //# Coupling
    else if ((L1(0)==L2(0))||(L1(1)==L2(1))) {
        return 0.25;
    }
    //# Repulsion
    else if ((L1(0)==L2(1))||(L1(0)==L2(1))) {
        return -.25;
    }
    //# Parent has heterotsygote
    else {
        return 0;
    }
}



// [[Rcpp::depends("RcppArmadillo")]]
double computeCM2(double L1, double L2, double M1, double M2) {
    //# Marker in the different chromosome -> 0.5
    
    double valuec;
    if (M1!=M2) {
        valuec = 0.5;
    }
    else if ((M1==M2)&&(std::abs(L1-L2) > 0.5)) {
        valuec = 0.5;
    }
    else {
        valuec = std::abs(L2 - L1);
    }
    valuec = 0.5*(1-std::pow(M_E,-2*valuec));
    
    return valuec;
}


// [[Rcpp::depends("RcppArmadillo")]]
double computeCrossVarM2(double Ds, double C, unsigned char type, unsigned int k) {
    // type: type of population, from Table 1 in Lehermeier et al. 2017. Values:
    //       0: DH F1 (k=1, k parameter not used)
    //       1: DH generation k (actually if k==1, the result is the same as type==0)
    //       2: RILs generation k
    //       3 (or any other value): DH and RILS generation Inf (k parameter not used)
    // k: "generation k when DH lines are generated (k = 1 for DH from F1)" From Lehermeier et al. 2017
    double cvar;
    if (type==0) {
        cvar = 4 * Ds * (1-2*C);
    } else if (type == 1) {
        cvar = 0;
        for (int r=1; r<=k; r++) {
            cvar += std::pow((1-2*C)/2, r);
        }
        cvar = 4*Ds*(cvar + std::pow((1-2*C)/2, k));
    } else if (type == 2) {
        cvar = 0;
        for (int r=1; r<=k; r++) {
            cvar += std::pow((1-2*C)/2, r);
        }
        cvar *= 4*Ds;
    } else {
        cvar = 4*Ds*(1 - 4*C / (1-2*C));
    }
    
    return cvar;
}


// [[Rcpp::depends("RcppArmadillo")]]
double calculatecrossvalueM2(arma::vec Parents1, arma::vec Parents2,arma::vec markereffects, arma::mat markermap, unsigned char type=0, unsigned int generation=1){
    //# Take columns of genotype according to cross parents
    
    // type: type of population, from Table 1 in Lehermeier et al. 2017. Values:
    //       0: DH F1 (k=1, k parameter not used)
    //       1: DH generation k (actually if k==1, the result is the same as type==0)
    //       2: RILs generation k
    //       3 (or any other value): DH and RILS generation Inf (k parameter not used)
    // k: "generation k when DH lines are generated (k = 1 for DH from F1)" From Lehermeier et al. 2017
    
    int nmarkers = Parents1.size();
    arma::mat genotype2(2, nmarkers);
    genotype2.row(0) =Parents1.t();
    genotype2.row(1) =Parents2.t();
    //# Create X matrix
    arma::mat X=zeros(nmarkers,nmarkers);
    
    //# Computing diagonal
    for (int i = 0; i < nmarkers-1; i++) {
        
        //# Same allele = 0
        if(genotype2(0,i) == genotype2(1,i)){
            X(i,i) = 0;
            
        } else {
            X(i,i)  = 1;
        }
    }
    double Ds;
    double C;
    double D;
    //# computing off-diagonal
    for (int j = 0; j < nmarkers-1; j++) {
        
        for (int l = j+1; l < nmarkers; l++) {
            //#do computation
            
            Ds =computeDsM2(genotype2.col(j), genotype2.col(l));
            C = computeCM2(markermap(j,1), markermap(l,1), markermap(j,0), markermap(l,0));
            D = computeCrossVarM2(Ds, C, type, generation);
            X(j,l) = D;
            X(l,j) = D;
        }
    }
    
    //Rcout << "Armadillo matrix is" << std::endl << X << std::endl;
    //# Variance for a cross
    
    double  crossvariance = as_scalar(markereffects.t()*X*markereffects);
    
    return crossvariance;
    
}




// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec getstatsM2(arma::mat Markers, arma::mat K, arma::vec markereffects, arma::mat P, arma::mat markermap,  unsigned char type, unsigned int generation) {
  double crossvalout;
  arma::vec output(3);
  arma::vec EBV=Markers*markereffects;
  output(1)=as_scalar(ones(1,P.n_rows)*P*EBV);
  
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
  
  output(0)=as_scalar(ones(1,P.n_rows)*(P*K*P.t()+diagmat(Psi))*ones(P.n_rows,1));
  arma::vec crossvalues(P.n_rows);
  arma::mat Parents(2,Markers.n_cols);
  for (int iter = 0; iter < P.n_rows; ++iter)
  {
    tempbool=find(P.row(iter)>0);
    if (tempbool.n_elem==1){
      Parents.row(0)=Markers.rows(tempbool);
      Parents.row(1)=Markers.rows(tempbool);
      
    } else {Parents=Markers.rows(tempbool);}
    
    crossvalout=calculatecrossvalueM2(Parents.row(0).t(), Parents.row(1).t(),markereffects, markermap, type, generation);
    crossvalues(iter)=crossvalout;
  }
  output(2)=sum(crossvalues);
  
  return output;
}




/////////////////////////////////////////////

// poliploidy is a parameter. divide varx with polip-1
// [[Rcpp::depends("RcppArmadillo")]]
arma::vec mapfunctM3(arma::vec x){
  arma::vec output(2);
  double ex=(1*x(0)*x(1)+.5*(x(0)*(1-x(1))+x(1)*(1-x(0))));
  
  double varx=(powf(1-ex,2)*x(0)*x(1)+powf(.5-ex,2)*(x(0)*(1-x(1))+x(1)*(1-x(0))))+powf(0-ex,2)*(1-x(1))*(1-x(0));
  output(0)=ex;
  output(1)=varx;

  return output;
}





// [[Rcpp::depends("RcppArmadillo")]]
double calculatecrossvalueM3(arma::vec parent1, arma::vec parent2,arma::vec markereffects){
  
  
  arma::mat p1p2 = join_rows(parent1, parent2);
  p1p2=p1p2.t();
  
  
  arma::mat scores(2,p1p2.n_cols);
  for (int iter = 0; iter < p1p2.n_cols; ++iter)
  {
    scores.col(iter)=mapfunctM3(p1p2.col(iter));
  }
  
  arma::mat scoresv=scores.row(1);
  
  //Rcout << "The value 1 is " << scores.is_colvec()<< std::endl;
  // Rcout << "The value 2 is " << markereffects.is_colvec()<< std::endl;
  
  // check if this is actually a crossproduct.
  double output=as_scalar(scoresv*pow(markereffects,2));
  
  return output;
}



// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]
arma::vec getstatsM3(arma::mat Markers, arma::mat K, arma::vec markereffects, arma::mat P) {
  arma::vec output(3);
  arma::vec EBV=Markers*markereffects;
  output(1)=as_scalar(ones(1,P.n_rows)*P*EBV);
  
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
  
  output(0)=as_scalar(ones(1,P.n_rows)*(P*K*P.t()+diagmat(Psi))*ones(P.n_rows,1));
  arma::vec crossvalues(P.n_rows);
  arma::mat Parents(2,Markers.n_cols);
  for (int iter = 0; iter < P.n_rows; ++iter)
  {
    tempbool=find(P.row(iter)>0);
    if (tempbool.n_elem==1){
      Parents.row(0)=Markers.rows(tempbool);
      Parents.row(1)=Markers.rows(tempbool);
      
    } else {Parents=Markers.rows(tempbool);}
    
    crossvalues(iter)=calculatecrossvalueM3(Parents.row(0).t(), Parents.row(1).t(),markereffects);
  }
  output(2)=sum(crossvalues);
 
  return output;
}






