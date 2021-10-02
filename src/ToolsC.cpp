#include <RcppEigen.h>
#include <Rcpp.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <random>
#include <ctime>
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace Eigen;
using namespace std;


// [[Rcpp::export]]
Eigen::MatrixXd SigmaXC(double rho, Eigen::MatrixXd SigmaX){
  
  int p = SigmaX.rows();
  
  for(int i=0; i<p; i++){
    for(int j=0; j<p; j++){
      SigmaX(i,j) = pow(rho, abs(i-j));
    }
  }
  
  return SigmaX;   // covariance matrix for X
}



//   One dimensional kernel function  ####
// Input: 
//   x: m*1 matrix
//   h: scalar
// Output:
//   kern: m*1 matrix

// [[Rcpp::export]]
Eigen::MatrixXd knC(Eigen::MatrixXd x, double h){
  
  Eigen::MatrixXd tmp = x/h;
  int m = x.rows();
  
  Eigen::MatrixXd kern = x;
  for(int i=0; i<m; i++){
    kern(i,0) = 0.75 * (1 - pow(tmp(i,0), 2)) / h * (tmp(i,0)<=1) * (tmp(i,0)>=-1);
  }
  
  return kern;
}




//  Locallinear for beta0   ####
// Input:
//   p, x0, hx: scalar
//   x: n*1 matrix
//   y: n*p matrix
// Output:
//   f: p*1 matrix
//

// [[Rcpp::export]]
Eigen::MatrixXd locallinear0C(double p, double x0, double hx, Eigen::MatrixXd x, Eigen::MatrixXd y){
  
  double eps = pow(10, -8);
  Eigen::MatrixXd u = x.array()-x0;
  Eigen::MatrixXd w = knC(u, hx);  // n*1
  double a = w.sum();
  double b = (w.array() * u.array()).sum();
  double d = (w.array() * u.array() * u.array()).sum();
  double deno = a*d - b*b;
  
  Eigen::MatrixXd f = y.row(0).transpose();
  for(int j=0; j<p; j++){
    double t1 = (w.array() * y.col(j).array()).sum();
    double t2 = (w.array() * y.col(j).array() * u.array()).sum();
    double num1 = d*t1 - b*t2;
    f(j,0) = num1 / (deno+eps);
  }
  
  return f;
  
}



//  Locallinear for beta0 and beta1   ####
// Input:
//   p, x0, hx: scalar
//   x: n*1 matrix
//   y: n*p matrix
// Output:
//   f1: p*1 matrix
//   f2: p*1 matrix
//
// [[Rcpp::export]]
Eigen::MatrixXd locallinear1C(double p, double x0, double hy, Eigen::MatrixXd x, Eigen::MatrixXd y){
  
  double eps = pow(10, -8);
  Eigen::MatrixXd u = x.array()-x0;
  Eigen::MatrixXd w = knC(u, hy);  // n*1
  double a = w.sum();
  double b = (w.array() * u.array()).sum();
  double d = (w.array() * u.array() * u.array()).sum();
  double deno = a*d - b*b;
  
  // Eigen::MatrixXd f1 = y.row(0).transpose();
  Eigen::MatrixXd f2 = y.row(0).transpose();
  for(int j=0; j<p; j++){
    double t1 = (w.array() * y.col(j).array()).sum();
    double t2 = (w.array() * y.col(j).array() * u.array()).sum();
    // double num1 = d*t1 - b*t2;
    // f1(j,0) = num1 / (deno+eps);
    double num2 = -b*t1 + a*t2;
    f2(j,0) = num2 / (deno+eps);
  }
  
  return f2;
  
}



//   Derivative of local polynomial   ####
// Input:
//   hy: scalar
//   xb: n*1 matrix
//   y: n*1 matrix
// Output:
//   f: n*1 matrix
//
// [[Rcpp::export]]
Eigen::MatrixXd dlogetaC(Eigen::MatrixXd xb, Eigen::MatrixXd y, double hy){
  
  int n = y.rows();
  
  Eigen::MatrixXd f = xb;
  for(int i=0; i<n; i++){
    double x0 = xb(i,0);
    Eigen::MatrixXd f2 = locallinear1C(1, x0, hy, xb, y);
    f.row(i) = f2;
  }
  
  return f;
}



//   Score function   ####
// Input:
//   hx, hy: scalar
//   beta: p*1 matrix
//   x: n*p matrix
//   y: n*1 matrix
// Output:
//   yout: n*p matrix
//
// [[Rcpp::export]]
Eigen::MatrixXd cseffC(Eigen::MatrixXd x, Eigen::MatrixXd y, Eigen::MatrixXd beta, double hx, double hy){
  
  int n = x.rows();
  int p = x.cols();
  Eigen::MatrixXd xb = x * beta;
  Eigen::MatrixXd mx0 = x;
  Eigen::MatrixXd yout = x;
  
  Eigen::MatrixXd gest = y;
  for(int i=0; i<n; i++){
    Eigen::MatrixXd xb0 = x.row(i) * beta;
    gest.row(i) = locallinear0C(1, xb0(0,0), hy, xb, y);
  }
  Eigen::MatrixXd vareps = y - gest;
  
  for(int i=0; i<n; i++){
    double x0 = xb(i,0);
    Eigen::MatrixXd mx = locallinear0C(p, x0, hx, xb, x);
    mx0.row(i) = x.row(i) - mx.transpose();
  }
  
  Eigen::MatrixXd fb = dlogetaC(xb, y, hy);  // derivative of g
  
  for(int i=0; i<n; i++){
    yout.row(i) = vareps(i,0) * mx0.row(i) * fb(i,0); // Eqn.(6)
  }
  
  return yout;
}



//   Weighted sum of score function   ####
// Input:
//   x: n*p matrix
//   ally: n*m matrix
//   s, m, hx, hy: scalar
//   tm: 1*m vector
//   beta: p*1 matrix
// Output:
//   yout: n*p matrix
//
// [[Rcpp::export]]
Eigen::MatrixXd allcseffC(Eigen::MatrixXd x, Eigen::MatrixXd ally, double m, Eigen::VectorXd tm, double s, 
                   Eigen::MatrixXd beta, double hx, double hy, double h){
  
  Eigen::MatrixXd delta = ally.row(1).transpose();
  for(int ss=0; ss<m; ss++){
    delta(ss,0) = tm(ss) - tm(s-1);
  }
  Eigen::MatrixXd kern = knC(delta, h);
  
  Eigen::MatrixXd yout = 0*x;
  for(int ss=0; ss<m; ss++){
    if(kern(ss,0)>0){
      Eigen::MatrixXd y = ally.col(ss);
      yout = yout + kern(ss,0)*cseffC(x, y, beta, hx, hy);
    }
  }
  
  return yout;
}




//   Weighted sum of score function   ####
// Input:
//   x: n*p matrix
//   ally: n*m matrix
//   s, m, hx, hy, h: scalar
//   tm: 1*m vector
//   beta: p*1 matrix
// Output:
//   fval: 1*p matrix
//
// [[Rcpp::export]]
Eigen::MatrixXd sumseffC(Eigen::MatrixXd x, Eigen::MatrixXd ally, double m, Eigen::VectorXd tm, double s,
                  Eigen::MatrixXd beta, double hx, double hy, double h){
  
  Eigen::MatrixXd yout = allcseffC(x, ally, m, tm, s, beta, hx, hy, h);
  int n = yout.rows();
  Eigen::MatrixXd fval = yout.colwise().sum() / n;  // Eqn.(11)
  
  return fval;
  
}



// Estimation of g
// Input:
//    h1: scalar
//    x:  n*p matrix
//    betaest:   p*m matrix
//    vec_ally:   nm*1 matrix
//    vec_allxb:  nm*1 matrix
//    gest:  n*m matrix
// Output:
//    gest: n*m matrix
//
// [[Rcpp::export]]
Eigen::MatrixXd gestC(double h1, Eigen::MatrixXd x, Eigen::MatrixXd betaest, Eigen::MatrixXd vec_allxb, 
               Eigen::MatrixXd vec_ally, Eigen::MatrixXd gest){
  
  int n = x.rows();
  int m = betaest.cols();
  
  for(int i=0; i<n; i++){
    for(int s=0; s<m; s++){
      Eigen::MatrixXd xb0 = x.row(i) * betaest.col(s);
      Eigen::MatrixXd a = locallinear0C(1, xb0(0,0), h1, vec_allxb, vec_ally);
      gest(i,s) = a(0,0);
    }
  }
  
  return gest;
}



//  Delete the ith row of X #####
// [[Rcpp::export]]
Eigen::MatrixXd delC(Eigen::MatrixXd x, Eigen::MatrixXd xdel, int n){
  
  int nx = x.rows();
  
  if(n==1){
    xdel = x.bottomRows(nx-1);
  }else{
    if(n==nx){
      xdel = x.topRows(nx-1);
    }else{
      xdel.topRows(n-1) = x.topRows(n-1);
      xdel.bottomRows(nx-n) = x.bottomRows(nx-n);
    }
  } 
  
  return xdel;
}



//   CV for h1   ####
// Input:
//   vh: nh vector
//   x: n*p matrix
//   ally: n*m matrix
//   betaest: p*m matrix
//   vec_allxb: nm*1 matrix
//   vec_y: nm*1 matrix
//   allxb, ally: n*m matrix
// Output:
//   CV: nh * 1 matrix
//
// [[Rcpp::export]]
Eigen::MatrixXd cvC(Eigen::VectorXd vh, Eigen::MatrixXd x, Eigen::MatrixXd ally, Eigen::MatrixXd betaest, Eigen::MatrixXd vec_allxb, 
             Eigen::MatrixXd vec_y, Eigen::MatrixXd vec_cv, Eigen::MatrixXd gest, Eigen::MatrixXd CV){
  
  int nh = CV.rows();
  int n = x.rows();
  int m = betaest.cols();
  for(int ii=0; ii<nh; ii++){
    
    double h= vh(ii);
    int ncv = 0;
    
    for(int i=0; i<n; i++){
      for(int s=0; s<m; s++){
        ncv = ncv + 1;
        Eigen::MatrixXd xb_cv = delC(vec_allxb, vec_cv, ncv);
        Eigen::MatrixXd y_cv = delC(vec_y, vec_cv, ncv);
        Eigen::MatrixXd xb0 = x.row(i) * betaest.col(s);
        Eigen::MatrixXd a = locallinear0C(1, xb0(0,0), h, xb_cv, y_cv);
        gest(i,s) = a(0,0);
      }
    }
    
    CV(ii,0) = (ally - gest).array().pow(2).mean();
  }
  
  return CV;
  
}


// [[Rcpp::export]]
Eigen::VectorXd Q( Eigen::MatrixXd M, Eigen::MatrixXd VX, Eigen::MatrixXd VY, Eigen::MatrixXd W) {
  return M *  (VX.transpose() * W * VX).inverse() * ( VX.transpose() * W * VY  );
}

//M %*% solve(t(vec.X) %*% weight %*% vec.X) %*% (t(vec.X) %*% weight %*% vec.Y)



// [[Rcpp::export]]
Eigen::MatrixXd kroneckerC( Eigen::MatrixXd M1, Eigen::MatrixXd M2) {
  
  return  kroneckerProduct(M1, M2);
  
}



// [[Rcpp::export]]
double TnC( Eigen::MatrixXd U, Eigen::MatrixXd angle, double n, double M) {
  
  return ( U.transpose() * angle * U ).trace() / (n * n * M);
  
}
// Tn[l] <- sum( diag(t(U) %*% angle %*% U) ) / (n * n * M) 


// [[Rcpp::export]]
double NMCTnC( Eigen::MatrixXd U, Eigen::MatrixXd En, Eigen::MatrixXd angle, Eigen::MatrixXd Cn, double n, double M) {
  
  Eigen::MatrixXd UEn = U.cwiseProduct(En);
  double term1 = ( UEn.transpose() * angle * UEn ).trace() / (n * n *M);
  double term2 = 2 * (UEn.transpose() * angle * Cn * UEn).trace() / (n * n * n * M);
  double term3 = (UEn.transpose() * Cn * angle * Cn * UEn).trace() / (n * n * n * n* M);
  
  return  term1 - term2 + term3;
}


// [[Rcpp::export]]
double TrueNMCTnC( Eigen::MatrixXd U, Eigen::MatrixXd En, Eigen::MatrixXd angle, double n, double M) {
  
  Eigen::MatrixXd UEn = U.cwiseProduct(En);
  double term1 = ( UEn.transpose() * angle * UEn ).trace() / (n * n *M);
  
  return  term1;
}


// [[Rcpp::export]]
double BootTnC( Eigen::MatrixXd NewU, Eigen::MatrixXd angle, double n, double M ){
  
  return ( NewU.transpose() * angle * NewU ).trace() / (n * n * M);
  
}
// 
// [[Rcpp::export]]
Eigen::MatrixXd CnC( Eigen::MatrixXd X, double n){
  Eigen::MatrixXd OmegaInv = (X.transpose() * X / n).inverse();
  return X * OmegaInv * X.transpose();
}


// [[Rcpp::export]]
double PatTnC(int n, int M, double h, Eigen::MatrixXd U0,
              Eigen::MatrixXd DisKF, Eigen::MatrixXd delta){
  
  Eigen::MatrixXd Q(M, 1);
  Eigen::MatrixXd V(M, M);
  for(int sm1=0; sm1 < M; sm1++){
    Eigen::MatrixXd U0ij = U0.col(sm1) * (U0.col(sm1).transpose());
    Q(sm1,0) = ( U0ij.cwiseProduct(DisKF) ).sum();
    for(int sm2=0; sm2 < M; sm2++){
      Eigen::MatrixXd Uij = (U0.col(sm1) .cwiseProduct(U0.col(sm2)) ) * ((U0.col(sm1) .cwiseProduct(U0.col(sm2))).transpose());
      V(sm1,sm2) = 2 * ( Uij.cwiseProduct(DisKF.cwiseProduct(DisKF)) ).sum();
    }
  }
  
  double Tn = n*sqrt(h)* (Q.cwiseProduct(delta).sum()
  ) / (delta * delta.transpose()).cwiseProduct(V).sum();
  
  return Tn;
}


