// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

#include <omp.h>                 // OpenMP 多线程支持
#include <RcppArmadillo.h>       // RcppArmadillo 提供矩阵和线性代数操作
#include <RcppEigen.h>           // RcppEigen 提供矩阵和线性代数操作
#include <iostream>              // C++ 标准输入输出流
#include <cmath>                 // 数学函数库

using namespace Eigen;
using namespace arma;

// [[Rcpp::export]]
VectorXd Cweight(double kh, int k){
  VectorXd w(k+1);
  int d,j;
  for(j=0;j<=k;j++){
    d = j;
    w(j) = (1/kh) * ((2 * kh - d) * ((2 * kh) > d) - (kh - d) * (kh > d));
  }
  return w;  
}

// [[Rcpp::export]]
double CP(int m, int n){
  int s = 1;
  int i;
  for (i=1;i<=m;i++){
    s = s * (n - (i - 1));
  }
  return s;
}

// [[Rcpp::export]]
double CadT(MatrixXd x){
  int M1 = x.rows();
  int M2 = x.cols();
  VectorXd a(M1);
  int i;
  for(i=0;i<M1;i++){
    a(i)=1;
  }
  double mf1 = 1 / CP(2, M1);
  double t1 = 0;
  double s1 = 0;
  for (i=0;i< M1;i++){
    double t2 = 0;
    double s2 = 0;
    VectorXd temp1 = x.row(i);
    int j;
    for (j=0;j< M1;j++){
      VectorXd temp2 = x.row(j);
      t2 = mf1 * (temp1.dot(temp2))*(temp1.dot(temp2));
      s2 = s2 + t2;
    }
    t1 = s2 - mf1 * (temp1.dot(temp1))*(temp1.dot(temp1));
    s1 = s1 + t1;
  }
  return s1;
}

// [[Rcpp::export]]
VectorXd Cbandpen1(MatrixXd x, int i, int j,int q){
  int M1 = x.rows();
  int M2 = x.cols();
  double temp1 = 0;
  double sband1 = 0;
  double temppen1 = 0;
  double penalty1 = 0;
  double mfb1 = 1 /CP(2, M1);
  double C1 = 1 /CP(2, M1);
  int n;
  for (n=0;n<(M2 - q);n++){
    temp1 = mfb1 * x(i, n) * x(j, n+q) * x(j, n) * x(i, n+q);
    temppen1 = C1 * (x(i, n)*x(i, n)) * (x(j, n+q)*x(j, n+q));
    sband1 = sband1 + temp1;
    penalty1 = penalty1 + temppen1;
  }
  VectorXd c(2);
  c(0)=sband1;
  c(1)=penalty1;
  return c;
}


// [[Rcpp::export]]
VectorXd Cbandpen2(MatrixXd x, int q){
  int M1 = x.rows();
  int M2 = x.cols();
  VectorXd temp2=VectorXd::Zero(2);
  VectorXd sbandpen2=VectorXd::Zero(2);
  int i;
  for (i=0;i<M1;i++){
    VectorXd temp3=VectorXd::Zero(2);
    VectorXd sbandpen3=VectorXd::Zero(2);
    int j;
    for (j=0;j<M1;j++){
      temp3 = Cbandpen1(x, i, j, q);
      sbandpen3 = sbandpen3 + temp3;
    }
    temp2 = sbandpen3 - Cbandpen1(x, i, i, q);
    sbandpen2 = sbandpen2 + temp2;
  }
  return sbandpen2;
}

// [[Rcpp::export]]
MatrixXd Cbandpen(MatrixXd x, int K){
  MatrixXd sbandpen4 = MatrixXd::Zero(K + 1, 2);
  sbandpen4.row(0) = Cbandpen2(x, 0);
  int l;
  for (l=0;l<K;l++){
    sbandpen4.row(l + 1)= 2 * Cbandpen2(x, l+1);
  }
  return sbandpen4;
}

// [[Rcpp::export]]
MatrixXd Cbandloss(MatrixXd x, int K){
  int M1 = x.rows();
  int M2 = x.cols();
  MatrixXd tempres = Cbandpen(x, K);
  VectorXd loss=VectorXd::Zero(K+1);
  double tmpAdT=CadT(x);
  
  int i;
  for (i=0;i<(K + 1);i++){
    loss(i) = (1/((double)M2)) * (tmpAdT - tempres.block(0,0,(i+1),1).sum()) + (1/((double)M1)) * (1/((double)M2)) * tempres.block(0,1,(i+1),1).sum();
  }
  
  // MatrixXd c=MatrixXd::Zero(K + 1, 1);
  // c.col(0)=loss;
  return loss;
}

// [[Rcpp::export]]
VectorXd Ccal_banding(MatrixXd X, int K){
  MatrixXd resprop = Cbandloss(X, K);
  VectorXd bandprop = resprop.col(0);
  VectorXd::Index minb;
  double IXband = bandprop.minCoeff(&minb);
  
  VectorXd c=VectorXd::Zero(1);
  c(0)=minb+1;
  return c;
}

// [[Rcpp::export]]
MatrixXd Ctaperloss(MatrixXd x, int K){
  int M1 = x.rows();
  int M2 = x.cols();
  MatrixXd tempres = Cbandpen(x, K);
  VectorXd loss=VectorXd::Zero(K+1);
  double tmpAdT=CadT(x);
  
  int Lt=(K / 2);
  double t0=(1/((double)M2)) * (tmpAdT- tempres(0, 0)) + (1/((double)M1)) * (1/((double)M2)) * tempres(0, 1);
  loss(0)=t0;
  
  int i;
  for (i=0;i<(K + 1);i++){
    if (i <Lt){
      VectorXd w=Cweight(i+1, K);
      VectorXd wgttape=VectorXd::Zero(K+1);
      VectorXd wgtpen=VectorXd::Zero(K+1);
      int j;
      for(j=0;j<K+1;j++)
      {
        wgttape(j)=1-(1-w(j))*(1-w(j));
        wgtpen(j)=w(j)*w(j);
      }
      loss(i+1) = (1/((double)M2)) * (tmpAdT - wgttape.dot(tempres.col(0))) +  (1/((double)M1)) * (1/((double)M2)) * wgtpen.dot(tempres.col(1));
    }
  }
  
  // MatrixXd c=MatrixXd::Zero(K + 1, 1);
  // c.col(0)=loss;
  return loss;
}

// [[Rcpp::export]]
VectorXd Ccal_tapering(MatrixXd X, int K){
  MatrixXd resprop = Ctaperloss(X, K);
  int Lt=(K/2);
  VectorXd tapeprop = resprop.block(0,0,Lt,0).col(0);
  VectorXd::Index mint;
  double IXtape = tapeprop.minCoeff(&mint);
  
  VectorXd c=VectorXd::Zero(1);
  c(0)=mint+1;
  return c;
}