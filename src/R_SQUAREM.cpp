//TODO:
//1. Exception Handling, try&catch 1.Done
//2. additional input to each of the functions
//3. Test squarem1 on R and compare with the result by original package
//4. If squarem1 passed,then implement the other 3+1 functions
//5. Search for cpp equivalent 	coef <- try(solve(qr(crossprod(U), LAPACK=TRUE, tol=1.e-14), rep(1,K+1)), silent=TRUE)
//6. Search for cpp equivalent 	coef <- try(solve(qr(U[,-(K+1)], LAPACK=TRUE, tol=1.e-14), -U[,K+1]), silent=TRUE)
//7. Allow user to specify the control variables

//Increase # of inputs: testingY
//GITHUB implementation
//Increase # of functions
//SEXP pointer?

#include <iostream>
#include <string>
#include <algorithm>
#include <Rcpp.h>
//#include <RInside.h>
#include <cmath>
#include <math.h>
//#include <stdarg.h>

using namespace Rcpp;

//Global Control Variables
const int K_=1;
const int method_=3;//1,2,3 indicates the types of step length to be used in squarem1,squarem2, 4,5 for "rre" and "mpe" in cyclem1 and cyclem2,  standing for reduced-rank ("rre") or minimal-polynomial ("mpe") extrapolation.
// K=1 must go with method=1,2 or 3
// K>1 must go with method=4 or 5
const int mstep=4;
const int maxiter=1500;
const bool square=true;
const bool trace=true;//currently set to be true for debugging purpose
const double stepmin0=1;
const double stepmax0=1;
const double kr=1;
const double objfninc=1;//0 to enforce monotonicity, Inf for non-monotonic scheme, 1 for monotonicity far from solution and allows for non-monotonicity closer to solution
const double tol=1e-7;

List squarem1(NumericVector par,Function fixptfn,Function objfn,NumericVector testingY);
/*
List squarem2(NumericVector par,Function fixptfn);
List cyclem1(NumericVector par,Function fixptfn,Function objfn=NULL);
List cyclem2(NumericVector par,Function fixptfn);
*/

// [[Rcpp::export]]
Rcpp::List cxxSQUAREM(NumericVector par,NumericVector testingY,Function fixptfn,Function objfn=NULL)
{
    //char fff_strings[1024];
    //va_list argptr;
    //std::cout<<argptr<<std::endl;
    //va_start(argptr,objfn);
    //List fff_test=argptr;
    if(Rf_isNull(objfn)){
        std::cout<<"No objective function supplied"<<std::endl;
    }
    
    List sqobj;
    if (objfn!=NULL) {
        if (K_ == 1){
            sqobj=squarem1(par, fixptfn, objfn,testingY);}
        else if(K_ > 1 || method_>3){
            //sqobj=cyclem1(par, fixptfn, objfn );
        }
        else{sqobj=NULL;}
    } else {
        if (K_ == 1){
            //sqobj = squarem2(par, fixptfn);
        }
        else if (K_ > 1 || method_>3){
            //sqobj = cyclem2(par, fixptfn);
        }
        else{sqobj=NULL;}
    }
    return sqobj;
}

// [[Rcpp::export]]
List squarem1(NumericVector par,Function fixptfn,Function objfn,NumericVector testingY){
    SEXP lold,lnew,p,p1,p2;//R data types
    NumericVector loldcpp,lnewcpp,pcpp,p1cpp,p2cpp,pnew;
    NumericVector q1,q2,sr2,sq2,sv2,srv;
    double sr2_scalar,sq2_scalar,sv2_scalar,srv_scalar,alpha,stepmin,stepmax;
    int iter,feval,leval;
    bool conv,extrap;
    stepmin=stepmin0;
    stepmax=stepmax0;
    //pcpp=p;
    
    if(trace){std::cout<<"Squarem-1"<<std::endl;}
    if(objfn==NULL){
        std::cout<<"Squarem2 should be used if objective function is not available!"<<std::endl;
        return 0;
    }
    iter=1;p=par;
    try{lold=objfn(p,testingY);leval=1;}
    catch(...){
        std::cout<<"Error in fixptfn function evaluation";
        return 1;
    }
    loldcpp=lold;
    lnew=lold;
    if(trace){std::cout<<"Objective fn: "<<loldcpp[0]<<std::endl;}
    feval=0;conv=true;
    
    
    while(feval<maxiter){
        //Step 1
        extrap = true;
        try{p1=fixptfn(p,testingY);feval++;}
        catch(...){
            std::cout<<"Error in fixptfn function evaluation";
            return 1;
        }
        pcpp=p;
        p1cpp=p1;
        q1=p1cpp-pcpp;//conversion from SEXP to vector for crossprod
        sr2=q1*q1;
        sr2_scalar=0;
        for (int i=0;i<sr2.length();i++){sr2_scalar+=sr2[i];}
        if(sqrt(sr2_scalar)<tol){break;}
        
        //Step 2
        try{p2=fixptfn(p1,testingY);feval++;}
        catch(...){
            std::cout<<"Error in fixptfn function evaluation";
            return 1;
        }
        p2cpp=p2;
        q2=p2cpp-p1cpp;
        sq2=q2*q2;
        sq2_scalar=0;
        for (int i=0;i<q2.length();i++){sq2_scalar+=sq2[i];}
        sq2_scalar=sqrt(sq2_scalar);
        if (sq2_scalar<tol){break;}
        sv2=q2-q1;
        sv2_scalar=0;
        for (int i=0;i<sv2.length();i++){sv2_scalar+=sv2[i]*sv2[i];}
        srv_scalar=0;
        for (int i=0;i<q2.length();i++){srv_scalar+=sv2[i]*q1[i];}
        //std::cout<<"sr2,sv2,srv="<<sr2_scalar<<","<<sv2_scalar<<","<<srv_scalar<<std::endl;//debugging

        //Step 3 Proposing new value
        switch (method_){
            case 1: alpha= -srv_scalar/sv2_scalar;
            case 2: alpha= -sr2_scalar/srv_scalar;
            case 3: alpha= sqrt(sr2_scalar/sv2_scalar);
            //default: {
            //    std::cout<<"Misspecification in method, when K=1, method should be either 1, 2 or 3!";
            //    break;}
        }
        
        alpha=std::max(stepmin,std::min(stepmax,alpha));
        //std::cout<<"alpha="<<alpha<<std::endl;//debugging
        pnew = pcpp + 2.0*alpha*q1 + alpha*alpha*(q2-q1);
        
        //Step 4 stabilization
        if(std::abs(alpha-1)>0.01){
            try{pnew=fixptfn(pnew,testingY);feval++;}
            catch(...){
                pnew=p2cpp;
                try{lnew=objfn(pnew,testingY);leval++;}
                catch(...){
                    lnew=lold;
                }
                if(alpha==stepmax){
                    stepmax=std::max(stepmax0,stepmax/mstep);
                }
                alpha=1;
                extrap=false;
                if(alpha==stepmax){stepmax=mstep*stepmax;}
                if(stepmin<0 && alpha==stepmin){stepmin=mstep*stepmin;}
                p=pnew;
                lnewcpp=lnew;
                if(!std::isnan(lnewcpp[0])){lold=lnew;}
                if(trace){std::cout<<"Objective fn: "<<lnewcpp[0]<<"  Extrapolation: "<<extrap<<"  Steplength: "<<alpha<<std::endl;}
                iter++;
                continue;//next round in while loop
            }
            
            if (isfinite(objfninc)){
                try{lnew=objfn(pnew,testingY);leval++;}
                catch(...){
                    pnew=p2;
                    try{lnew=objfn(pnew,testingY);leval++;}
                    catch(...){
                        std::cout<<"Error in objfn function evaluation";
                        return 1;
                    }
                    if(alpha==stepmax){
                        stepmax=std::max(stepmax0,stepmax/mstep);
                    }
                    alpha=1;
                    extrap=false;
                }
            }else{lnew=lold;}
            lnewcpp=lnew;
            if (lnewcpp[0]>loldcpp[0]+objfninc) {
                pnew=p2;
                try{lnew=objfn(pnew,testingY);leval++;}
                catch(...){
                    std::cout<<"Error in objfn function evaluation";
                    return 1;
                }
                if(alpha==stepmax){
                    stepmax=std::max(stepmax0,stepmax/mstep);
                }
                alpha=1;
                extrap=false;
            }
        }else{//same as above, when stablization is not performed.
            if (isfinite(objfninc)){
                try{lnew=objfn(pnew,testingY);leval++;}
                catch(...){
                    pnew=p2;
                    try{lnew=objfn(pnew,testingY);leval++;}
                    catch(...){
                        std::cout<<"Error in objfn function evaluation";
                        return 1;
                    }
                    if(alpha==stepmax){
                        stepmax=std::max(stepmax0,stepmax/mstep);
                    }
                    alpha=1;
                    extrap=false;
                }
            }else{lnew=lold;}
            lnewcpp=lnew;
            if (lnewcpp[0]>loldcpp[0]+objfninc) {
                pnew=p2;
                try{lnew=objfn(pnew,testingY);leval++;}
                catch(...){
                    std::cout<<"Error in objfn function evaluation";
                    return 1;
                }
                if(alpha==stepmax){
                    stepmax=std::max(stepmax0,stepmax/mstep);
                }
                alpha=1;
                extrap=false;
            }
        }
        
        
        if(alpha==stepmax){stepmax=mstep*stepmax;}
        if(stepmin<0 && alpha==stepmin){stepmin=mstep*stepmin;}
        
        p=pnew;
        lnewcpp=lnew;
        if(!std::isnan(lnewcpp[0])){lold=lnew;}
        loldcpp=lold;
        if(trace){std::cout<<"Objective fn: "<<lnewcpp[0]<<"  Extrapolation: "<<extrap<<"  Steplength: "<<alpha<<std::endl;}
        
        iter++;
        //std::cout<<"leval="<<leval<<std::endl;//debugging
    }
    
    if (feval >= maxiter){conv=false;}
    if (!isfinite(objfninc)){lold=objfn(p,testingY);leval++;}
    
    return(List::create(Named("par")=p,
                        Named("value.objfn")=lold,
                        Named("iter")=iter,
                        Named("fpevals")=feval,
                        Named("objfevals")=leval,
                        Named("convergence")=conv));
}





List squarem2(NumericVector par,Function fixptfn){
    List sqobj;
    return sqobj;
}
List cyclem1(NumericVector par,Function fixptfn,Function objfn){
    List sqobj;
    return sqobj;
}
List cyclem2(NumericVector par,Function fixptfn,Function objfn){
    List sqobj;
    return sqobj;
}

List squaremtest(NumericVector par,Function fixptfn,Function objfn){
    List sqobj;
    return sqobj;
}

//main() used for debugging in Xcode
int main(){
    std::cout<<"Hello "<<tol<<std::endl;
    double tol1=tol*10000000;
    std::cout<<"Hello "<<tol1<<std::endl;
    return 0;
}

