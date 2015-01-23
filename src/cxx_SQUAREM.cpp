//TODO:
//1. example of fixptfn and objfn need to be added
//2. create class object
//3. cyclem equivalent
//4. combine squarem1&2 by defining the default objfn

#include <iostream>
#include <string>
#include <algorithm>
//#include <Rcpp.h>//A pure c++ SQUAREM
#include <cmath>
#include <math.h>
#include <vector>
#include <string>
#include <cstring>
//#include <stdarg.h>

using namespace std;//Rcpp;

struct SquaremControl{
    int K=1;
    int method=3;//1,2,3 indicates the types of step length to be used in squarem1,squarem2, 4,5 for "rre" and "mpe" in cyclem1 and cyclem2,  standing for reduced-rank ("rre") or minimal-polynomial ("mpe") extrapolation.
    // K=1 must go with method=1,2 or 3
    // K>1 must go with method=4 or 5
    double mstep=4;
    int maxiter=1500;
    bool square=true;
    bool trace=true;//currently set to be true for debugging purpose
    double stepmin0=1;
    double stepmax0=1;
    double kr=1;
    double objfninc=1;//0 to enforce monotonicity, Inf for non-monotonic scheme, 1 for monotonicity far from solution and allows for non-monotonicity closer to solution
    double tol=1e-7;
} SquaremDefault;


struct SquaremOutput{
    std::vector<double> par;
    double valueobjfn;
    int iter;
    int pfevals;
    int objfevals;
    bool convergence;
} sqobj,sqobjnull;

vector<double> fixptfn(std::vector<double> par,std::vector<double> testingY){
    std::vector<double> parnew=par;
    return parnew;
};
double objfn(std::vector<double> par,std::vector<double> testingY){
    double objvalue=testingY[0];
    return objvalue;
};

SquaremOutput squarem1(std::vector<double> par,std::vector<double> testingY);
/*
List squarem2(std::vector<double> par,Function fixptfn);
List cyclem1(std::vector<double> par,Function fixptfn,Function objfn=NULL);
List cyclem2(std::vector<double> par,Function fixptfn);
*/


// [[Rcpp::export]]
SquaremOutput cxxSQUAREM(std::vector<double> par,std::vector<double> testingY)
{
    if(SquaremDefault.K == 1){
        sqobj=squarem1(par,testingY);
    }
    else{
        sqobj=sqobjnull;
    }
    //int check_result=memcmp(sqobj,sqobjnull,sizeof(sqobj));
    //if(check_result==0){std::cout<<"Error in fixed-point/objective function evaluation"<<std::endl;}
    return sqobj;
}


// [[Rcpp::export]]
SquaremOutput squarem1(std::vector<double> par,std::vector<double> testingY){
    //std::vector<double> p,p1,p2;//R data types
    double loldcpp,lnewcpp;
    std::vector<double> pcpp,p1cpp,p2cpp,pnew;
    std::vector<double> q1,q2,sr2,sq2,sv2,srv;
    double sr2_scalar,sq2_scalar,sv2_scalar,srv_scalar,alpha,stepmin,stepmax;
    int iter,feval,leval;
    bool conv,extrap;
    stepmin=SquaremDefault.stepmin0;
    stepmax=SquaremDefault.stepmax0;
    if(SquaremDefault.trace){std::cout<<"Squarem-1"<<std::endl;}
    
    iter=1;pcpp=par;
    try{loldcpp=objfn(pcpp,testingY);leval=1;}
    catch(...){
        std::cout<<"Error in fixptfn function evaluation";
        return sqobjnull;
    }
    lnewcpp=loldcpp;
    if(SquaremDefault.trace){std::cout<<"Objective fn: "<<loldcpp<<std::endl;}
    feval=0;conv=true;
    
    long int parvectorlength=pcpp.size();
    
    while(feval<SquaremDefault.maxiter){
        //Step 1
        extrap = true;
        try{p1cpp=fixptfn(pcpp,testingY);feval++;}
        catch(...){
            std::cout<<"Error in fixptfn function evaluation";
            return sqobjnull;
        }
        sr2_scalar=0;
        for (int i=0;i<parvectorlength;i++){sr2_scalar+=pow(p1cpp[i]-pcpp[i],2);}
        if(sqrt(sr2_scalar)<SquaremDefault.tol){break;}

        
        //Step 2
        try{p2cpp=fixptfn(p1cpp,testingY);feval++;}
        catch(...){
            std::cout<<"Error in fixptfn function evaluation";
            return sqobjnull;
        }
        sq2_scalar=0;
        for (int i=0;i<parvectorlength;i++){sq2_scalar+=pow(p2cpp[i]-p1cpp[i],2);}
        sq2_scalar=sqrt(sq2_scalar);
        if (sq2_scalar<SquaremDefault.tol){break;}
        sv2_scalar=0;
        for (int i=0;i<parvectorlength;i++){sv2_scalar+=pow(p2cpp[i]-2*p1cpp[i]+pcpp[i],2);}
        srv_scalar=0;
        for (int i=0;i<parvectorlength;i++){srv_scalar+=(p2cpp[i]-2*p1cpp[i]+pcpp[i])*(p1cpp[i]-pcpp[i]);}
        //std::cout<<"sr2,sv2,srv="<<sr2_scalar<<","<<sv2_scalar<<","<<srv_scalar<<std::endl;//debugging

        //Step 3 Proposing new value
        switch (SquaremDefault.method){
            case 1: alpha= -srv_scalar/sv2_scalar;
            case 2: alpha= -sr2_scalar/srv_scalar;
            case 3: alpha= sqrt(sr2_scalar/sv2_scalar);
        }
        
        alpha=std::max(stepmin,std::min(stepmax,alpha));
        //std::cout<<"alpha="<<alpha<<std::endl;//debugging
        for (int i=0;i<parvectorlength;i++){pnew[i]=pcpp[i]+2.0*alpha*(p1cpp[i]-pcpp[i])+pow(alpha,2)*(p2cpp[i]-2*p1cpp[i]+pcpp[i]);}
        //pnew = pcpp + 2.0*alpha*q1 + alpha*alpha*(q2-q1);
        
        //Step 4 stabilization
        if(std::abs(alpha-1)>0.01){
            try{pnew=fixptfn(pnew,testingY);feval++;}
            catch(...){
                pnew=p2cpp;
                try{lnewcpp=objfn(pnew,testingY);leval++;}
                catch(...){
                    lnewcpp=loldcpp;
                }
                if(alpha==stepmax){
                    stepmax=std::max(SquaremDefault.stepmax0,stepmax/SquaremDefault.mstep);
                }
                alpha=1;
                extrap=false;
                if(alpha==stepmax){stepmax=SquaremDefault.mstep*stepmax;}
                if(stepmin<0 && alpha==stepmin){stepmin=SquaremDefault.mstep*stepmin;}
                pcpp=pnew;
                if(!std::isnan(lnewcpp)){loldcpp=lnewcpp;}
                if(SquaremDefault.trace){std::cout<<"Objective fn: "<<lnewcpp<<"  Extrapolation: "<<extrap<<"  Steplength: "<<alpha<<std::endl;}
                iter++;
                continue;//next round in while loop
            }
            
            if (isfinite(SquaremDefault.objfninc)){
                try{lnewcpp=objfn(pnew,testingY);leval++;}
                catch(...){
                    pnew=p2cpp;
                    try{lnewcpp=objfn(pnew,testingY);leval++;}
                    catch(...){
                        std::cout<<"Error in objfn function evaluation";
                        return sqobjnull;
                    }
                    if(alpha==stepmax){
                        stepmax=std::max(SquaremDefault.stepmax0,stepmax/SquaremDefault.mstep);
                    }
                    alpha=1;
                    extrap=false;
                }
            }else{lnewcpp=loldcpp;}
            if (lnewcpp>loldcpp+SquaremDefault.objfninc) {
                pnew=p2cpp;
                try{lnewcpp=objfn(pnew,testingY);leval++;}
                catch(...){
                    std::cout<<"Error in objfn function evaluation";
                    return sqobjnull;
                }
                if(alpha==stepmax){
                    stepmax=std::max(SquaremDefault.stepmax0,stepmax/SquaremDefault.mstep);
                }
                alpha=1;
                extrap=false;
            }
        }else{//same as above, when stablization is not performed.
            if (isfinite(SquaremDefault.objfninc)){
                try{lnewcpp=objfn(pnew,testingY);leval++;}
                catch(...){
                    pnew=p2cpp;
                    try{lnewcpp=objfn(pnew,testingY);leval++;}
                    catch(...){
                        std::cout<<"Error in objfn function evaluation";
                        return sqobjnull;
                    }
                    if(alpha==stepmax){
                        stepmax=std::max(SquaremDefault.stepmax0,stepmax/SquaremDefault.mstep);
                    }
                    alpha=1;
                    extrap=false;
                }
            }else{lnewcpp=loldcpp;}
            if (lnewcpp>loldcpp+SquaremDefault.objfninc) {
                pnew=p2cpp;
                try{lnewcpp=objfn(pnew,testingY);leval++;}
                catch(...){
                    std::cout<<"Error in objfn function evaluation";
                    return sqobjnull;
                }
                if(alpha==stepmax){
                    stepmax=std::max(SquaremDefault.stepmax0,stepmax/SquaremDefault.mstep);
                }
                alpha=1;
                extrap=false;
            }
        }
        
        
        if(alpha==stepmax){stepmax=SquaremDefault.mstep*stepmax;}
        if(stepmin<0 && alpha==stepmin){stepmin=SquaremDefault.mstep*stepmin;}
        
        pcpp=pnew;
        if(!std::isnan(lnewcpp)){loldcpp=lnewcpp;}
        if(SquaremDefault.trace){std::cout<<"Objective fn: "<<lnewcpp<<"  Extrapolation: "<<extrap<<"  Steplength: "<<alpha<<std::endl;}
        
        iter++;
        //std::cout<<"leval="<<leval<<std::endl;//debugging
    }
    if (feval >= SquaremDefault.maxiter){conv=false;}
    if (!isfinite(SquaremDefault.objfninc)){loldcpp=objfn(pcpp,testingY);leval++;}
    
    //assigning values
    sqobj.par=pcpp;
    sqobj.valueobjfn=loldcpp;
    sqobj.iter=iter;
    sqobj.pfevals=feval;
    sqobj.objfevals=leval;
    sqobj.convergence=conv;
    return(sqobj);
}


/*TODO
List squarem2(std::vector<double> par,Function fixptfn){
    List sqobj;
    return sqobj;
}
List cyclem1(std::vector<double> par,Function fixptfn,Function objfn){
    List sqobj;
    return sqobj;
}
List cyclem2(std::vector<double> par,Function fixptfn,Function objfn){
    List sqobj;
    return sqobj;
}

List squaremtest(std::vector<double> par,Function fixptfn,Function objfn){
    List sqobj;
    return sqobj;
}
*/
//main() used for debugging in Xcode
int main(){
    std::cout<<"Hello "<<SquaremDefault.tol<<std::endl;
    return 0;
}

