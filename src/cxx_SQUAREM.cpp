//TODO:
//1. example of fixptfn and objfn need to be added 1. Done
//2. create class object
//3. cyclem equivalent
//4. combine squarem1&2 by defining the default objfn
//5. matrix/vector operation without loop
#include <iostream>
#include <algorithm>
#include <cmath>
#include <math.h>
#include <vector>
#include <numeric>


using namespace std;

//Global Control Variable
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

//Output Struct
struct SquaremOutput{
    std::vector<double> par;
    double valueobjfn;
    int iter=0;
    int pfevals=0;
    int objfevals=0;
    bool convergence=false;
} sqobj,sqobjnull;

vector<double> fixptfn(std::vector<double> par);
double objfn(std::vector<double> par);
SquaremOutput squarem1(std::vector<double> par);


//Wrapper function
SquaremOutput cxxSQUAREM(std::vector<double> par)
{
    if(SquaremDefault.K == 1){
        sqobj=squarem1(par);
    }
    else{
        sqobj=sqobjnull;
    }
    return sqobj;
}

//actual acceleration function
SquaremOutput squarem1(std::vector<double> par){
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
    
    iter=1;pcpp=par;pnew=par;
    try{loldcpp=objfn(pcpp);leval=1;}
    catch(...){
        std::cout<<"Error in fixptfn function evaluation";
        return sqobjnull;
    }
    lnewcpp=loldcpp;
    if(SquaremDefault.trace){std::cout<<"Objective fn: "<<loldcpp<<std::endl;}
    feval=0;conv=true;
    
    const long int parvectorlength=pcpp.size();
    
    while(feval<SquaremDefault.maxiter){
        //Step 1
        extrap = true;
        try{p1cpp=fixptfn(pcpp);feval++;}
        catch(...){
            std::cout<<"Error in fixptfn function evaluation";
            return sqobjnull;
        }
        
        sr2_scalar=0;
        for (int i=0;i<parvectorlength;i++){sr2_scalar+=pow(p1cpp[i]-pcpp[i],2);}
        if(sqrt(sr2_scalar)<SquaremDefault.tol){break;}

        //Step 2
        try{p2cpp=fixptfn(p1cpp);feval++;}
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
            try{pnew=fixptfn(pnew);feval++;}
            catch(...){
                pnew=p2cpp;
                try{lnewcpp=objfn(pnew);leval++;}
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
                try{lnewcpp=objfn(pnew);leval++;}
                catch(...){
                    pnew=p2cpp;
                    try{lnewcpp=objfn(pnew);leval++;}
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
                try{lnewcpp=objfn(pnew);leval++;}
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
                try{lnewcpp=objfn(pnew);leval++;}
                catch(...){
                    pnew=p2cpp;
                    try{lnewcpp=objfn(pnew);leval++;}
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
                try{lnewcpp=objfn(pnew);leval++;}
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
    if (!isfinite(SquaremDefault.objfninc)){loldcpp=objfn(pcpp);leval++;}
    
    //assigning values
    sqobj.par=pcpp;
    sqobj.valueobjfn=loldcpp;
    sqobj.iter=iter;
    sqobj.pfevals=feval;
    sqobj.objfevals=leval;
    sqobj.convergence=conv;
    return(sqobj);
}



//main() used for demostration

std::vector<int> Other_input1 {162,267,271,185,111,61,27,8,3,1};

int main(){
    std::cout<<"Hi, this is a demostration using Poisson mixture!"<<std::endl;
    std::vector<double> par_initial {0.5,3,1};//some random starting point,
    SquaremOutput SQ_result=cxxSQUAREM(par_initial);
    cout<<"Initial input vector:"<<endl;
    for(int i=0;i<par_initial.size();i++) cout<<par_initial[i]<<" ";
    cout<<"\nOutput vector by SQUAREM:"<<endl;
    for(int i=0;i<par_initial.size();i++) cout<<SQ_result.par[i]<<" ";
    cout<<"\nNumber of evaluations on fixed point function:"<<endl;
    cout<<SQ_result.pfevals<<endl;
    cout<<"Number of evaluations on objective function:"<<endl;
    cout<<SQ_result.objfevals<<endl;
    return 0;
}


//Fixed point function and objective function that take only the vector for EM acceleration as input, while others are defined as global variables(const within the scope of cxxSQUAREM)

vector<double> fixptfn(std::vector<double> par){
    std::vector<double> parnew=par;
    
    std::vector<double> zi(Other_input1.size());
    for (int i=0;i<Other_input1.size();i++){
        zi[i]=par[0]*exp(-par[1])*pow(par[1],i)/(par[0]*exp(-par[1])*pow(par[1],i)+(1-par[0])*exp(-par[2])*pow(par[2],i));
    }
    
    
    double temp1,temp2;
    //parnew[0]=1;
    temp1=0,temp2=0;
    for (int i=0;i<Other_input1.size();i++){
        temp1+=Other_input1[i]*zi[i];
        temp2+=Other_input1[i];
    }
    parnew[0]=temp1/temp2;
    
    
    temp1=0,temp2=0;
    for (int i=0;i<Other_input1.size();i++){
        temp1+=Other_input1[i]*zi[i]*i;
        temp2+=Other_input1[i]*zi[i];
    }
    parnew[1]=temp1/temp2;
    
    temp1=0,temp2=0;
    for (int i=0;i<Other_input1.size();i++){
        temp1+=Other_input1[i]*(1-zi[i])*i;
        temp2+=Other_input1[i]*(1-zi[i]);
    }
    parnew[2]=temp1/temp2;
    
    return parnew;
}


double objfn(std::vector<double> par){
    double objvalue=0;
    vector<double> loglik(Other_input1.size());
    for (int i=0;i<Other_input1.size();i++){
    loglik[i]=Other_input1[i]*log(par[0]*exp(-par[1])*pow(par[1],i)/exp(lgamma(i+1))+
        (1-par[0])*exp(-par[2])*pow(par[2],i)/exp(lgamma(i+1)));
        objvalue-=loglik[i];
    }
    return objvalue;
};
