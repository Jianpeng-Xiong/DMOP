// ========================================================|//
// The 14 test functions are for cec2018 competition on    |
// dynamic multiobjective optimisation. This document is   |
// free to disseminate for academic use.                   |
// --------------------------------------------------------|//
// The "time" term in the test suite is defined as:        |
//          t=1/nt*floor(tau/taut)                         |
// where - nt:    severity of change                       |
//       - taut:  frequency of change                      |
//       - tau:   current generation counter               |
// --------------------------------------------------------|//
// Any questions can be directed to                        |
//    Dr. Shouyong Jiang at math4neu@gmail.com.            |
//                                                         |


//This code was downloaded directly from the official website by Xiong. 
#ifndef __CEC2018_DF_H_
#define __CEC2018_DF_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <algorithm>    // std::iter_swap

using namespace std;

#define pi 3.141592653589793238462643383279502884197169399375105
#define max(a,b) (a)>(b)? (a):(b)

double helper_sum(vector<double> &x, int start, int end, double y)
// including start index wile excluding end index
{
    double s=0;
    for (int i=start; i<end; i++)
    {
        s+=pow(x[i]-y,2);
    }
    return s;
}

vector<double> cec2018_DF(char* probID, vector<double> x, int tau, int taut, int nt){
    // INPUT:
    //       probID: test problem identifier (i.e. 'DF1')
    //       x:      variable vector
    //       tau:    current generation counter
    //       taut:   frequency of change
    //       nt:     severity of change
    //
    // RETURN:       objective vector
    //
    
    
    // the first change occurs after T0 generations, that is, the
    //generation at which a change occurs is (T0+1), (T0+taut+1), etc.
    int T0=50;
    
    // calculate time instant
    int tau_tmp=max(tau+taut-(T0+1),0);
    double t=1/nt*floor(tau_tmp/taut);
    
    int n=x.size();
    
    if (!strcmp(probID, "DF1"))
    {
        double G=fabs(sin(0.5*pi*t));
        double H=0.75*sin(0.5*pi*t)+1.25;
        double g=1+helper_sum(x, 1, n, G);
        
        vector<double> f(2,0);
        f[0]=x[0];
        f[1]=g*(1-pow(x[0]/g, H));
        return f;
    }
    else if (!strcmp(probID, "DF2"))
    {
        double G=fabs(sin(0.5*pi*t));
        int r=floor((n-1)*G);
        vector<double> y(x);
        std::iter_swap(y.begin(), y.begin()+r);
        double g=1+helper_sum(y, 1, n, G);
        
        vector<double> f(2,0);
        f[0]=x[0];
        f[1]=g*(1-pow(x[0]/g, 0.5));
        return f;
    }
    else if (!strcmp(probID, "DF3"))
    {
        double G=sin(0.5*pi*t);
        double H=G+1.5;
        double g=1+helper_sum(x, 1, n, G+pow(x[0],H));
        
        vector<double> f(2,0);
        f[0]=x[0];
        f[1]=g*(1-pow(x[0]/g, H));
        return f;
    }
    else if (!strcmp(probID, "DF4"))
    {
        double a=sin(0.5*pi*t);
        double b=1+fabs(cos(0.5*pi*t));
        double c=max(fabs(a),a+b);
        double H=1.5+a;
        double g=1;
        for (int i=1; i<n; i++)
            g+=pow(x[i]-a*pow(x[0]/c,2)/(i+1),2);
        
        vector<double> f(2,0);
        f[0]=g*pow(fabs(x[0]-a),H);
        f[1]=g*pow(fabs(x[0]-a-b),H);
        return f;
    }
    else if (!strcmp(probID, "DF5"))
    {
        double G=sin(0.5*pi*t);
        int w=floor(10*G);
        double g=1+helper_sum(x, 1, n, G);
        
        vector<double> f(2,0);
        f[0]=g*(x[0]+0.02*sin(w*pi*x[0]));
        f[1]=g*(1-x[0]+0.02*sin(w*pi*x[0]));
        return f;
    }
    else if (!strcmp(probID, "DF6"))
    {
        double G=sin(0.5*pi*t);
        double a=0.2+2.8*fabs(G);
        double g=1;
        for (int i=1; i<n; i++){
            double y=x[i]-G;
            g+=fabs(G)*pow(y,2)-10*cos(2*pi*y)+10;
        }
    
        vector<double> f(2,0);
        f[0]=g*pow(x[0]+0.1*sin(3*pi*x[0]),a);
        f[1]=g*pow(1-x[0]+0.1*sin(3*pi*x[0]),a);
        return f;
    }
    else if (!strcmp(probID, "DF7"))
    {
        double a=5*cos(0.5*pi*t);
        double tmp=1.0/(1+exp(a*(x[0]-2.5)));
        double g=1+helper_sum(x, 1, n, tmp);
        
        vector<double> f(2,0);
        f[0]=g*(1+t)/x[0];
        f[1]=g*x[0]/(1+t);
        return f;
    }
    else if (!strcmp(probID, "DF8"))
    {
        double G=sin(0.5*pi*t);
        double a=2.25+2*cos(2*pi*t);
        double b=100*pow(G,2);
        double tmp=G*sin(4*pi*pow(x[0],b))/(1+fabs(G));
        double g=1+helper_sum(x, 1, n, tmp);
        
        vector<double> f(2,0);
        f[0]=g*(x[0]+0.1*sin(3*pi*x[0]));
        f[1]=g*pow(1-x[0]+0.1*sin(3*pi*x[0]),a);
        return f;
    }
    else if (!strcmp(probID, "DF9"))
    {
        double N=1+floor(10*fabs(sin(0.5*pi*t)));
        double g=1;
        for (int i=1; i<n; i++){
            double tmp=x[i]-cos(4*t+x[0]+x[i-1]);
            g+=pow(tmp,2);
        }
        
        vector<double> f(2,0);
        f[0]=g*(1-x[0]+max(0, (0.1+0.5/N)*sin(2*N*pi*x[0])));
        f[1]=g*(1-x[0]+max(0, (0.1+0.5/N)*sin(2*N*pi*x[0])));
        return f;
    }
    else if (!strcmp(probID, "DF10"))
    {
        double G=sin(0.5*pi*t);
        double H=2.25+2*cos(0.5*pi*t);
        double tmp=sin(2*pi*(x[0]+x[1]))/(1+fabs(G));
        double g=1+helper_sum(x, 2, n, tmp);
        
        vector<double> f(3,0);
        f[0]=g*pow(sin(0.5*pi*x[0]),H);
        f[1]=g*pow(sin(0.5*pi*x[1])*cos(0.5*pi*x[0]),H);
        f[2]=g*pow(cos(0.5*pi*x[1])*cos(0.5*pi*x[0]),H);
        return f;
    }
    else if (!strcmp(probID, "DF11"))
    {
        double G=fabs(sin(0.5*pi*t));
        double g=1+G+helper_sum(x, 2, n, 0.5*G*x[0]);
        
        vector<double> f(3,0);
        double y1=pi*G/6+(pi/2-pi*G/3)*x[0];
        double y2=pi*G/6+(pi/2-pi*G/3)*x[1];
        
        f[0]=g*sin(y1);
        f[1]=g*sin(y2)*cos(y1);
        f[2]=g*cos(y2)*cos(y1);
        return f;
    }
    else if (!strcmp(probID, "DF12"))
    {
        double k=10*sin(pi*t);
        double g=1+helper_sum(x, 2, n, sin(t*x[0]));
        g+=fabs(sin(floor(k*(2*x[0]-1))*pi/2)*sin(floor(k*(2*x[1]-1))*pi/2));
        
        vector<double> f(3,0);
        f[0]=g*cos(0.5*pi*x[1])*cos(0.5*pi*x[0]);
        f[1]=g*sin(0.5*pi*x[1])*cos(0.5*pi*x[0]);
        return f;
        f[2]=g*sin(0.5*pi*x[0]);
    }
    else if (!strcmp(probID, "DF13"))
    {
        double G=sin(0.5*pi*t);
        int p=floor(6*G);
        double g=1+helper_sum(x, 2, n, G);
        
        vector<double> f(3,0);
        f[0]=g*cos(0.5*pi*x[0]);
        f[1]=g*cos(0.5*pi*x[1]);
        f[2]=g*pow(sin(0.5*pi*x[0]),2)+sin(0.5*pi*x[0])*pow(cos(p*pi*x[0]),2) +
        pow(sin(0.5*pi*x[1]),2)+sin(0.5*pi*x[1])*pow(cos(p*pi*x[1]),2);
        return f;
    }
    else if (!strcmp(probID, "DF14"))
    {
        double G=sin(0.5*pi*t);
        double g=1+helper_sum(x, 2, n, G);
        
        vector<double> f(3,0);
        double y=0.5+G*(x[1]-0.5);
        f[0]=g*(1-y+0.05*sin(6*pi*y));
        f[1]=g*(1-x[1]+0.05*sin(6*pi*x[1]))*(y+0.05*sin(6*pi*y));
        f[2]=g*(x[1]+0.05*sin(6*pi*x[1]))*(y+0.05*sin(6*pi*y));
        return f;
    }
    
    cout<<"Error: no test problem matched."<<endl;
    
    vector<double> f(2,0);
    return f;
}
#endif
