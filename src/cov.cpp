#include <RcppEigen.h>
#include <map>          // to map kernels to integers for the switch
#include <string>       // to read in the kernel name
#include <vector>       // to use vectors
#include <algorithm>    // to get the intersect, sort, lower_bound, upper_bound

typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

using Rcpp::List;
using Rcpp::as;

const double SQRT2 = sqrt(2);
const double invSqrt2pi=  1./(sqrt(2.*M_PI));

// [[Rcpp::export]]
Rcpp::List crawcov(const Rcpp::List & Lt,
                   const Rcpp::List & Ly,// assume that mean has been subtracted from y
                   const Eigen::VectorXd weig)
{
    int n = Lt.size();

    // construct raw cov and arrange it according to clocpoly
    int m = 0;
    for(int i = 0; i < n; i++)
    {
        const MapVecd& Ti = as<MapVecd>(Lt[i]);
        int Ni = Ti.size();
        if(Ni > 1)
        {
            m = m + Ni*(Ni-1);
        }
    }


    Eigen::MatrixXd tpairs(m,2);
    Eigen::VectorXd cxx(m);
    Eigen::VectorXd w(m);

    int k = 0;
    for(int i = 0; i < n; i++)
    {
        const MapVecd& Ti = as<MapVecd>(Lt[i]);
        int Ni = Ti.size();
        if(Ni > 1)
        {
            const MapVecd& Yi = as<MapVecd>(Ly[i]);

            for(int u = 0; u < Ni; u++)
            {
                for(int v = 0; v < Ni; v++)
                {

                    if(u != v)
                    {
                        tpairs(k,0) = Ti(u);
                        tpairs(k,1) = Ti(v);
                        cxx(k) = Yi(u) * Yi(v);
                        w(k) = weig(i);
                        k++;
                    }
                }
            }
        }
    }


    List outList;
    outList["cxx"] = cxx;
    outList["tpairs"] = tpairs;
    outList["w"] = w;

    return Rcpp::wrap(outList);
}


/* The following is adapted from  RmullwlskUniversal.cpp in fdapace package,
 * developed by Hans-Georg Mueller and Jane-Ling Wang (ORGANIZATION: University of California, Davis)

 License: BSD_3_clause + file LICENSE

 Copyright (c) <YEAR>, <COPYRIGHT HOLDER>

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in
 the documentation and/or other materials provided with the
 distribution.

 Neither the name of the <ORGANIZATION> nor the names of its
 contributors may be used to endorse or promote products derived
 from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */



typedef std::pair<double, unsigned int> valIndPair;
bool comparePair(const valIndPair& l, const valIndPair& r) {
    return(l.first < r.first);
}


inline Eigen::VectorXd get_w(const int kcode, const Eigen::MatrixXd& llx, const Eigen::VectorXd& w)
{
    Eigen::VectorXd temp(w.size());

    switch (kcode){
    case 1: // Epan
        temp=  ((1-llx.row(0).array().pow(2))*(1- llx.row(1).array().pow(2))).array() *
            ((9./16)*w).transpose().array();
        break;
    case 2 : // Rect
        temp=(w.array())*.25 ;
        break;
    case 3:  // triangular
        temp = (1-llx.row(0).array().abs()) * (1-llx.row(1).array().abs()) * (w.transpose().array());
        break;
    case 4:  // quartic
        temp = (w.transpose().array()) *
            ((1.-llx.row(0).array().pow(2)).array().pow(2)).array() *
            ((1.-llx.row(1).array().pow(2)).array().pow(2)).array() * (225./256.);
        break;
    case 5: // triweight
        temp = (w.transpose().array()) *
            ((1.-llx.row(0).array().pow(2)).array().pow(3)).array() *
            ((1.-llx.row(1).array().pow(2)).array().pow(3)).array() * (1225./1024.);
        break;
    case 6: // tricube
        temp = (w.transpose().array()) *
            ((1.-llx.row(0).array().pow(3)).array().pow(3)).array() *
            ((1.-llx.row(1).array().pow(3)).array().pow(3)).array() * (4900./6561.);
        break;
    case 7: // cosine
        temp = (w.transpose().array()) *
            (((M_PI/2)*llx.row(0).array()).cos()).array() *
            (((M_PI/2)*llx.row(1).array()).cos()).array() * (M_PI*M_PI/4);
        break;
    case 21: // gauss
        temp = ((-.5*(llx.row(0).array().pow(2))).exp()) * invSqrt2pi  *
            ((-.5*(llx.row(1).array().pow(2))).exp()) * invSqrt2pi  *
            (w.transpose().array());
        break;
    case 22: // logistic
        temp = w.transpose().array() / ((2+llx.row(0).array().exp()+(-llx.row(0).array()).exp()) *
            (2+llx.row(1).array().exp()+(-llx.row(1).array()).exp()));
        break;
    case 23: // sigmoid
        temp = ((4./(M_PI*M_PI)) * w.transpose().array()) / ((llx.row(0).array().exp()+(-llx.row(0).array()).exp()) *
            (llx.row(1).array().exp()+(-llx.row(1).array()).exp()));
        break;
    case 24:  // silverman
        //return exp(-abs(z)/SQRT2) * sin(abs(z)/SQRT2+M_PI/4);
        temp = (0.5 * w.transpose().array()) *
            ((llx.row(0).array().abs() / (-SQRT2)).exp() * ((M_PI/4)+(llx.row(0).array().abs() / SQRT2)).sin()) *
            ((llx.row(1).array().abs() / (-SQRT2)).exp() * ((M_PI/4)+(llx.row(1).array().abs() / SQRT2)).sin());
        break;
    default: Rcpp::stop("Unsupported kernel");

    }

    return temp;
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd csmoothcov( const Eigen::Map<Eigen::VectorXd> & bw,
                            const std::string kernel_type,
                            const Eigen::Map<Eigen::MatrixXd> & tPairs, // 2 by m matrix
                            const Eigen::Map<Eigen::MatrixXd> & cxxn,
                            const Eigen::Map<Eigen::VectorXd> & win,
                            const Eigen::Map<Eigen::VectorXd> & xgrid,
                            const Eigen::Map<Eigen::VectorXd> & ygrid,
                            const bool & autoCov,
                            const double delta){
    // Assumes the first row of tPairs is sorted in increasing order.
    // tPairs : xin (in MATLAB code)
    // cxxn : yin (in MATLAB code)
    // xgrid: out1 (in MATLAB code)
    // ygrid: out2 (in MATLAB code)
    // autoCov : boolean / cause the function to return the autocovariance.



    // map kernels, code <= 20 is reserved for compactly supported kernel on [-1,1]
    std::map<std::string,int> supported_kernels;
    supported_kernels["epanechnikov"]    = 1;
    supported_kernels["rectangular"]    = 2;
    supported_kernels["triangular"]   = 3;
    supported_kernels["quartic"]    = 4;
    supported_kernels["triweight"]    = 5;
    supported_kernels["tricube"]   = 6;
    supported_kernels["cosine"] = 7;
    supported_kernels["gauss"]   = 21;
    supported_kernels["logistic"] = 22;
    supported_kernels["sigmoid"] = 23;
    supported_kernels["silverman"] = 24;

    // The following test is here for completeness, we mightwant to move it up a
    // level (in the wrapper) in the future.

    // If the kernel_type key exists set kcode appropriately
    int kcode = 0;
    if ( supported_kernels.count( kernel_type ) != 0){
        kcode = supported_kernels.find( kernel_type )->second; //Set kernel choice
    } else {
        // otherwise use "epan"as the kernel_type
        Rcpp::warning("Kernel_type argument was not set correctly; Epanechnikov kernel used.");
        kcode = supported_kernels.find( "epanechnikov" )->second;;
    }

    // Check that we do not have zero weights // Should do a try-catch here
    // Again this might be best moved a level-up.
    if ( !(win.all()) ){  //
        Rcpp::Rcout << "Cases with zero-valued windows are not yet implemented" << std::endl;
        return (tPairs);
    }

    // ProfilerStart("sort.log");
    // Start the actual smoother here
    const unsigned int xgridN = xgrid.size();
    const unsigned int ygridN = ygrid.size();
    const unsigned int n = tPairs.cols();

    // For sorted x1
    Eigen::VectorXd x1(tPairs.row(0).transpose());
    const double* tDat = x1.data();

    Eigen::MatrixXd mu(xgrid.size(), ygrid.size());
    mu.setZero();


    for (unsigned int i = 0; i != xgridN; ++i) {
        const double xl = xgrid(i) - bw(0) - 1e-6,
            xu = xgrid(i) + bw(0) + 1e-6;

        unsigned int indl = std::lower_bound(tDat, tDat + n, xl) - tDat,
            indu = std::upper_bound(tDat, tDat + n, xu) - tDat;

        // sort the y index
        std::vector<valIndPair> yval(indu - indl);
        for (unsigned int k = 0; k < yval.size(); ++k){
            yval[k] = std::make_pair(tPairs(1, k + indl), k + indl);
        }
        std::sort<std::vector<valIndPair>::iterator>(yval.begin(), yval.end(), comparePair);

        std::vector<valIndPair>::iterator ylIt = yval.begin(),
            yuIt = yval.begin();

        for (unsigned int j = !autoCov? 0: i; j != ygridN; ++j) {
            double del = xgrid(i) - ygrid(j);
            if (del < 0) del = -del;
            if(del > delta)
            {
                mu(i,j) = R_PosInf;
                continue;
            }

            const double yl = ygrid(j) - bw(1) - 1e-6,
                yu = ygrid(j) + bw(1) + 1e-6;


            //locating local window (LOL) (bad joke)
            std::vector <unsigned int> indx;

            //if the kernel is compactly supported
            if ( kcode <= 20) {
                // Search the lower and upper bounds increasingly.
                ylIt = std::lower_bound(ylIt, yval.end(), valIndPair(yl, 0), comparePair);
                yuIt = std::upper_bound(yuIt, yval.end(), valIndPair(yu, 0), comparePair);

                // The following works nice for the Gaussian
                //  but for very small samples it complains
                //} else {
                //  ylIt = yval.begin();
                //  yuIt = yval.end();
                //}

                for (std::vector<valIndPair>::iterator y = ylIt; y != yuIt; ++y){
                    indx.push_back(y->second);
                }
            }
            else { //When we finally get c++11 we will use std::iota
                for( unsigned int y = 0; y != n; ++y){
                    indx.push_back(y);
                }
            }

            unsigned int indxSize = indx.size();


            Eigen::VectorXd lw(indxSize);
            Eigen::VectorXd ly(indxSize);
            Eigen::MatrixXd lx(2,indxSize);

            for (unsigned int u = 0; u !=indxSize; ++u){
                lx.col(u) = tPairs.col(indx[u]);
                lw(u) = win(indx[u]);
                ly(u) = cxxn(indx[u]);
            }

            // check enough points are in the local window
            unsigned int meter=1;
            for (unsigned int u =0; u< indxSize; ++u) {
                for (unsigned int t = u + 1; t < indxSize; ++t) {
                    if ( (lx(0,u) !=  lx(0,t) ) || (lx(1,u) != lx(1,t) ) ) {
                        meter++;
                    }
                }
                if (meter >= 3) {
                    break;
                }
            }

            //compute weight matrix
            if (meter >=  3) {
                Eigen::VectorXd temp(indxSize);
                Eigen::MatrixXd llx(2, indxSize );
                llx.row(0) = (lx.row(0).array() - xgrid(i))/bw(0);
                llx.row(1) = (lx.row(1).array() - ygrid(j))/bw(1);

                temp = get_w(kcode,llx,lw);

                // make the design matrix
                Eigen::MatrixXd X(indxSize ,3);
                X.setOnes();
                X.col(1) = lx.row(0).array() - xgrid(i);
                X.col(2) = lx.row(1).array() - ygrid(j);
                Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWX(X.transpose() * temp.asDiagonal() *X);
                // The solver should stop if the value is NaN. See the HOLE example in gcvlwls2dV2.
                Eigen::VectorXd beta = ldlt_XTWX.solve(X.transpose() * temp.asDiagonal() * ly);
                mu(i,j)=beta(0);
            }
            else if(meter < 3){
                Rcpp::stop("No enough points in local window, please increase bandwidth using userBwCov.");

            }
        }
    }

    if(autoCov){
        return ( Eigen::MatrixXd(mu.triangularView<Eigen::StrictlyUpper>().transpose()) + Eigen::MatrixXd(mu.triangularView<Eigen::Upper>()));
    } else {
        // ProfilerStop();
        return ( mu );
    }
}

