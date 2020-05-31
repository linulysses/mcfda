#include <RcppEigen.h>
#include <map>      
#include <string>       
#include <vector>       
#include <algorithm>    

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


typedef std::pair<double, unsigned int> Vip;
bool compare(const Vip& p1, const Vip& p2) {
    return(p1.first < p2.first);
}


inline Eigen::VectorXd get_w(const int kcode, const Eigen::MatrixXd& xy, const Eigen::VectorXd& w)
{
    Eigen::VectorXd temp(w.size());
    
    switch (kcode){
    case 1: // Epan
        temp=  ((1-xy.col(0).array().pow(2))*(1- xy.col(1).array().pow(2))).array() *
            ((9./16)*w).array();
        break;
    case 2 : // Rect
        temp=(w.array())*.25 ;
        break;
    case 3:  // triangular
        temp = (1-xy.col(0).array().abs()) * (1-xy.col(1).array().abs()) * (w.array());
        break;
    case 4:  // quartic
        temp = (w.array()) *
            ((1.-xy.col(0).array().pow(2)).array().pow(2)).array() *
            ((1.-xy.col(1).array().pow(2)).array().pow(2)).array() * (225./256.);
        break;
    case 5: // triweight
        temp = (w.array()) *
            ((1.-xy.col(0).array().pow(2)).array().pow(3)).array() *
            ((1.-xy.col(1).array().pow(2)).array().pow(3)).array() * (1225./1024.);
        break;
    case 6: // tricube
        temp = (w.transpose().array()) *
            ((1.-xy.col(0).array().pow(3)).array().pow(3)).array() *
            ((1.-xy.col(1).array().pow(3)).array().pow(3)).array() * (4900./6561.);
        break;
    case 7: // cosine
        temp = (w.array()) *
            (((M_PI/2)*xy.col(0).array()).cos()).array() *
            (((M_PI/2)*xy.col(1).array()).cos()).array() * (M_PI*M_PI/16);
        break;
    case 21: // gauss
        temp = ((-.5*(xy.col(0).array().pow(2))).exp()) * (invSqrt2pi*2) *
            ((-.5*(xy.col(1).array().pow(2))).exp()) * (w.array());
        break;
    case 22: // logistic
        temp = w.array() / ((2+xy.col(0).array().exp()+(-xy.col(0).array()).exp()) *
            (2+xy.col(1).array().exp()+(-xy.col(1).array()).exp()));
        break;
    case 23: // sigmoid
        temp = ((4./(M_PI*M_PI)) * w.array()) / ((xy.col(0).array().exp()+(-xy.col(0).array()).exp()) *
            (xy.col(1).array().exp()+(-xy.col(1).array()).exp()));
        break;
    case 24:  // silverman
        temp = (0.25 * w.array()) *
            ((xy.col(0).array().abs() / (-SQRT2)).exp() * ((M_PI/4)+(xy.col(0).array().abs() / SQRT2)).sin()) *
            ((xy.col(1).array().abs() / (-SQRT2)).exp() * ((M_PI/4)+(xy.col(1).array().abs() / SQRT2)).sin());
        break;
    default: Rcpp::stop("Unsupported kernel");
    
    }
    
    return temp;
}

inline double locpoly(const Eigen::MatrixXd& xy, 
                      const Eigen::VectorXd& z, 
                      const Eigen::VectorXd& w, 
                      const Eigen::VectorXd& h, 
                      const double xi, 
                      const double yj,
                      const int kcode)
{
    
    Eigen::MatrixXd xy_tmp(z.size(),2);
    xy_tmp.col(0) = (xy.col(0).array() - xi) / h[0];
    xy_tmp.col(1) = (xy.col(1).array() - yj) / h[1];
    
    Eigen::VectorXd temp = get_w(kcode,xy_tmp,w);
    
    Eigen::MatrixXd X(z.size() ,3);
    X.setOnes();
    X.col(1) = xy.col(0).array() - xi;
    X.col(2) = xy.col(1).array() - yj;
    Eigen::LDLT<Eigen::MatrixXd> ldlt_XTWX(X.transpose() * temp.asDiagonal() *X);
    Eigen::VectorXd b = ldlt_XTWX.solve(X.transpose() * temp.asDiagonal() * z);
    return b(0);
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]

Eigen::MatrixXd csmoothcov( const Eigen::Map<Eigen::VectorXd> & h,
                            const std::string kernel,
                            const Eigen::Map<Eigen::MatrixXd> & xy, // m by 2 matrix, first column is sorted in increasing order
                            const Eigen::Map<Eigen::MatrixXd> & z, // C(x,y)
                            const Eigen::Map<Eigen::VectorXd> & w, // weight for each observation
                            Eigen::Map<Eigen::VectorXd> & xgrid,
                            Eigen::Map<Eigen::VectorXd> & ygrid,
                            const double delta // the proportion of width of snippet 
)
{
    // kcode <= 20 is reserved for compactly supported kernel on [-1,1]
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
    
    int kcode = 1;
    if ( supported_kernels.count( kernel ) != 0){
        kcode = supported_kernels.find( kernel )->second;
    }
    else
    {
        Rcpp::stop("The specified kernel is not supported");
    }
    
    
    const unsigned int xn = xgrid.size();
    unsigned int yn = ygrid.size();
    const unsigned int n = xy.rows();
    
    //Rcpp::Rcout << "xn=" << xn << ", yn=" << yn << ", n=" << n << ", kcode=" << kcode << std::endl;
    
    bool symmetric = true;
    if(xn!=yn) symmetric = false;
    else
    {
        for(unsigned int k=0; k < xn; k++)
        {
            if(xgrid(k) != ygrid(k))
            {
                symmetric = false;
                break;
            }
        }
    }
    
    //Rcpp::Rcout << "symmetric" << symmetric << std::endl;
    
    Eigen::VectorXd x(xy.col(0));
    const double* xDat = x.data();
    
    Eigen::MatrixXd chat(xn,yn);// = MatrixXd::Zero(xn,yn);
    chat.setOnes();
    
    
    // for each value in the first coordinate
    for (unsigned int i = 0; i < xn; ++i) {
        
        if(kcode > 20) // noncompact kernel
        {
            if(symmetric)
            {
                for( unsigned int j = 0; j <=i; ++j){
                    chat(i,j) = locpoly(xy,z,w,h,xgrid(i),ygrid(j),kcode);
                }
            }
            else
            {
                for( unsigned int j = 0; j < yn; ++j){
                    chat(i,j) = locpoly(xy,z,w,h,xgrid(i),ygrid(j),kcode);
                }
            }
        }
        else // compactly supported kernel
        {
            const double xl = xgrid(i) - h(0);
            const double xu = xgrid(i) + h(0);
            
            unsigned int idxl = std::lower_bound(xDat, xDat + n, xl) - xDat;
            unsigned int idxu = std::upper_bound(xDat, xDat + n, xu) - xDat;
            
            std::vector<Vip> yval(idxu - idxl);
            for (unsigned int k = 0; k < yval.size(); ++k){
                yval[k] = std::make_pair(xy(k + idxl,1), k + idxl);
            }
            std::sort<std::vector<Vip>::iterator>(yval.begin(), yval.end(), compare);
            
            std::vector<Vip>::iterator ylIt = yval.begin();
            std::vector<Vip>::iterator yuIt = yval.begin();
            
            unsigned int jend  = yn-1;
            if(symmetric) jend = i;
            
            for (unsigned int j = 0; j <= jend; ++j) {
                double del = xgrid(i) - ygrid(j);
                if (del < 0) del = -del;
                if(del > delta)
                {
                    chat(i,j) = R_PosInf;
                    continue;
                }
                
                const double yl = ygrid(j) - h(1);
                const double yu = ygrid(j) + h(1);
                
                
                std::vector <unsigned int> idx;
                
                ylIt = std::lower_bound(ylIt, yval.end(), Vip(yl, 0), compare);
                yuIt = std::upper_bound(yuIt, yval.end(), Vip(yu, 0), compare);
                
                for (std::vector<Vip>::iterator y = ylIt; y != yuIt; ++y){
                    idx.push_back(y->second);
                }
                
                if(idx.size() < 5)
                {
                    Rcpp::stop("Less than 5 points in the local window.");
                }
                
                
                // prepare data for 2D local polynomial regression
                unsigned int m = idx.size();
                
                //Rcpp::Rcout << "m="  << m << std::endl;
                
                Eigen::VectorXd w_tmp(m);
                Eigen::VectorXd z_tmp(m);
                Eigen::MatrixXd xy_tmp(m,2);
                
                for (unsigned int k = 0; k < m; ++k){
                    xy_tmp.row(k) = xy.row(idx[k]);
                    w_tmp(k) = w(idx[k]);
                    z_tmp(k) = z(idx[k]);
                }
                
                //Rcpp::Rcout << "h="  << h(0) << std::endl;
                
                chat(i,j) = locpoly(xy_tmp,z_tmp,w_tmp,h,xgrid(i),ygrid(j),kcode);
                
            }
        }
    }
    
    if(symmetric)
        return Eigen::MatrixXd(chat.triangularView<Eigen::StrictlyLower>().transpose()) + 
            Eigen::MatrixXd(chat.triangularView<Eigen::Lower>());
    else 
        return chat;
}