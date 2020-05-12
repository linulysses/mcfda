#include <RcppEigen.h>

using Eigen::VectorXd;
using Rcpp::as;
using Rcpp::List;
typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
typedef Eigen::Map<Eigen::VectorXd> MapVecd;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


Eigen::MatrixXd expM(const Eigen::MatrixXd & X, const Eigen::MatrixXd & S)
{
    Eigen::MatrixXd R(X);
    R += S;
    R.diagonal().array() = 0;

    Eigen::VectorXd D = X.diagonal();
    Eigen::VectorXd Z = S.diagonal().array() / D.array();
    Z = Z.array().exp();
    D = D.array() * Z.array();
    R.diagonal() = D;

    return(R);

    //    (X-diag(diag(X))) + (S-diag(diag(S))) + diag(diag(X)*exp(diag(S)/diag(X)))

}

// [[Rcpp::export]]
Eigen::MatrixXd RgradQ( const Rcpp::List & Lt,
                       const Rcpp::List & Lr,
                       const Eigen::MatrixXd& L,
                       const Rcpp::List & B,
                       const Eigen::MatrixXd& U,
                       const Eigen::MatrixXd& V,
                       const Eigen::MatrixXd& W,
                       const double lam,
                       const Eigen::VectorXd& weight)
{
    int n = Lt.size();
    int p = L.rows();
    const Eigen::MatrixXd & C = L * L.transpose();
    Eigen::MatrixXd H = Eigen::MatrixXd::Constant(p, p, 0);

    for(int i=0; i<n; i++)
    {
        const MapVecd& Ti = as<MapVecd >(Lt[i]);
        int Ni = Ti.size();

        if(Ni > 1)
        {
            const MapMatd & Ri = as<MapMatd >(Lr[i]);
            const MapMatd & Bi = as<MapMatd >(B[i]);
            Eigen::MatrixXd tmp = Eigen::MatrixXd(-4*(Bi.transpose() * Ri * Bi * L));

            tmp = tmp + 4*(Bi.transpose() * Bi * C * Bi.transpose() * Bi * L);

            for(int j=0; j < Ni; j++)
            {
                const VectorXd & bij = Bi.row(j); // column vector
                tmp = tmp - 4 * (bij * bij.transpose() * C * bij * bij.transpose() * L);
                tmp = tmp + 4 * (Ri(j,j) * bij * bij.transpose() * L);
            }

            tmp = tmp * weight(i);

            H += tmp;
        }
    }

    //Rcpp::Rcout << "H=" << H << std::endl;

    // Rcpp::Rcout << "U=" << U << std::endl;
    // Rcpp::Rcout << "C=" << C << std::endl;
    // Rcpp::Rcout << "V=" << V << std::endl;
    // Rcpp::Rcout << "W=" << W << std::endl;
    // Rcpp::Rcout << "L=" << L << std::endl;
    // Rcpp::Rcout << "lam=" << lam << std::endl;
    // Rcpp::Rcout << "XX=" << 2 * lam * (U * C * W * L + W * C * U * L) + 4 * lam * V * C * V * L << std::endl;

    H +=  2 * lam * (U * C * W * L + W * C * U * L) + 4 * lam * V * C * V * L;

    return(H);
}


// [[Rcpp::export]]

Rcpp::List RhessianQ(const Rcpp::List & Lt,
                    const Rcpp::List & Lr,
                    const Eigen::MatrixXd& S,
                    const Rcpp::List & B,
                    const Eigen::MatrixXd& U,
                    const Eigen::MatrixXd& V,
                    const Eigen::MatrixXd& W,
                    const double lam,
                    const Eigen::VectorXd& weight,
                    const Eigen::MatrixXd& X)
{
    int n = Lt.size();
    const Eigen::MatrixXd & L = expM(X,S);

    int p = L.rows();
    const Eigen::MatrixXd & C = L * L.transpose();
    int d = p * (p+1) / 2;
    Eigen::MatrixXd H = Eigen::MatrixXd::Constant(d, d, 0);

    Eigen::MatrixXd G = RgradQ(Lt,Lr,L,B,U,V,W,lam,weight);

    double tmp1,tmp2,tmp3;
    Eigen::MatrixXd Eim;
    int j,k;
    double coef;

    List outList;

    for(int i=0; i<n; i++)
    {
        const MapVecd& Ti = as<MapVecd >(Lt[i]);
        int Ni = Ti.size();



        if(Ni > 1)
        {
            const MapMatd & Ri = as<MapMatd >(Lr[i]);
            const MapMatd & Bi = as<MapMatd >(B[i]);

            const Eigen::MatrixXd & E1 = Bi.transpose() * Ri * Bi;
            const Eigen::MatrixXd & E2 = Bi.transpose() * Bi;

            for(int j1=0; j1 < p; j1++)
            {
                for(int j2=0; j2 <= j1; j2++)
                {
                    j = j1*(j1+1)/2 + j2;
                    for(int k1=0; k1 < p; k1++)
                    {
                        for(int k2=0; k2 <= k1; k2++)
                        {
                            k = k1*(k1+1)/2 + k2;

                            if(j1==j2) coef = std::exp(S(j1,j2)/X(j1,j2));
                            else coef = 1;
                            if(k1==k2) coef *= std::exp(S(k1,k2)/X(k1,k2));

                            coef *= weight(i);

                            tmp1 = L.col(k2).transpose() * E2 * L.col(j2);
                            tmp2 = E2.row(j1).dot(L.col(k2)) * E2.row(k1).dot(L.col(j2));


                            if(j2==k2)
                            {
                                // Term I
                                H(j,k) -= coef*E1(j1,k1);

                                // Term II
                                tmp3 = E2.row(j1) * C * E2.col(k1);
                                H(j,k) += coef*(E2(j1,k1)*tmp1 + tmp2 + tmp3);

                            }
                            else
                            {
                                // Term I is zero
                                // Term II
                                H(j,k) += coef*(E2(j1,k1)*tmp1 + tmp2);
                            }

                            for(int m=0; m < Ni; m++)
                            {
                                Eim = Bi.row(m).transpose() * Bi.row(m);

                                tmp1 = L.col(k2).transpose() * Eim * L.col(j2);
                                tmp2 = (Eim.row(j1).dot(L.col(k2))) * (Eim.row(k1).dot(L.col(j2)));

                                if(j2 == k2)
                                {
                                    //Term III
                                    tmp3 = Eim.row(j1) * C * Eim.col(k1);
                                    H(j,k) -= coef*(Eim(j1,k1)*tmp1 + tmp2 + tmp3);

                                    //Term IV
                                    H(j,k) += coef*Ri(m,m)*Eim(j1,k1);

                                }
                                else
                                {
                                    //Term III
                                    H(j,k) -= coef*(Eim(j1,k1)*tmp1 + tmp2);

                                    //Term IV is zero
                                }
                            }

                        }
                    }
                }
            }
        }
    }

    H = 4*H;

    for(j=0; j<p; j++)
    {
        int jj = (j+2)*(j+1)/2 - 1;
        H(jj,jj) += G(j,j) * std::exp(S(j,j)/X(j,j)) / X(j,j);
    }

    for(int j1=0; j1 < p; j1++)
    {
        for(int j2=0; j2 <= j1; j2++)
        {
            j = j1*(j1+1)/2 + j2;
            for(int k1=0; k1 < p; k1++)
            {
                for(int k2=0; k2 <= k1; k2++)
                {
                    k = k1*(k1+1)/2 + k2;

                    if(j1==j2) coef = std::exp(S(j1,j2)/X(j1,j2));
                    else coef = 1;

                    if(j2==k2)
                    {
                        //Term I
                        //tmp1 = U(j1,k1)*(L.row(k2) * W * L.col(j2));
                        tmp1 = (double)(L.col(k2).transpose() * W * L.col(j2));
                        //if(j==0 && k==0 && j1==0 && j2==0 && k1==0 && k2==0)
                        //{
                        //    tmp1 = (double)(L.col(k2).transpose() * W * L.col(j2));
                            // Rcpp::Rcout << "tmp1=" << tmp1 << std::endl;
                            // Eigen::MatrixXd fuck = L.col(k2).transpose() * W * L.col(j2);
                            // Rcpp::Rcout << "fuck="
                            //            << fuck << std::endl;
                            // double fucktmp = L.col(k2).transpose() * W * L.col(j2);
                            // Rcpp::Rcout << "fucktmp="
                            //             << fucktmp << std::endl;
                            // Rcpp::Rcout << "L.col(k2).transpose() * W * L.col(j2)="
                            //             << L.col(k2).transpose() * W * L.col(j2) << std::endl;
                        //}
                        tmp1 = U(j1,k1) * tmp1;
                        tmp2 = U.row(j1).dot(L.col(k2)) * W.row(k1).dot(L.col(j2));
                        tmp3 = U.row(j1) * C * W.col(k1);
                        H(j,k) += 2*lam*coef*(tmp1+tmp2+tmp3);

                        // if(j==0 && k==0)
                        // {
                        //     Rcpp::Rcout << "U(j1,k1)=" << U(j1,k1) << std::endl;
                        //     Rcpp::Rcout << "L.col(k2).transpose() * W * L.col(j2)="
                        //                 << L.col(k2).transpose() * W * L.col(j2) << std::endl;
                        //     Rcpp::Rcout << "W=" << W << std::endl;
                        //     Rcpp::Rcout << "L=" << L << std::endl;
                        //     Rcpp::Rcout << "tmp1=" << tmp1 << std::endl;
                        //     Rcpp::Rcout << "tmp2=" << tmp2 << std::endl;
                        //     Rcpp::Rcout << "tmp3=" << tmp3 << std::endl;
                        // }

                        //Term II
                        tmp1 = L.col(k2).transpose() * U * L.col(j2);
                        tmp1 = W(j1,k1) * tmp1;
                        //tmp1 = W(j1,k1) * (L.col(k2).transpose() * U * L.col(j2));
                        tmp2 = W.row(j1).dot(L.col(k2)) * U.row(k1).dot(L.col(j2));
                        tmp3 = W.row(j1) * C * U.col(k1);
                        H(j,k) += 2*lam*coef*(tmp1+tmp2+tmp3);

                        //Term III
                        //tmp1 = V(j1,k1) * (L.col(k2).transpose() * V * L.col(j2));
                        tmp1 = L.col(k2).transpose() * V * L.col(j2);
                        tmp1 = V(j1,k1) * tmp1;
                        tmp2 = V.row(j1).dot(L.col(k2)) * V.row(k1).dot(L.col(j2));
                        tmp3 = V.row(j1) * C * V.col(k1);
                        H(j,k) += 4*lam*coef*(tmp1+tmp2+tmp3);
                    }
                    else
                    {
                        //Term I
                        //tmp1 = U(j1,k1)*(L.row(k2) * W * L.col(j2));
                        tmp1 = L.col(k2).transpose() * W * L.col(j2);
                        tmp1 = U(j1,k1) * tmp1;
                        tmp2 = U.row(j1).dot(L.col(k2)) * W.row(k1).dot(L.col(j2));
                        H(j,k) += 2*lam*coef*(tmp1+tmp2);

                        //Term II
                        tmp1 = L.col(k2).transpose() * U * L.col(j2);
                        tmp1 = W(j1,k1) * tmp1;
                        //tmp1 = W(j1,k1) * (L.col(k2).transpose() * U * L.col(j2));
                        tmp2 = W.row(j1).dot(L.col(k2)) * U.row(k1).dot(L.col(j2));
                        H(j,k) += 2*lam*coef*(tmp1+tmp2);

                        //Term III
                        //tmp1 = V(j1,k1) * (L.col(k2).transpose() * V * L.col(j2));
                        tmp1 = L.col(k2).transpose() * V * L.col(j2);
                        tmp1 = V(j1,k1) * tmp1;
                        tmp2 = V.row(j1).dot(L.col(k2)) * V.row(k1).dot(L.col(j2));
                        H(j,k) += 4*lam*coef*(tmp1+tmp2);
                    }
                }
            }
        }
    }
    outList["G"] = G;
    outList["H"] = H;
    return Rcpp::wrap(outList);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

// /*** R
// n <- 10
// Lt <- lapply(1:n, function(t) rnorm(5))
// Lr <- lapply(1:n, function(t) rnorm(5))
// X <- matrix(rnorm(5*5),5,5)
// S <- matrix(rnorm(5*5),5,5)
// Rhessian(Lt,Lr,X,X,X,S,X)
// */
