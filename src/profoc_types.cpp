
#include <profoc_types.h>

namespace Rcpp
{

    template <>
    std::map<std::string, arma::colvec> as(SEXP matsexp)
    {
        Rcpp::NumericMatrix mat(matsexp);

        std::vector<std::string> cn = Rcpp::as<std::vector<std::string>>(Rcpp::colnames(mat));

        std::map<std::string, arma::colvec> mymap;

        for (int n = 0; n < mat.ncol(); n++)
        {
            mymap[cn[n]] = mat.column(n);
        }

        return mymap;
    }

    template <>
    SEXP wrap(const std::map<std::string, arma::colvec> &mymap)
    {
        Rcpp::NumericMatrix mat(mymap.begin()->second.n_elem, mymap.size());
        // Get all keys of the map
        std::vector<std::string> keys;
        for (auto const &x : mymap)
        {
            keys.push_back(x.first);
        }
        // Get all values of the map
        std::vector<arma::colvec> values;
        for (auto const &x : mymap)
        {
            values.push_back(x.second);
        }
        // Fill the matrix
        for (int n = 0; n < mat.ncol(); n++)
        {
            for (int m = 0; m < mat.nrow(); m++)
            {
                mat(m, n) = values[n](m);
            }
        }
        // Set column names
        Rcpp::CharacterVector colnames(keys.size());
        for (unsigned int n = 0; n < keys.size(); n++)
        {
            colnames[n] = keys[n];
        }
        Rcpp::colnames(mat) = colnames;

        return Rcpp::wrap(mat);
    }
}
