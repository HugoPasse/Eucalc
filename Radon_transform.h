//#include "Embedded_cubical_complex.h"
#include <vector>

#include <stdexcept>
#include <limits>


#ifndef RADON_TRANSFORM_H_
#define RADON_TRANSFORM_H_

class Radon_transform{
    public :
        std::vector<double> T;
        std::vector<double> Values;

        std::vector<double> singular_T;
        std::vector<double> singular_values;

        std::size_t dichotomie(std::vector<double> table, double t, std::size_t begin, std::size_t end){
            if(end - begin <= 1){
                return begin;
            }else{
                std::size_t index = (begin + end) / 2;
                if(table[index] > t){
                    return dichotomie(table, t, begin, index);
                }else{
                    return dichotomie(table, t, index, end);
                }
            }
        }

        double get_singular_critical_value(double t, bool* is_value_critical){
            std::size_t index = dichotomie(singular_T, t, 0, singular_T.size());
            if(std::abs(singular_T[index] - t) <= std::numeric_limits<double>::epsilon()){
                *is_value_critical = true;
                return singular_values[index];
            }else if(index < singular_T.size() - 1 && std::abs(singular_T[index + 1] - t) <= std::numeric_limits<double>::epsilon()){
                *is_value_critical = true;
                return singular_values[index + 1];
            }
            *is_value_critical = false;
            return 0;
        }

        Radon_transform(std::vector<double> _T, std::vector<double> _Values, std::vector<double> _singular_T, std::vector<double> _singular_Values){
            T = _T;
            Values = _Values;
            singular_T = _singular_T;
            singular_values = _singular_Values;
        }
   
        double evaluate(double t){
            bool is_value_critical = false;
            double val = get_singular_critical_value(t, &is_value_critical);

            if(is_value_critical){
                return val;
            }

            if(t < T[0]){
                return 0;
            }else{
                int index = dichotomie(T,t,0,T.size());
                return Values[index];
            }
        }

        std::vector<std::vector<double>> get_attributes(){
            std::vector<std::vector<double>> result;

            result.push_back(T);
            result.push_back(Values);
            result.push_back(singular_T);
            result.push_back(singular_values);

            return result;
        }
};

#endif //define RADON_TRANSFORM_H_