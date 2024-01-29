#import "utils/quantile.h"
#import <iostream>
int main(){
    std::vector<double> test1{1,4,2,3,9,6,7,8};
    std::vector<double> probs{0.25, 0.75};

    std::vector<double> output = Quantile(test1,probs);

    for (auto j : output){
        std::cout << j << "\n";
    }
    return 0;
};



