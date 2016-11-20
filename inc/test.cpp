#include <iostream>
#include "polygon.hpp"
#include <vector>
#include "../../commontool/common_tool.hpp"

int main()
{
    std::vector<std::vector<double> > Q = {
        {2.0,1.3},
        {2.4,1.7},
        {2.8,1.8},
        {3.4,1.2},
        {3.7,1.6},
        {3.4,2.0},
        {4.1,3.0},
        {5.3,2.6},
        {5.4,1.2},
        {4.9,0.8},
        {2.9,0.7},
        {2.0,1.3}
    };
    dyy::poly::ConvexHull CH(Q);
    CH.print();

    auto vectexes = CH.get_points();
    std::cout << "v size: " << vectexes.size() << std::endl;

    for(auto &q : vectexes){
        put_vector(q);
    }


    return 0;
}
