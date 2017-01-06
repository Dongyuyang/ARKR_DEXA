#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <vector>
#include <iostream>
#include <boost/foreach.hpp>
#include <iostream>
#include <stdlib.h>
#include <cstdlib>
#include <time.h>
#include <algorithm>
#include "model.hpp"
#include "are_levels_ok.hpp"

#define D 5

namespace alo = boost::geometry::index::detail::rtree::utilities;

void rt::creat_rtree_p_w(std::vector<Point> &nodes, std::vector<Weight> &vectors)
{
  // creat rtree
  int dimension = nodes[0].value.size();
  if(dimension == 2){
    for(int i = 0;i < nodes.size(); ++i){
      point p1(nodes[i].value[0], nodes[i].value[1]);
      rtree.insert(p1);
    }
    for(int i = 0; i < vectors.size(); i++){
      point p2(vectors[i].value[0], vectors[i].value[1]);
      rtree_w.insert(std::make_pair(p2,vectors[i].id));
    }
  }

  if(dimension == 3){
    for(int i = 0;i < nodes.size(); ++i){
      point p1(nodes[i].value[0], nodes[i].value[1], nodes[i].value[2]);
      rtree.insert(p1);
    }
    for(int i = 0 ;i < vectors.size(); ++i){
      point p2(vectors[i].value[0], vectors[i].value[1], vectors[i].value[2]);
      rtree_w.insert(std::make_pair(p2,vectors[i].id));
    }
  }
  if(dimension > 3){
    //S
    point p1;
    //ddd
    for(int i = 0;i < nodes.size(); ++i){
      bg::set<0>(p1, nodes[i].value[0]);
      bg::set<1>(p1, nodes[i].value[1]);
      bg::set<2>(p1, nodes[i].value[2]);
      bg::set<3>(p1, nodes[i].value[3]);
      bg::set<4>(p1, nodes[i].value[4]);
      /*bg::set<5>(p1, nodes[i].value[5]);
      bg::set<6>(p1, nodes[i].value[6]);
      bg::set<7>(p1, nodes[i].value[7]);
      bg::set<8>(p1, nodes[i].value[8]);
      bg::set<9>(p1, nodes[i].value[9]);
      bg::set<10>(p1, nodes[i].value[10]);
      bg::set<11>(p1, nodes[i].value[11]);
      bg::set<12>(p1, nodes[i].value[12]);
      bg::set<13>(p1, nodes[i].value[13]);
      bg::set<14>(p1, nodes[i].value[14]);
      bg::set<15>(p1, nodes[i].value[15]);
      bg::set<16>(p1, nodes[i].value[16]);
      bg::set<17>(p1, nodes[i].value[17]);
      bg::set<18>(p1, nodes[i].value[18]);
      bg::set<19>(p1, nodes[i].value[19]);
      bg::set<20>(p1, nodes[i].value[20]);
      bg::set<21>(p1, nodes[i].value[21]);
      bg::set<22>(p1, nodes[i].value[22]);
      bg::set<23>(p1, nodes[i].value[23]);
      bg::set<24>(p1, nodes[i].value[24]);
      bg::set<25>(p1, nodes[i].value[25]);*/
      rtree.insert(p1);
    }
      
    //W
    point p2;
    for(int i = 0;i < vectors.size(); ++i){
      bg::set<0>(p2, vectors[i].value[0]);
      bg::set<1>(p2, vectors[i].value[1]);
      bg::set<2>(p2, vectors[i].value[2]);
      bg::set<3>(p2, vectors[i].value[3]);
      bg::set<4>(p2, vectors[i].value[4]);
      /*bg::set<5>(p2, vectors[i].value[5]);
      bg::set<6>(p2, vectors[i].value[6]);
      bg::set<7>(p2, vectors[i].value[7]);
      bg::set<8>(p2, vectors[i].value[8]);
      bg::set<9>(p2, vectors[i].value[9]);
      bg::set<10>(p2, vectors[i].value[10]);
      bg::set<11>(p2, vectors[i].value[11]);
      bg::set<12>(p2, vectors[i].value[12]);
      bg::set<13>(p2, vectors[i].value[13]);
      bg::set<14>(p2, vectors[i].value[14]);
      bg::set<15>(p2, vectors[i].value[15]);
      bg::set<16>(p2, vectors[i].value[16]);
      bg::set<17>(p2, vectors[i].value[17]);
      bg::set<18>(p2, vectors[i].value[18]);
      bg::set<19>(p2, vectors[i].value[19]);
      bg::set<20>(p2, vectors[i].value[20]);
      bg::set<21>(p2, vectors[i].value[21]);
      bg::set<22>(p2, vectors[i].value[22]);
      bg::set<23>(p2, vectors[i].value[23]);
      bg::set<24>(p2, vectors[i].value[24]);
      bg::set<25>(p2, vectors[i].value[25]);*/
      rtree_w.insert(std::make_pair(p2,vectors[i].id));
    }
  }
}
