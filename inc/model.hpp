#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
//#include "model.hpp"

#define D 3

/* Point: the product p (id,vector)*/
class Point
{
public:
  int id;
  std::vector<double> value;
  Point(){};
  Point(int i,std::vector<double> &v){id = i; value = v;};
};

/* Weight: the user perference w (id,vector).*/
class Weight
{
public:
  int id;
  std::vector<double> value;
  Weight(){};
  Weight(int i, std::vector<double> v){id = i; value = v;};
};

/*rtree*/
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<float, D, bg::cs::cartesian> point;
typedef std::pair<point, int> value;
/*rtree class*/
class rt{
public:
  rt(){}
  void creat_rtree_p_w(std::vector<Point> &nodes, std::vector<Weight> &vectors);

public:
  bgi::rtree< point, bgi::rstar<16> > rtree;
  bgi::rtree< value, bgi::rstar<16> > rtree_w;

};




