/*based on polygon Boost library*/
#ifndef _PLOYGON_HPP_
#define _PLOYGON_HPP_

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>
#include <vector>
#include <iostream>

/*change dim jianghan*/
#define D 3

namespace dyy{
namespace poly{

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
typedef bg::model::point<float, D, bg::cs::cartesian> point_poly;
typedef bg::model::polygon<point_poly> polygon;

class ConvexHull
{
public:
    ConvexHull(const std::vector<std::vector<double> > &Q)
    {
        make_convex_hull(Q);
    }

    /*get polygon convex hull*/
    polygon get(){return m_hull;}

    /*get vectexes transfer to std::vector*/
    std::vector<std::vector<double> > get_points()
    {
        std::vector<std::vector<double> > vectexes;
        auto points = m_hull.outer();
        for(auto &p : points){
            std::vector<double> v(D);
            /*change dim jianghan*/
            v[0] = bg::get<0>(p);
            v[1] = bg::get<1>(p);
            v[2] = bg::get<2>(p);
            vectexes.push_back(v);
        }
        return vectexes;
    }

    void print(){std::cout << "hull: " << boost::geometry::dsv(m_hull) << std::endl;}

private:
    polygon make_polygon(const std::vector<std::vector<double> > &Q)
    {
        polygon poly;
        for(auto &q : Q){
            //D times
            point_poly p;
            /*change dim jianghan*/
            bg::set<0>(p,q[0]);
            bg::set<1>(p,q[1]);
            bg::set<2>(p,q[2]);
            bg::append(poly,p);
        }
        return poly;
    }

    void make_convex_hull(const std::vector<std::vector<double> > &Q)
    {
        polygon poly = make_polygon(Q);
        boost::geometry::convex_hull(poly, m_hull);
    }

private:
    polygon m_hull;
};

}}// dyy::poly


#endif /*_POLYGON_HPP_*/
