// Boost.AGeomAeAtry IndexAOA
//
// R-tree levels validating visitor implementation
//
// Copyright (c) 2011-2013 Adam Wulkiewicz, Lodz, Poland.
//
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#include <boost/geometry/index/detail/rtree/node/node.hpp>
#include <boost/foreach.hpp>

namespace bgi = boost::geometry::index;

namespace boost { namespace geometry { namespace index { namespace detail { namespace rtree { namespace utilities {

namespace visitors {

template <typename Value, typename Options, typename Box, typename Allocators>
class traversal_rtree
    : public rtree::visitor<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag, true>::type
{
    typedef typename rtree::internal_node<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag>::type internal_node;
    typedef typename rtree::leaf<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag>::type leaf;

public:
  int m_threshold;
  int m_dimension;
  double m_qs_lower;
  double m_qs_upper;
  int m_small_rank;
  int access_non_leaf;
  int access_leaf;
  bool m_all_in;
  bool stop;
  int result;
  std::vector<double> m_weight;
  std::vector<std::vector<double> > m_query_points;
  std::vector<std::vector<double> > m_mbr;
  std::vector<double> q_scores;

  inline traversal_rtree(int d,int t):m_small_rank(0),m_all_in(false), result(0),stop(false), access_non_leaf(0),access_leaf(0)
  {
    m_dimension = d;
    m_threshold = t;
  }

  /*non leaf node*/
    inline void operator()(internal_node const& n)
    {
        //IO
        access_non_leaf++;

        typedef typename rtree::elements_type<internal_node>::type elements_type;
        elements_type const& elements = rtree::elements(n);
        if(stop)
            return;

        if(m_all_in){
            m_small_rank += m_query_points.size() * (elements.size() - 1);
            if(m_small_rank >= m_threshold){
                stop = true;
                return;
            }
            for ( typename elements_type::const_iterator it = elements.begin();
                  it != elements.end() ; ++it){
                rtree::apply_visitor(*this, *it->second);
            }

        }
        else{
            for ( typename elements_type::const_iterator it = elements.begin();
                  it != elements.end() ; ++it)
                {
                    double lower_score =
                        boost::geometry::get<0>(it->first.min_corner()) * m_weight[0] +
                        boost::geometry::get<1>(it->first.min_corner()) * m_weight[1] +
                        boost::geometry::get<2>(it->first.min_corner()) * m_weight[2] +
                        boost::geometry::get<3>(it->first.min_corner()) * m_weight[3] +
                        boost::geometry::get<4>(it->first.min_corner()) * m_weight[4] ;

                    double upper_score =
                        boost::geometry::get<0>(it->first.max_corner()) * m_weight[0] +
                        boost::geometry::get<1>(it->first.max_corner()) * m_weight[1] +
                        boost::geometry::get<2>(it->first.max_corner()) * m_weight[2] +
                        boost::geometry::get<3>(it->first.max_corner()) * m_weight[3] +
                        boost::geometry::get<4>(it->first.max_corner()) * m_weight[4] ;


                    if (upper_score < m_qs_lower){
                        m_small_rank += m_query_points.size();
                        if(m_small_rank >= m_threshold){
                            stop = true;
                            return;
                        }
                        m_all_in = true;
                        rtree::apply_visitor(*this, *it->second);
                        m_all_in = false;
                    }
                    else if(lower_score > m_qs_upper){

                    }
                    else{
                        rtree::apply_visitor(*this, *it->second);
                    }

                }
        }
    }


    /*leaf node*/
    inline void operator()(leaf const& n)
    {
        //IO
        access_leaf++;

        typedef typename rtree::elements_type<leaf>::type elements_type;
        elements_type const& elements = rtree::elements(n);

        if(stop)
            return;

        if(m_all_in){
            m_small_rank += m_query_points.size() * (elements.size() - 1);
            if(m_small_rank >= m_threshold){
                stop = true;
                return;
            }
        }
        else{
            for ( typename elements_type::const_iterator it = elements.begin();
                  it != elements.end() ; ++it)
                {
                    double leaf_score =
                        boost::geometry::get<0>(it) * m_weight[0] +
                        boost::geometry::get<1>(it) * m_weight[1] +
                        boost::geometry::get<2>(it) * m_weight[2] +
                        boost::geometry::get<3>(it) * m_weight[3] +
                        boost::geometry::get<4>(it) * m_weight[4] ;

                    for(auto q_score : q_scores){
                        if(leaf_score < q_score){
                            m_small_rank++;
                            if(m_small_rank >= m_threshold){
                                stop = true;
                                return;
                            }
                        }
                    }
                }
        }
    }
};

} // namespace visitors
class Results
{
public:
    int access_leaf;
    int access_non_leaf;
    int rank;
};

template <typename Rtree> inline
Results traversal_rtree(Rtree const& tree, const std::vector<std::vector<double> > &qs, int t, int d, const std::vector<double> &w, const std::vector<std::vector<double> > &quplow, bool uplow)
{
    typedef utilities::view<Rtree> RTV;
    RTV rtv(tree);
    visitors::traversal_rtree<
        typename RTV::value_type,
        typename RTV::options_type,
        typename RTV::box_type,
        typename RTV::allocators_type
        > v(d,t);
    v.m_query_points = qs;
    v.m_weight = w;

    if(!uplow){
        auto mbr_qs = get_mbr(qs);
        v.m_mbr = mbr_qs;
        v.m_qs_lower = inner_product(mbr_qs[0],w);
        v.m_qs_upper = inner_product(mbr_qs[1],w);
    }else{
        v.m_qs_lower = inner_product(w, quplow[0]);
        v.m_qs_upper = inner_product(w, quplow[1]);
    }

    std::vector<double> qs_scores(qs.size());
    for(int i = 0; i < qs_scores.size();i++)
        qs_scores[i] = inner_product(qs[i],w);

    v.q_scores = qs_scores;
    rtv.apply_visitor(v);

    Results rs;
    if(v.stop)
        rs.rank = 0;
    else
        rs.rank = v.m_small_rank;

    rs.access_leaf = v.access_leaf;
    rs.access_non_leaf = v.access_non_leaf;
    return rs;
}

}}}}}} // namespace boost::geometry::index::detail::rtree::utilities







