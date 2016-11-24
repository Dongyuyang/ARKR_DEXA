// Boost.Geometry Index
//
// R-tree levels validating visitor implementation
//
// Copyright (c) 2011-2013 Adam Wulkiewicz, Lodz, Poland.
//
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_INDEX_DETAIL_RTREE_UTILITIES_ARE_LEVELS_OK_HPP
#define BOOST_GEOMETRY_INDEX_DETAIL_RTREE_UTILITIES_ARE_LEVELS_OK_HPP

#include <boost/geometry/index/detail/rtree/node/node.hpp>


namespace boost { namespace geometry { namespace index { namespace detail { namespace rtree { namespace utilities {

namespace visitors {

template <typename Value, typename Options, typename Box, typename Allocators>
class are_levels_ok
    : public rtree::visitor<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag, true>::type
{
    typedef typename rtree::internal_node<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag>::type internal_node;
    typedef typename rtree::leaf<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag>::type leaf;
  typedef typename rtree::elements_type<internal_node>::type elements_type;

public:
    inline are_levels_ok()
        : result(true),
          m_leafs_level((std::numeric_limits<size_t>::max)()),
          m_current_level(0),number(0),
          is_ignore(false),is_single_w(false),
          num_E(0),flag(2)
    {}

    inline void operator()(internal_node const& n)
    {
        if(stop)
            return;

        typedef typename rtree::elements_type<internal_node>::type elements_type;
        elements_type const& elements = rtree::elements(n);

        if(m_all_in){
            num_E -= m_q_num;
            num_E += elements.size() * m_q_num;
            if(num_E >= k){
                flag = -1;
                stop = true;
                return;
            }
            for ( typename elements_type::const_iterator it = elements.begin();
                  it != elements.end() ; ++it){
                rtree::apply_visitor(*this, *it->second);
            }
        }else{
            for ( typename elements_type::const_iterator it = elements.begin();
                  it != elements.end() ; ++it)
                {
                    //version 2 not finished
                    m_up =
                        boost::geometry::get<0>(it->first.max_corner()) * w_up[0] +
                        boost::geometry::get<1>(it->first.max_corner()) * w_up[1] +
                        boost::geometry::get<2>(it->first.max_corner()) * w_up[2] ;

                    m_down =
                        boost::geometry::get<0>(it->first.min_corner()) * w_down[0] +
                        boost::geometry::get<1>(it->first.min_corner()) * w_down[1] +
                        boost::geometry::get<2>(it->first.min_corner()) * w_down[2] ;

                    if(m_up < q_down){
                        counter++;
                        num_E += m_q_num;
                        if(num_E >= k){
                            flag = -1;
                            return;
                        }
                        m_all_in = true;
                        rtree::apply_visitor(*this, *it->second);
                        m_all_in = false;
                    }
                    else if(m_down > q_up){
                        counter++;
                    }else{
                        rtree::apply_visitor(*this, *it->second);
                    }
                }//for
        }//else
    }

    inline void operator()(leaf const& n)
    {
        if(stop)
            return;
        typedef typename rtree::elements_type<leaf>::type elements_type;
        elements_type const& elements = rtree::elements(n);

        if(m_all_in){
            num_E += m_q_num * elements.size();
            if(num_E >= k){
                flag = -1;
                stop = true;
                return;
            }
        }else{
            for ( typename elements_type::const_iterator it = elements.begin();
                  it != elements.end() ; ++it)
                {
                    p_up =
                        boost::geometry::get<0>(it) * w_up[0] +
                        boost::geometry::get<1>(it) * w_up[1] +
                        boost::geometry::get<2>(it) * w_up[2] ;

                    if(p_up < q_down){
                        num_E++;
                        if(num_E >= k){
                            flag = -1;
                            return;
                        }
                    }
                }//for
        }

    }

    bool result;
    //add by tou
    bool is_ignore;
    bool rtk;
    int flag;
    int k,number,num_E;
    double q_score;
    std::vector<double> w;

    //new2
    bool is_single_w;
    double q_down;
    double q_up;
    std::vector<double> w_up;
    std::vector<double> w_down;
    bool is_frist;
    double m_up,m_down,p_up,p_down;
    int m_q_num;
    bool m_all_in = false;
    bool stop = false;
    int counter = 0;
private:
    size_t m_leafs_level;
    size_t m_current_level;
};

} // namespace visitors

template <typename Rtree> inline
int  are_levels_ok(Rtree const& tree,
                   const std::vector<double> &vector_up,
                   const std::vector<double> &vector_down,
                   const std::vector<double> &qs_up,
                   const std::vector<double> &qs_down,
                   int kk,
                   int q_num)
{
    typedef utilities::view<Rtree> RTV;
    RTV rtv(tree);

    visitors::are_levels_ok<
        typename RTV::value_type,
        typename RTV::options_type,
        typename RTV::box_type,
        typename RTV::allocators_type
        > v;

    v.rtk = true;
    v.k = kk;
    v.w_up = vector_up;
    v.w_down = vector_down;
    v.m_q_num = q_num;

    double qq_up = 0, qq_down = 0;
    for(int i = 0; i < qs_up.size();i++){
        qq_up += qs_up[i] * vector_up[i];
        qq_down += qs_down[i] * vector_down[i];
    }

    v.q_up = qq_up;
    v.q_down = qq_down;

    if((vector_up[0] == vector_down[0]) && (vector_up[1] == vector_down[1]))
        v.is_single_w = true;

    rtv.apply_visitor(v);

    if(v.flag == 2){
        if(v.num_E < kk)
            v.flag = 1;
        else
            v.flag = 0;
    }
    //return v.result;
    return v.flag;
    //return 0;

}

//new version 2

}}}}}} // namespace boost::geometry::index::detail::rtree::utilities

#endif // BOOST_GEOMETRY_INDEX_DETAIL_RTREE_UTILITIES_ARE_LEVELS_OK_HPP
