// Boost.Geometry Index
//
// R-tree levels validating visitor implementation
//
// Copyright (c) 2011-2013 Adam Wulkiewicz, Lodz, Poland.
//
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef _ARE_LEVELS_OK_HPP_
#define _ARE_LEVELS_OK_HPP_

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
        : 
          m_leafs_level((std::numeric_limits<size_t>::max)()),
          m_current_level(0),
          m_small_rank(0),flag(2)
    {}

    inline void operator()(internal_node const& n)
    {
        if(stop)
            return;

        typedef typename rtree::elements_type<internal_node>::type elements_type;
        elements_type const& elements = rtree::elements(n);

        if(m_all_in){
            m_small_rank +=
                (elements.size() - 1) * m_q_cluster_bound[current_cindex][2];
            if(m_small_rank >= m_threshold){
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
                    double upper_score =
                        boost::geometry::get<0>(it->first.max_corner()) * w_up[0] +
                        boost::geometry::get<1>(it->first.max_corner()) * w_up[1] +
                        boost::geometry::get<2>(it->first.max_corner()) * w_up[2] ;

                    double lower_score =
                        boost::geometry::get<0>(it->first.min_corner()) * w_down[0] +
                        boost::geometry::get<1>(it->first.min_corner()) * w_down[1] +
                        boost::geometry::get<2>(it->first.min_corner()) * w_down[2] ;

                    if(single_mode){
                        if(upper_score < m_q_cluster_bound[current_cindex][0]){
                            m_small_rank += m_q_cluster_bound[current_cindex][2];
                            if(m_small_rank >= m_threshold){
                                flag = -1;
                                stop = true;
                                return;
                            }
                            m_all_in = true;
                            rtree::apply_visitor(*this, *it->second);
                            m_all_in = false;
                        }

                        if(upper_score >= m_q_cluster_bound[current_cindex][0] &&
                               lower_score <= m_q_cluster_bound[current_cindex][1])
                                {
                                    rtree::apply_visitor(*this, *it->second);
                                }
                    } else {
                        for(int nc = 0; nc < m_q_cluster_bound.size(); nc++){
                            if(upper_score < m_q_cluster_bound[nc][0]){
                                m_small_rank += m_q_cluster_bound[nc][2];
                                if(m_small_rank >= m_threshold){
                                    flag = -1;
                                    stop = true;
                                    return;
                                }
                                m_all_in = true;
                                current_cindex = nc;
                                rtree::apply_visitor(*this, *it->second);
                                m_all_in = false;
                            }

                            if(upper_score >= m_q_cluster_bound[nc][0] &&
                               lower_score <= m_q_cluster_bound[nc][1])
                                {
                                    single_mode = true;
                                    current_cindex = nc;
                                    rtree::apply_visitor(*this, *it->second);
                                    single_mode = false;
                                }
                        }
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
            m_small_rank +=
                m_q_cluster_bound[current_cindex][2] * (elements.size() - 1);
            if(m_small_rank >= m_threshold){
                flag = -1;
                stop = true;
                return;
            }
        }else{
            for ( typename elements_type::const_iterator it = elements.begin();
                  it != elements.end() ; ++it)
                {
                    double p_up =
                        boost::geometry::get<0>(it) * w_up[0] +
                        boost::geometry::get<1>(it) * w_up[1] +
                        boost::geometry::get<2>(it) * w_up[2] ;

                    bool temp_break = false;
                    if(single_mode){
                        for(auto &q_score : m_q_cluster_scores[current_cindex]){
                            if(p_up < q_score){
                                m_small_rank++;
                                if(m_small_rank >= m_threshold){
                                    flag = -1;
                                    stop = true;
                                    temp_break = true;
                                    break;
                                }
                            }
                        }
                    } else {
                        for(auto &&q_c : m_q_cluster_scores){
                            for(auto q_score : q_c){
                                if(p_up < q_score){
                                    m_small_rank++;
                                    if(m_small_rank >= m_threshold){
                                        flag = -1;
                                        stop = true;
                                        temp_break = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if(temp_break)
                        break;
                }//for
        }

    }

    int flag;
    int m_threshold,m_small_rank;
    int current_cindex;
    std::vector<double> w;
    std::vector<double> w_up,w_down;
    bool m_all_in = false, stop = false, single_mode = false;
    std::vector<std::vector<double> > m_q_cluster_bound;
    std::vector<std::vector<double> > m_q_cluster_scores;

    int counter = 0;
private:
    size_t m_leafs_level, m_current_level;
};

} // namespace visitors

template <typename Rtree> inline
int  are_levels_ok(Rtree const& tree,
                   int t,
                   const std::vector<double> &vector_up,
                   const std::vector<double> &vector_down,
                   const std::vector<std::vector<std::vector<double> > > &q_clusters
                   )
{
    typedef utilities::view<Rtree> RTV;
    RTV rtv(tree);

    visitors::are_levels_ok<
        typename RTV::value_type,
        typename RTV::options_type,
        typename RTV::box_type,
        typename RTV::allocators_type
        > v;

    v.m_threshold = t;
    v.w_up = vector_up;
    v.w_down = vector_down;

    v.m_q_cluster_bound.resize(q_clusters.size());
    v.m_q_cluster_scores.resize(q_clusters.size());
    int cindex = 0;
    for(auto &c : q_clusters){
        auto mbr = get_mbr(c);
        v.m_q_cluster_bound[cindex].push_back(inner_product(mbr[0],vector_down));
        v.m_q_cluster_bound[cindex].push_back(inner_product(mbr[1],vector_up));
        v.m_q_cluster_bound[cindex].push_back(c.size());
        for(int i = 0; i < c.size(); i++)
          v.m_q_cluster_scores[cindex].push_back(inner_product(c[i],vector_down));
        cindex++;
    }
    v.current_cindex = q_clusters.size();

    /*apply*/
    rtv.apply_visitor(v);

    if(v.flag == 2){
        if(v.m_small_rank < t)
            v.flag = 1;
        else
            v.flag = 0;
    }

    return v.flag;
    //return 0;

}

//new version 2

}}}}}} // namespace boost::geometry::index::detail::rtree::utilities


#endif
