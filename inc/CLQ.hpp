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
class CLQ
    : public rtree::visitor<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag, true>::type
{
    typedef typename rtree::internal_node<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag>::type internal_node;
    typedef typename rtree::leaf<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag>::type leaf;

public:
    int access_non_leaf;int access_leaf;
    int m_threshold;
    int m_small_rank; int current_cindex = 0;int now_size;
    double m_qs_lower;double m_qs_upper;
    bool m_all_in;bool stop;
    std::vector<double> m_weight;
    std::vector<std::vector<double> > m_query_points;
    std::vector<double> q_scores;
    std::vector<std::vector<double> > m_q_cluster_bound;
    std::vector<std::vector<double> > m_q_cluster_scores;


    inline CLQ(int t):
        m_small_rank(0),m_all_in(false),stop(false),
        access_non_leaf(0),access_leaf(0),m_threshold(t) {}

  /*non leaf node*/
    inline void operator()(internal_node const& n)
    {
        //access_non_leaf++;
        access_non_leaf += sizeof(n);
        typedef typename rtree::elements_type<internal_node>::type elements_type;
        elements_type const& elements = rtree::elements(n);
        if(stop)
            return;

        if(m_all_in){
            m_small_rank +=
                m_q_cluster_bound[current_cindex][2] * (elements.size() - 1);
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
                        boost::geometry::get<2>(it->first.min_corner()) * m_weight[2] ;
                    double upper_score =
                        boost::geometry::get<0>(it->first.max_corner()) * m_weight[0] +
                        boost::geometry::get<1>(it->first.max_corner()) * m_weight[1] +
                        boost::geometry::get<2>(it->first.max_corner()) * m_weight[2] ;


                    for(int nc = 0; nc < m_q_cluster_bound.size(); nc++){
                        if(upper_score < m_q_cluster_bound[nc][0]){
                            m_small_rank += m_q_cluster_bound[nc][2];
                            if(m_small_rank >= m_threshold){
                                stop = true;
                                return;
                            }
                            m_all_in = true;
                            current_cindex = nc;
                            rtree::apply_visitor(*this, *it->second);
                            m_all_in = false;
                        }
                    }

                    if(upper_score >= m_qs_lower && lower_score <= m_qs_upper){
                        m_all_in = false;
                        current_cindex = -1;
                        rtree::apply_visitor(*this, *it->second);
                    }



                }
        }
    }


    /*leaf node*/
    inline void operator()(leaf const& n)
    {
        //access_leaf++;
        access_leaf += sizeof(n);
        typedef typename rtree::elements_type<leaf>::type elements_type;
        elements_type const& elements = rtree::elements(n);
        if(stop)
            return;

        if(m_all_in){
            m_small_rank +=
                m_q_cluster_bound[current_cindex][2] * (elements.size() - 1);
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
                        boost::geometry::get<2>(it) * m_weight[2] ;

                    bool temp_break = false;
                    if(current_cindex != -1){
                        for(auto q_score : m_q_cluster_scores[current_cindex]){
                            if(leaf_score < q_score){
                                m_small_rank++;
                                if(m_small_rank >= m_threshold){
                                    stop = true;
                                    temp_break = true;
                                    break;
                                }
                            }
                        }
                    }
                    else{
                        for(auto &q_c : m_q_cluster_scores){
                            for(auto &q_score : q_c){
                                if(leaf_score < q_score){
                                    m_small_rank++;
                                    if(m_small_rank >= m_threshold){
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
                }
        }
    }
};

} // namespace visitors


/*class Results
 {
public:
    int access_leaf;
    int access_non_leaf;
    int rank;
    };*/

template <typename Rtree> inline
Results CLQ(Rtree const& tree,
            const std::vector<std::vector<double> > &qs,
            int t,
            const std::vector<double> &w,
            const std::vector<std::vector<std::vector<double> > > &q_clusters
            )
{
    typedef utilities::view<Rtree> RTV;
    RTV rtv(tree);
    visitors::CLQ<
        typename RTV::value_type,
        typename RTV::options_type,
        typename RTV::box_type,
        typename RTV::allocators_type
        > v(t);
    v.m_query_points = qs;
    v.m_weight = w;

    auto qmbr = get_mbr(qs);
    v.m_qs_lower = inner_product(qmbr[0],w);
    v.m_qs_upper = inner_product(qmbr[1],w);

    v.m_q_cluster_bound.resize(q_clusters.size());
    v.m_q_cluster_scores.resize(q_clusters.size());

    int cindex = 0;
    for(auto &c : q_clusters){
        auto mbr = get_mbr(c);
        v.m_q_cluster_bound[cindex].push_back(inner_product(mbr[0],w));
        v.m_q_cluster_bound[cindex].push_back(inner_product(mbr[1],w));
        v.m_q_cluster_bound[cindex].push_back(c.size());
        for(int i = 0; i < c.size(); i++)
          v.m_q_cluster_scores[cindex].push_back(inner_product(c[i],w));
        cindex++;
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







