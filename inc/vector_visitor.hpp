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
#include "are_levels_ok.hpp"


namespace bgi = boost::geometry::index;
namespace alo = boost::geometry::index::detail::rtree::utilities;
typedef bg::model::point<float, 3, bg::cs::cartesian> point;

namespace boost { namespace geometry { namespace index { namespace detail { namespace rtree { namespace utilities {
namespace visitors {

template <typename Value, typename Options, typename Box, typename Allocators>
class vector_visitor
    : public rtree::visitor<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag, true>::type
{
    typedef typename rtree::internal_node<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag>::type internal_node;
    typedef typename rtree::leaf<Value, typename Options::parameters_type, Box, Allocators, typename Options::node_tag>::type leaf;

public:
    inline vector_visitor()
      : result(true), is_ignore(false)
    {}

    inline void operator()(internal_node const& n)
    {
      w_non_leaf ++;
      typedef typename rtree::elements_type<internal_node>::type elements_type;
      elements_type const& elements = rtree::elements(n);

      for ( typename elements_type::const_iterator it = elements.begin();
	    it != elements.end() ; ++it)
	{
	  if(is_ignore){
	    rtree::apply_visitor(*this, *it->second);
	  }
	  else{
	    /*w_min*/
        /*change dim jianghan*/
	    w_min[0] = boost::geometry::get<0>(it->first.min_corner());
	    w_min[1] = boost::geometry::get<1>(it->first.min_corner());
	    w_min[2] = boost::geometry::get<2>(it->first.min_corner());
	    /*w_max*/
        /*change dim jianghan*/
	    w_max[0] = boost::geometry::get<0>(it->first.max_corner());
	    w_max[1] = boost::geometry::get<1>(it->first.max_corner());
	    w_max[2] = boost::geometry::get<2>(it->first.max_corner());

        if(m_buffer.size() == k)
            thresold = std::prev(m_buffer.end())->first;

	    int temp = alo::are_levels_ok(rtr,w_max,w_min,qs_mbr[1],qs_mbr[0],thresold,m_qs.size());

	    if(temp == 1){
	      is_ignore = true;
	      rtree::apply_visitor(*this, *it->second);
	      is_ignore = false;
	    }
	    if(temp == 0){
	      is_ignore = false;
	      rtree::apply_visitor(*this, *it->second);
	    }
	  }
	}//for

    }

    inline void operator()(leaf const& n)
    {
      w_leaf++;
      typedef typename rtree::elements_type<leaf>::type elements_type;
      elements_type const& elements = rtree::elements(n);

      for ( typename elements_type::const_iterator it = elements.begin();
	    it != elements.end() ; ++it)
        {
            /*change dim jianghan*/
            w[0] = boost::geometry::get<0>(it->first);
            w[1] = boost::geometry::get<1>(it->first);
            w[2] = boost::geometry::get<2>(it->first);

            if(m_buffer.size() >= k)
                thresold = std::prev(m_buffer.end())->first;

            /*min max*/
            std::vector<double> q;
            double min = std::numeric_limits<double>::max();
            for(auto &_q : m_qs){
                double q_s = inner_product(_q,w);
                if(q_s < min){
                    min = q_s;
                    q = _q;
                }
            }

            auto rs =
                alo::traversal_rtree(rtr, {q},
                                     thresold, m_qs[0].size(),w,{{1,2}},false);
            non_leaf_times += rs.access_non_leaf;
            compute_times += rs.access_leaf;

            if(rs.rank != 0){
                if(m_buffer.size() >= k)
                    update_buffer(m_buffer,rs.rank,it->second);
                else
                    m_buffer.insert(std::pair<int,int>(rs.rank,it->second));
            }

	}

    }

private:
  bool result;
public:
  std::vector<std::vector<double> > m_qs;
  bool is_ignore;
  bgi::rtree< point, bgi::rstar<16> > rtr;
  int k;
  int thresold;
  std::multimap<int,int> m_buffer;
  std::vector<std::vector<double> > qs_mbr;
  std::vector<double> w_min, w_max, w;
  double non_leaf_times = 0;
  double compute_times = 0;
  double w_non_leaf = 0;
  double w_leaf = 0;
};

} // namespace visitors

template <typename Rtree> inline
std::multimap<int,int> vector_visitor(
                                      Rtree const& tree,
                                      const std::vector<std::vector<double> > &qs,
                                      bgi::rtree< point, bgi::rstar<16> > &rte,
                                      int current_rank,
                                      int result_num)
{
  typedef utilities::view<Rtree> RTV;
  RTV rtv(tree);
  visitors::vector_visitor<
    typename RTV::value_type,
    typename RTV::options_type,
    typename RTV::box_type,
    typename RTV::allocators_type
    > v;

  v.qs_mbr = get_mbr(qs);
  v.m_qs = qs;
  v.rtr = rte;
  v.thresold = current_rank;
  v.k = result_num;
  int dimension = qs[0].size();
  v.w_min.resize(dimension);
  v.w_max.resize(dimension);
  v.w.resize(dimension);

  rtv.apply_visitor(v);

  std::cout << "DTM: P non leaf times: " << v.non_leaf_times  << std::endl;
  std::cout << "DTM: P leaf times: " << v.compute_times  << std::endl;

  std::cout << "DTM: W non leaf times: " << v.w_non_leaf << std::endl;
  std::cout << "DTM: W leaf times: " << v.w_leaf  << std::endl;


  return v.m_buffer;
}

}}}}}} // namespace boost::geometry::index::detail::rtree::utilities



