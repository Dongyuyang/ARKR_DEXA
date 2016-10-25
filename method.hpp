#include <vector>
#include <iterator>
#include <limits>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "traversal_rtree.hpp"

/*count the number of p rank better than q*/
static int counter_rank_p_q(double p_score, const std::vector<double> &q_scores)
{
  int counter = 0;
  for(auto q_score : q_scores){
      if(p_score < q_score)
	counter += 1;
  }
  return counter;
}


/*Naive top-k rank*/
int agg_rank_navie(const std::vector<Point> &p, const Weight &w, const std::vector<std::vector<double> > &qs)
{
  std::vector<double> q_scores(qs.size());
  for(int i = 0; i < qs.size();i++)
    q_scores[i] = inner_product(qs[i], w.value);

  int agg_rank = 0;
  for(int i = 0; i < p.size();i++){
    double p_score = inner_product(p[i].value,w.value);
    agg_rank += counter_rank_p_q(p_score,q_scores);
  }
  return agg_rank;
}

/*Naive arkr*/
std::multimap<int,int> naive_arkr(const std::vector<Point> &p, const std::vector<Weight> &w, const std::vector<std::vector<double> > &qs, int k)
{
  std::multimap<int,int> id_rank;
  for(int i = 0; i < w.size();i++){
    int c_rank = agg_rank_navie(p,w[i],qs);
    id_rank.insert(std::pair<int,int>(c_rank,w[i].id));
  }

  std::multimap<int,int> result;
  std::multimap<int,int>::iterator it = id_rank.begin();
  int time = k;
  while( time-- >= 1 )
      result.insert(*it++);

  return result;
}

void print_map(const std::multimap<int,int>& map) {
  for ( std::multimap<int,int>::const_iterator it = map.begin(); it != map.end(); it++) {
    std::cout << "ARank: " << it->first;
    std::cout << " ,Id: " << it->second <<std::endl;
  }
}


/*TPM*/
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;
namespace alo = boost::geometry::index::detail::rtree::utilities;
typedef bg::model::point<float, 3, bg::cs::cartesian> r_point;

std::multimap<int,int> tpm(bgi::rtree< r_point, bgi::rstar<16> > &rtree, const std::vector<Point> &p, const std::vector<Weight> &w, const std::vector<std::vector<double> > &qs, int k, const std::vector<std::vector<double> > &quplow, bool uplow )
{
  std::multimap<int,int> buffer;
  int dimension = qs[0].size();
  double compute_times = 0;
  double non_leaf = 0;
  int thresold = std::numeric_limits<int>::max();
  for(int i = 0; i < w.size(); i++){
      if(buffer.size() == k)
          int thresold = std::prev(buffer.end())->first;
      auto rs =
          alo::traversal_rtree(rtree,qs,thresold,dimension,w[i].value,quplow,uplow);
      compute_times += rs.access_leaf;
      non_leaf += rs.access_non_leaf;
      if(rs.rank == 0){
          continue;
      }
      else{
          if(buffer.size() == k)
              update_buffer(buffer,rs.rank,w[i].id);
          else
              buffer.insert(std::pair<int,int>(rs.rank,w[i].id));
      }
  }

  std::cout << "TPM: non leaf times: " << non_leaf  << std::endl;
  std::cout << "TPM: leaf times: " << compute_times  << std::endl;
  return buffer;
}

static std::vector<double> most_sim_w(const std::vector<Weight> &w, int direcet_index, int dim)
{
  std::vector<double> direct(dim,0);
  direct[direcet_index] = 1;
  double max = -1;
  int id = 0;
  for(int i = 0; i < w.size();i++){
    double sim = cos_sim(w[i].value,direct);
    if(sim > max){
      max = sim;
      id = w[i].id;
    }
  }

  return w[id].value;

}

static std::vector<std::vector<double> > most_sim_w_all_dim(const std::vector<Weight> &w, int dim)
{
  std::vector<std::vector<double> > result(dim);
  for(int i = 0; i < dim; i++)
    result[i] = most_sim_w(w,i,dim);
  return result;
}


std::vector<std::vector<double> > get_Q_low_up(const std::vector<Weight> &w, const std::vector<std::vector<double> > &q, int dim)
{
    auto w_t = most_sim_w_all_dim(w,dim);
    std::vector<std::vector<double> > uppers(dim);
    std::vector<std::vector<double> > lowers(dim);
    for(int i = 0; i < dim; i++){
        double max = -1, min = 10000000;
        int upper_q = 0, lower_q = 0;
        for(int j = 0; j < q.size();j++){
            double score = inner_product(w_t[i],q[j]);
            if(score > max){
                max = score;
                upper_q = j;
            }

            if(score < min){
                min = score;
                lower_q = j;
            }
        }
        uppers[i] = q[upper_q];
        lowers[i] = q[lower_q];
    }
    return {get_mbr(lowers)[0], get_mbr(uppers)[1]};
}














