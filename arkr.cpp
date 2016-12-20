#include "../commonTool/catch.h"
#include "../commonTool/init_data.hpp"
#include "../commonTool/common_tool.hpp"
#include "inc/model.hpp"
#include "inc/method.hpp"
#include "inc/vector_visitor.hpp"
#include "inc/polygon.hpp"

#define D 3
int main(int argc, char* argv[])
{
  int k = atoi(argv[6]);
  int times = atoi(argv[5]);
  int q_num = atoi(argv[2]);
  int pw_num = atoi(argv[1]) * 1000;
  double q_range_low = atof(argv[3]);
  double q_range_up = atof(argv[4]);
  double tpm_time = 0, dtm_time = 0, nai_time = 0;



  for(int n = 0; n < times; n++){
      
      std::vector<std::vector<double> > p(pw_num);
      std::vector<std::vector<double> > w(pw_num);
      std::vector<std::vector<double> > qs(q_num);

      randinit(p,D,0,1);
      //read_file(p, D, "dianping.dat");
      //normal_init(p,D,0.5,1,0,1);
      //read_point_file(p,6,"data/houseD6.csv");
      //cluster_init(p,D,0,1,1);

      randinit_w(w,D,0,1);
      randinit(qs,D,q_range_low,q_range_up);
      //randinit(qt,D,0.1,0.2);

      std::vector<Point> points;
      std::vector<Weight> weights;

      for(int i = 0; i < p.size();i++){
          Point temp_p(i,p[i]);
          points.push_back(temp_p);
          //put_vector(p[i]);
      }
      for(int i = 0; i < w.size();i++){
          Weight temp_w(i,w[i]);
          weights.push_back(temp_w);
          //put_vector(w[i]);
      }



      /*cluster*/
      /*CATCH clucost;
      clucost.catch_time();
      int clu_number = std::pow(qs.size(),1/3.);
      std::cout<< "number" << clu_number << std::endl;
      auto c_Q = kmeans_Q(qs,D,clu_number);
      for(auto && c : c_Q){
          for(auto & cc : c){
              put_vector(cc);
          }
          std::cout << "next cluser." << std::endl;
      }
      clucost.catch_time();
      std::cout << "time: " << clucost.get_cost(2) << std::endl;*/

      /*Navie*/
      CATCH naviecost;
      naviecost.catch_time();
      auto navie_result = naive_arkr(points,weights,qs,k);
      naviecost.catch_time();

      /*init rtree_p,w*/
      rt rr;
      rr.creat_rtree_p_w(points,weights);

      /*TPM*/
      /*init rtree*/
      CATCH tpmcost;
      tpmcost.catch_time();
      auto tpm_result = tpm(rr.rtree, points, weights,qs,k,{},false);
      tpmcost.catch_time();


      /*DTM*/
      std::multimap<int,int> buffer;
      int dimension = qs[0].size();
      int current_rank = std::numeric_limits<int>::max();
      CATCH dtmcost;
      dtmcost.catch_time();
      namespace alo = boost::geometry::index::detail::rtree::utilities;
      /*CH*/
      dyy::poly::ConvexHull CH(qs);
      auto new_q = filter_Q_min_max(qs, CH.get_points(), 0, D);

      std::cout << "old_q: " << qs.size() << ",new_q: " << new_q.size() << std::endl;
      auto dtm_result =
          alo::vector_visitor(rr.rtree_w,new_q,rr.rtree,current_rank,k);
      dtmcost.catch_time();


      dtm_time += dtmcost.get_cost(2);
      tpm_time += tpmcost.get_cost(2);
      nai_time += naviecost.get_cost(2);

      /*report*/
      /*std::cout << "naive: " << std::endl;
      print_map(navie_result);
      std::cout << "TPM: " << std::endl;
      print_map(tpm_result);
      std::cout << "DTM: " << std::endl;
      print_map(dtm_result);*/

  } //times loop

  /*average time report*/
  std::cout << "naiv: cpu cost is " << (double) nai_time / times << " millisecond(s)" << std::endl;
  std::cout << "DTM: cpu cost is " << (double) tpm_time / times << " millisecond(s)" << std::endl;
  std::cout << "CDTM: cpu cost is " << (double) dtm_time / times << " millisecond(s)" << std::endl;
}
