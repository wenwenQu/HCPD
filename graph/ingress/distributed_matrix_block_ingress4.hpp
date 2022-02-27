/**  
 * Copyright (c) 2009 Carnegie Mellon University.
 *     All rights reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing,
 *  software distributed under the License is distributed on an "AS
 *  IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
 *  express or implied.  See the License for the specific language
 *  governing permissions and limitations under the License.
 *
 * For more about this software visit:
 *
 *      http://www.graphlab.ml.cmu.edu
 *
 */

#ifndef GRAPHLAB_DISTRIBUTED_MATRIX_RR_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_MATRIX_RR_INGRESS_HPP

#include <boost/functional/hash.hpp>

#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/ingress/distributed_ingress_base.hpp>
#include <graphlab/graph/distributed_graph.hpp>

#include <set>


#include <graphlab/macros_def.hpp>
namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class distributed_graph;

  /**
   * \brief Ingress object assigning edges using randoming hash function.
   */
  template<typename VertexData, typename EdgeData>
  class distributed_matrix_block_ingress :
    public distributed_ingress_base<VertexData, EdgeData> {
  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData   edge_data_type;

    typedef distributed_ingress_base<VertexData, EdgeData> base_type;

    typedef typename graph_type::vertex_record vertex_record;
    typedef typename graph_type::mirror_type mirror_type;

    typedef typename buffered_exchange<vertex_id_type>::buffer_type
        vertex_id_buffer_type;

    /// The rpc interface for this object
    dc_dist_object<distributed_matrix_block_ingress> matrix_rpc;
    /// The underlying distributed graph object that is being loaded
    graph_type& graph;

    /// threshold to divide high-degree and low-degree vertices
    size_t threshold;

    //@topoX
      size_t theta;
      size_t etheta;
      size_t ceng;

    bool standalone;

    typedef typename base_type::edge_buffer_record edge_buffer_record;
    typedef typename buffered_exchange<edge_buffer_record>::buffer_type
        edge_buffer_type;

    typedef typename base_type::vertex_buffer_record vertex_buffer_record;
    typedef typename buffered_exchange<vertex_buffer_record>::buffer_type
        vertex_buffer_type;

    std::vector<edge_buffer_record> matrix_edges;

    //@hybrid
    hopscotch_map<vertex_id_type, size_t> in_degree_set;

    /* ingress exchange */
    buffered_exchange<edge_buffer_record> matrix_edge_exchange;
    buffered_exchange<vertex_buffer_record> matrix_vertex_exchange;




    /// detail vertex record for the second pass coordination.

      typedef typename base_type::vertex_negotiator_record
              vertex_negotiator_record;

      //@wenwen: topox codes
      typedef typename boost::unordered_map<vertex_id_type, std::vector<vertex_id_type> > target_hash_table_type;
      std::queue<vertex_id_type> p_queue;

      target_hash_table_type iht;
      target_hash_table_type oht;
      typedef typename boost::unordered_map<vertex_id_type,int> flag_hash_table_type;
      flag_hash_table_type fht;
      flag_hash_table_type flaght;
      flag_hash_table_type countht;
      typedef typename boost::unordered_map<vertex_id_type, size_t> degree_hash_table_type;
      degree_hash_table_type dht;
      mutex point_lock;
      //mutex lock2;

      typedef typename boost::unordered_map<vertex_id_type, size_t > vertex_hash_table_type;
      vertex_hash_table_type vht;
      typedef typename boost::unordered_map<vertex_id_type, std::vector<size_t> > vertex_hash_table_local_type;
      vertex_hash_table_local_type vht_local;
      size_t* count;         //璁板綍鍚勮妭鐐圭幇鏈夌殑鍒嗛厤鍒扮殑杈圭殑鏁扮洰锛屽浣曞垎甯冨紡锛�
      typedef typename std::pair<vertex_id_type, size_t> hash_pair_type;
      buffered_exchange<hash_pair_type> vht_exchange;
      typedef typename buffered_exchange<hash_pair_type>::buffer_type vht_buffer_type;
   
  public:

      //@matrix
    std::vector<std::vector<int> > myMatrix;
	int matrixSize, maxValue;

	//@ginder
    /** Type of the master location hash table: [vertex-id, location-of-master] */
    typedef typename boost::unordered_map<vertex_id_type, procid_t>
        master_hash_table_type;
    typedef typename std::pair<vertex_id_type, procid_t>
        master_pair_type;
    typedef typename buffered_exchange<master_pair_type>::buffer_type
        master_buffer_type;
    typedef typename buffered_exchange<hopscotch_map<vertex_id_type, size_t> >::buffer_type degree_buffer_type;
    /** master hash table (mht): location mapping of low-degree vertices  */
    master_hash_table_type mht;
    buffered_exchange<master_pair_type> mht_exchange;

    buffered_exchange<hopscotch_map<vertex_id_type, size_t> > in_degree_exchange;


      distributed_matrix_block_ingress(distributed_control& dc,
      graph_type& graph, size_t threshold = 100, size_t theta = 100,size_t etheta=1000,size_t ceng=1) :
      base_type(dc, graph), matrix_rpc(dc, this),
      graph(graph), threshold(threshold), theta(theta), etheta(etheta),ceng(ceng),
#ifdef _OPENMP
      matrix_edge_exchange(dc, omp_get_max_threads()),
		vht_exchange(dc, omp_get_max_threads()),
        matrix_vertex_exchange(dc, omp_get_max_threads())
#else
      matrix_edge_exchange(dc),
      vht_exchange(dc),
      matrix_vertex_exchange(dc)
#endif
		, mht_exchange(dc)
		,in_degree_exchange(dc)
	{
    	std::fstream fin(graph.partitionfile.c_str());
    	fin >> matrixSize >> maxValue;
    	myMatrix.resize(matrixSize);
    	for(int i = 0;i<matrixSize;i++){
    		myMatrix[i].resize(matrixSize);
    		for(int j = 0; j< matrixSize; j++){
    			fin >> myMatrix[i][j];
			}
    	}

        count = (size_t *) malloc((matrix_rpc.numprocs() / ceng)*sizeof (size_t));
    } // end of constructor

    ~distributed_matrix_block_ingress() {
          free (count);
      }


      /** Add an edge to the ingress object using random hashing assignment.
       *  This function acts as the first phase for SNAP graph to deliver edges
       *  via the hashing value of its target vertex.
       */
      void add_edge(vertex_id_type source, vertex_id_type target,
                    const EdgeData& edata) {
          size_t numproc = matrix_rpc.numprocs();
          const edge_buffer_record record(source, target, edata);
          const size_t sh = source%numproc;
          const size_t th = target%numproc;
#ifdef _OPENMP
          if (sh!=th)
			matrix_edge_exchange.send(sh, record, omp_get_thread_num());
		matrix_edge_exchange.send(th, record, omp_get_thread_num());
#else
          if (sh != th)
              matrix_edge_exchange.send(sh, record);                  //锟斤拷一锟斤拷锟竭凤拷锟斤拷source锟斤拷target锟斤拷锟节碉拷锟斤拷锟斤拷锟节碉拷锟铰ｏ拷锟斤拷锟轿伙拷锟酵伙拷诘悖伙拷锟斤拷锟揭伙拷锟�
          matrix_edge_exchange.send(th, record);
#endif

      } // end of add edge



      /* add vdata */
      void add_vertex(vertex_id_type vid, const VertexData& vdata) {
          const vertex_buffer_record record(vid, vdata);
          const procid_t owning_proc = standalone ? 0 :
                                       vid % matrix_rpc.numprocs();
#ifdef _OPENMP
          matrix_vertex_exchange.send(owning_proc, record, omp_get_thread_num());
#else
          matrix_vertex_exchange.send(owning_proc, record);
#endif
      } // end of add vertex

      void set_queue(const size_t& v, const size_t& owning,const size_t& vp,const size_t& matrixcount){
          if (flaght[v] == 0){
              point_lock.lock();
              p_queue.push(v);
              if (fht[v] == -1){
                  fht[v] = owning;
                  countht[v] = matrixcount+1;
              }
              point_lock.unlock();
          }
          if (flaght[vp] == 0){
              flaght[vp] = 1;
          }
          if (fht[vp] == -1){
              fht[vp] = owning;
          }

      }
      void add_queue(const size_t& v,const size_t& owning,const size_t& vp,const size_t& matrixcount){
          size_t owningmachine = v%matrix_rpc.numprocs();
          if (owningmachine != matrix_rpc.procid())
              matrix_rpc.remote_call(owningmachine, &distributed_matrix_block_ingress<VertexData, EdgeData>::set_queue, v, owning, vp,matrixcount);
          else
              set_queue(v, owning,vp,matrixcount);
      }


      void finalize() {

          graphlab::timer ti;
          int ceng = 1;
          size_t nprocs = matrix_rpc.numprocs();
          size_t cnodes = nprocs / ceng;
          procid_t l_procid = matrix_rpc.procid();
          size_t nedges = 0;

          matrix_rpc.full_barrier();

          if (l_procid == 0) {
              memory_info::log_usage("start finalizing");
              logstream(LOG_EMPH) << "matrix finalizing ..."
                                  << " #vertices=" << graph.local_graph.num_vertices()
                                  << " #edges=" << graph.local_graph.num_edges()
                                  << " threshold=" << threshold
                                  << " theta=" << theta
                                  << "etheta=" <<etheta
                                  << std::endl;
          }


          /**************************************************************************/
          /*                                                                        */
          /*                       Flush any additional data                        */
          /*                                                                        */
          /**************************************************************************/
          matrix_edge_exchange.flush(); matrix_vertex_exchange.flush();

          /**
           * Fast pass for redundant finalization with no graph changes.
           */
          {
              size_t changed_size = matrix_edge_exchange.size() + matrix_vertex_exchange.size();
              matrix_rpc.all_reduce(changed_size);
              if (changed_size == 0) {
                  logstream(LOG_INFO) << "Skipping Graph Finalization because no changes happened..." << std::endl;
                  return;
              }
          }


          /**************************************************************************/
          /*                                                                        */
          /*                       Prepare matrix ingress                           */
          /*                                                                        */
          /**************************************************************************/
          {
              edge_buffer_type edge_buffer;
              procid_t proc;
              nedges = matrix_edge_exchange.size();

              matrix_edges.reserve(nedges);
              if (standalone) { /* fast pass for standalone */
                  proc = -1;
                  while(matrix_edge_exchange.recv(proc, edge_buffer))
                      foreach(const edge_buffer_record& rec, edge_buffer){
                      if (rec.target%nprocs == l_procid){
                          matrix_edges.push_back(rec);
                          iht[rec.target].push_back(rec.source); //in-degree hash table
                          if (rec.source%nprocs == l_procid)
                              oht[rec.source].push_back(rec.target); //out-degree hash table
                      }
                      else{
                          oht[rec.source].push_back(rec.target);
                      }
                      fht[rec.source] = -1;
                      fht[rec.target] = -1;
                      flaght[rec.target] = 0;
                      flaght[rec.source] = 0;

                  }

                  matrix_edge_exchange.clear();
              } else {


                  proc = -1;
                  while(matrix_edge_exchange.recv(proc, edge_buffer)) {

                      foreach(const edge_buffer_record& rec, edge_buffer){         //鎸夌洰鐨勭偣瀛樺偍璇ヨ竟锛屽苟鏍规嵁杈规瀯寤篿ht鍜宱ht
                          if (rec.target%nprocs == l_procid){
                              matrix_edges.push_back(rec);
                              iht[rec.target].push_back(rec.source);
                              if (rec.source%nprocs == l_procid)
                                  oht[rec.source].push_back(rec.target);
                          }
                          else{
                              oht[rec.source].push_back(rec.target);
                          }
                          fht[rec.source] = -1;                                 //fht浣滅敤锛歮apping proc hash table
                          fht[rec.target] = -1;
                          flaght[rec.target] = 0;
                          flaght[rec.source] = 0;
                          countht[rec.target] = 0;

                      }


                  }
                  matrix_edge_exchange.clear();
                  matrix_edge_exchange.barrier(); // barrier before reusing
#ifdef TUNING
                  if(l_procid == 0) {
            logstream(LOG_INFO) << "iet,oet building in parallel: "
                                << ti.current_time()
                                << " secs"
                                << std::endl;
          }
#endif


                  for (target_hash_table_type::iterator iter = iht.begin(); iter != iht.end(); iter++){
                      in_degree_set[iter->first] = iht[iter->first].size();
                  }


                  for (procid_t p = 0; p < nprocs; ++p) {
                      if (p != l_procid)
                          in_degree_exchange.send(p, in_degree_set);
                  }

                  // synchronize on in_degree_set
                  in_degree_exchange.flush();
                  degree_buffer_type degree_buffer;
                  proc = -1;

                  typedef hopscotch_map<vertex_id_type, size_t> tesf;
                  while(in_degree_exchange.recv(proc, degree_buffer)) {
                      foreach(tesf& ot, degree_buffer)
                      for(hopscotch_map<vertex_id_type, size_t>::iterator it = ot.begin(); it!=ot.end(); it++){
                          in_degree_set[it->first] = it->second;
                      }
                  }
                  in_degree_exchange.clear();

                  //size_t point[1000]; // 褰撳墠绛夊緟璁＄畻鐨勭偣鐨勯槦鍒�
                  //size_t topoX_count=0;
                  //int point_count = -1;
                    /*
                  procid_t min_proc = 0;
                  if (fht[iter->first] == -1){

                      int temp2 = count[0];
                      for (size_t ix = 1; ix < cnodes; ix++){
                          if (temp2>count[ix]){
                              temp2 = count[ix];
                              min_proc = ix;
                          }
                      }
                      fht[iter->first] = min_proc;                  //绗竴娆℃洿鏂癴ht锛屾鏃跺彧鏇存柊鏈湴fht
                  }
                     */

                  for (size_t ix = 0; ix < cnodes; ix++)
                      count[ix] = 0;

                  //boost::unordered_map<vertex_id_type, std::vector<size_t> > high_node_neighbor;

                  for (target_hash_table_type::iterator iter = iht.begin(); iter != iht.end(); iter++){
                	  if (in_degree_set[iter->first] > threshold){
                		  /*
                		  int cnt = 0;
                		  for (size_t j = 0; j < iht[iter->first].size(); j++){
                			  if(in_degree_set[iht[iter->first][j]] <= threshold){
                				  //high_node_neighbor[iter->first].push_back(iht[iter->first][j]);
                				  cnt+=in_degree_set[iht[iter->first][j]];
                			  }
                		  }
                		  for (size_t j = 0; j < oht[iter->first].size(); j++){
                			  if(in_degree_set[oht[iter->first][j]] <= threshold){
                			      //high_node_neighbor[iter->first].push_back(oht[iter->first][j]);
                				  cnt+=in_degree_set[oht[iter->first][j]];
                			  }
                		  }
                		  */

                		  //todo the proc_id in current col can be re
                		  int col = iter->first % matrixSize;
                		  int temp2 = count[myMatrix[0][col]];
                		  size_t min_proc = myMatrix[0][col];
                		  for(int ix=1;ix<matrixSize;ix++){
                			  if (temp2>count[myMatrix[ix][col]]){
                				  temp2 = count[myMatrix[ix][col]];
                			      min_proc = myMatrix[ix][col];
                			  }
                		  }

                		  count[min_proc] += iht[iter->first].size();

                		  flaght[iter->first] = 1;
                		  fht[iter->first] = min_proc;

                		  for (procid_t p = 0; p < nprocs; ++p) {
							  if (p != l_procid)
								  mht_exchange.send(p, master_pair_type(iter->first, min_proc));
							  else
								  mht[iter->first] = min_proc;
                		  }

                          for (size_t j = 0; j < iht[iter->first].size(); j++){
                              if (fht[iht[iter->first][j]] == -1 && in_degree_set[iht[iter->first][j]] <=threshold){
                            	  //vht_local[iht[rec.target][j]].push_back(rec.target);
                                  fht[iht[iter->first][j]] = min_proc;
                                  add_queue(iht[iter->first][j], min_proc, iter->first, -1);
                              }
                          }
                          for (size_t j = 0; j < oht[iter->first].size(); j++){
                              if (fht[oht[iter->first][j]] == -1 && in_degree_set[oht[iter->first][j]] <=threshold){
                            	  //vht_local[oht[rec.target][j]].push_back(rec.target);
                                  fht[oht[iter->first][j]] = min_proc;
                                  add_queue(oht[iter->first][j], min_proc, iter->first, -1);
                              }
                          }
                	  }
                  }






                  //float avg_cnt = 0;
                  /*
                  for (size_t i = 0; i < matrix_edges.size(); i++){
                      procid_t owning_proc = 0;
                      edge_buffer_record& rec = matrix_edges[i];


                      if (in_degree_set[rec.source] > threshold && in_degree_set[rec.target] > threshold){

                          int row, col;
                          col = (rec.source) % matrixSize;
                          row = (rec.target) % matrixSize;
                          if (col == row) col = rand() % matrixSize;


                          owning_proc = myMatrix[row][col];
                       */
                          /*
                          if (mht.find(rec.target) == mht.end()) {
                              // update mht
                              for (procid_t p = 0; p < nprocs; ++p) {
                                 if (p != l_procid)
                                     mht_exchange.send(p, master_pair_type(rec.target, owning_proc));
                                 else
                                     mht[rec.target] = owning_proc;
                              }

                              //count[owning_proc] += dht[rec.target];
                              //avg_cnt = (avg_cnt*nprocs+dht[rec.target])/dht[rec.target];

                              // update vht of high-degree's neighbor(must be low degree)
                              count[owning_proc] += in_degree_set[rec.target];
                              for (size_t j = 0; j < iht[rec.target].size(); j++){
                                  if (fht[iht[rec.target][j]] == -1 && in_degree_set[iht[rec.target][j]] <=threshold){
                                	  //vht_local[iht[rec.target][j]].push_back(rec.target);
                                      fht[iht[rec.target][j]] = owning_proc;
                                      add_queue(iht[rec.target][j], owning_proc, rec.target, -1);
                                  }
                              }
                              for (size_t j = 0; j < oht[rec.target].size(); j++){
                                  if (fht[oht[rec.target][j]] == -1 && in_degree_set[oht[rec.target][j]] <=threshold){
                                	  //vht_local[oht[rec.target][j]].push_back(rec.target);
                                      fht[oht[rec.target][j]] = owning_proc;
                                      add_queue(oht[rec.target][j], owning_proc, rec.target, -1);
                                  }
                              }

			  	  	  	  }
			  	  	  	  */
			  /*
			  if (fht[rec.source] == -1 && in_degree_set[rec.source] <=threshold){
				  fht[rec.source] = owning_proc;
				  add_queue(rec.source, owning_proc, rec.target, -1);
			  }
*/
                      	  /*
                          if (owning_proc != l_procid){

                              matrix_edge_exchange.send(owning_proc%cnodes, rec);


                              // set re-sent edges as empty for skipping
                              matrix_edges[i] = edge_buffer_record();
                              --nedges;
                          }

                      }
                  }
*/
                  /*
                  {
                	  //procid_t min_neighbor;
			  	  	  int cnt = 0;
                	  for (vertex_hash_table_local_type::iterator iter = vht_local.begin(); iter != vht_local.end(); iter++){
                		  std::vector<size_t>& vec = iter->second;
                    	  int temp2 = count[mht[vec[0]]];
				  	  	  size_t min_neighbor = vec[0];
				  	  	  assert(mht.find(vec[0])!=mht.end());
                		  for (int i = 1; i< vec.size(); i++){
					  	  	  cnt++;
					  	  	  //std::cout << vec[i] << " " ;
					  	  	  size_t cur = vec[i];
					  	  	  assert(mht.find(cur)!=mht.end());
                			  if (temp2>count[mht[cur]]){
                				  temp2 = count[mht[cur]];
                				  min_neighbor = cur;
                		      }
				  	  	  }

				  	  	  //std::cout << ">>"<< vec.size() << "|" << min_neighbor << std::endl;
				  	  	  assert(mht.find(min_neighbor)!=mht.end());
				  	  	  count[mht[min_neighbor]]++;
                		  fht[iter->first] = mht[min_neighbor];
                		  add_queue(iter->first, mht[min_neighbor], min_neighbor, -1);
                	  }
			  	  	  std::cout << "@@@@@@@@@@@@" << cnt << std::endl;
                  }
*/
                  for (size_t ix = 0; ix < cnodes; ix++)
                       std::cout << count[ix] << " " 	;
                  std::cout << std::endl;


                  mht_exchange.flush();
                  master_buffer_type mht_buffer;
                  proc = -1;
                  while (mht_exchange.recv(proc, mht_buffer)){
                      foreach(const master_pair_type& pair, mht_buffer){
                          mht[pair.first] = pair.second;
                          //logstream(LOG_INFO) << "first = " << pair.first << " second = " << pair.second << std::endl;
                      }
                  }
                  mht_exchange.clear();





                  double last_time = 0.0+ti.current_time();
                  //for (size_t ix = 0; ix < cnodes; ix++)
                      //count[ix] = 0;
                  for (target_hash_table_type::iterator iter = oht.begin(); iter != oht.end(); iter++){
                      if(iht.find(iter->first) == iht.end()){
                          procid_t min_proc = 0;
                          iht[iter->first];
                          /*
                          if (fht[iter->first] == -1){
                              int temp2 = count[0];
                              for (size_t ix = 1; ix < cnodes; ix++){
                                  if (temp2>count[ix]){
                                      temp2 = count[ix];
                                      min_proc = ix;
                                  }
                              }
                              fht[iter->first] = min_proc;                  //绗竴娆℃洿鏂癴ht锛屾鏃跺彧鏇存柊鏈湴fht
                          }

                          for (size_t ix = 0; ix < nprocs; ix++){
                              if (ix != l_procid)
                                  vht_exchange.send(ix, hash_pair_type(iter->first, min_proc ));
                              else
                                  vht[iter->first] = min_proc;
                          }
                           */
                      }
                  }

                  for (target_hash_table_type::iterator iter = iht.begin(); iter != iht.end(); iter++){
                      //logstream(LOG_INFO) << "in ok!" << std::endl;
                      //   topoX_rpc.cout() << iht[iter->first].size() << std::endl;
                      if (l_procid == 0&&ti.current_time()-last_time>2) {
                          for (size_t ix = 0; ix < cnodes; ix++)
                              logstream(LOG_INFO) << count[ix] <<" ";
                          logstream(LOG_INFO) << std::endl;
                          last_time = 0.0+ti.current_time();
                      }
                      dht[iter->first] = iht[iter->first].size();//@wenwen: can be replace by in_degree_set
                      if (!p_queue.empty()){
                          //logstream(LOG_INFO) << "queue ok!" << std::endl;
                          while (!p_queue.empty()){
                              point_lock.lock();
                              size_t p_point = p_queue.front(); //@wenwen: p_point maybe high-degree
                              p_queue.pop();
                              point_lock.unlock();
                              if (flaght[p_point] == 0){ // @wenwen锛� flag hash table

                                  flaght[p_point] += 1; // @wenwen: check whether is visited

                                  procid_t min_proc = fht[p_point]; //@wenwen: current -1

                                  if (countht[p_point] > etheta){ //@wenwen: move to the smallest proc
                                      int temp2 = count[0];
                                      min_proc = 0;
                                      for (size_t ix = 1; ix < cnodes; ix++){
                                          if (temp2>count[ix]){
                                              temp2 = count[ix];
                                              min_proc = ix;
                                          }
                                      }
                                      fht[p_point] = min_proc;
                                      countht[p_point] = 0;
                                  }
                                  for (size_t ix = 0; ix < nprocs; ix++){
                                      if (ix != l_procid)
                                          vht_exchange.send(ix, hash_pair_type(p_point, min_proc));
                                      else
                                          vht[p_point] = min_proc; //@wenwen: vertex hash table
                                  }
                                  dht[p_point] = iht[p_point].size();
                                  if (dht[p_point] <= threshold)
                                      count[min_proc] += dht[p_point]; //@wenwen: we want to place p_point to min_proc
                                  if (dht[p_point] <= threshold&&oht[p_point].size() <= theta){

                                      for (size_t j = 0; j < dht[p_point]; j++){
                                          //if (fht[iht[p_point][j]] == -1 && in_degree_set[iht[p_point][j]] <= threshold){
                                    	  if (fht[iht[p_point][j]] == -1){
                                              fht[iht[p_point][j]] = min_proc;
                                              add_queue(iht[p_point][j], min_proc, p_point, countht[p_point]);
                                          }
                                      }
                                      for (size_t j = 0; j < oht[p_point].size(); j++){
                                          //if (fht[oht[p_point][j]] == -1 && in_degree_set[oht[p_point][j]] <= threshold){
                                    	  if (fht[oht[p_point][j]] == -1){
                                              fht[oht[p_point][j]] = min_proc;
                                              add_queue(oht[p_point][j], min_proc, p_point, countht[p_point]);
                                          }
                                      }
                                  }



                              }

                          }
                          if (flaght[iter->first] == 0 && dht[iter->first]<=threshold){

                              flaght[iter->first]+= 1;

                              procid_t min_proc = 0;
                              if (fht[iter->first] == -1){

                                  int temp2 = count[0];
                                  for (size_t ix = 1; ix < cnodes; ix++){
                                      if (temp2>count[ix]){
                                          temp2 = count[ix];
                                          min_proc = ix;
                                      }
                                  }
                                  fht[iter->first] = min_proc;                  //绗竴娆℃洿鏂癴ht锛屾鏃跺彧鏇存柊鏈湴fht
                              }
                              else {
                                  if (countht[iter->first] > etheta){
                                      int temp2 = count[0];
                                      for (size_t ix = 1; ix < cnodes; ix++){
                                          if (temp2>count[ix]){
                                              temp2 = count[ix];
                                              min_proc = ix;
                                          }
                                      }
                                      fht[iter->first] = min_proc;
                                      countht[iter->first] = 0;
                                  }
                                  else min_proc = fht[iter->first];
                              }
                              for (size_t ix = 0; ix < nprocs; ix++){
                                  if (ix != l_procid)
                                      vht_exchange.send(ix, hash_pair_type(iter->first, min_proc));
                                  else
                                      vht[iter->first] = min_proc;
                              }
                              if (dht[iter->first] <= threshold)
                                  count[min_proc] += dht[iter->first];
                              if (dht[iter->first] <= threshold&&oht[iter->first].size() <= theta){

                                  for (size_t j = 0; j < dht[iter->first]; j++){
                                      //if (fht[iht[iter->first][j]] == -1 && in_degree_set[iht[iter->first][j]] <= threshold){
                                	  if (fht[iht[iter->first][j]] == -1){
                                          fht[iht[iter->first][j]] = min_proc;
                                          add_queue(iht[iter->first][j], min_proc, iter->first, countht[iter->first]);
                                      }
                                  }

                                  for (size_t j = 0; j < oht[iter->first].size(); j++){
                                      //if (fht[oht[iter->first][j]] == -1 && in_degree_set[oht[iter->first][j]] <= threshold){
                                	  if (fht[oht[iter->first][j]] == -1){
                                          fht[oht[iter->first][j]] = min_proc;
                                          add_queue(oht[iter->first][j], min_proc, iter->first, countht[iter->first]);
                                      }
                                  }

                              }
                          }


                      }
                      else{
                          if (flaght[iter->first] == 0 && dht[iter->first] <= threshold){

                              flaght[iter->first] += 1;

                              procid_t min_proc = 0;
                              if (fht[iter->first] == -1){

                                  int temp2 = count[0];
                                  for (size_t ix = 1; ix < cnodes; ix++){
                                      if (temp2>count[ix]){
                                          temp2 = count[ix];
                                          min_proc = ix;
                                      }
                                  }
                                  fht[iter->first] = min_proc;                  //绗竴娆℃洿鏂癴ht锛屾鏃跺彧鏇存柊鏈湴fht
                              }
                              else {
                                  if (countht[iter->first] > etheta){
                                      int temp2 = count[0];
                                      for (size_t ix = 1; ix < cnodes; ix++){
                                          if (temp2>count[ix]){
                                              temp2 = count[ix];
                                              min_proc = ix;
                                          }
                                      }
                                      fht[iter->first] = min_proc;
                                      countht[iter->first] = 0;
                                  }
                                  else min_proc = fht[iter->first];
                              }
                              for (size_t ix = 0; ix < nprocs; ix++){
                                  if (ix != l_procid)
                                      vht_exchange.send(ix, hash_pair_type(iter->first, min_proc ));
                                  else
                                      vht[iter->first] = min_proc;
                              }
                              if (dht[iter->first] <= threshold)
                                  count[min_proc] += dht[iter->first];
                              if (dht[iter->first] <= threshold&&oht[iter->first].size() <= theta){

                                  for (size_t j = 0; j < dht[iter->first]; j++){
                                      //if (fht[iht[iter->first][j]] == -1 && in_degree_set[iht[iter->first][j]] <= threshold){
                                	  if (fht[iht[iter->first][j]] == -1){
                                          fht[iht[iter->first][j]] = min_proc;
                                          add_queue(iht[iter->first][j], min_proc, iter->first, countht[iter->first]);
                                      }
                                  }

                                  for (size_t j = 0; j < oht[iter->first].size(); j++){
                                      //if (fht[oht[iter->first][j]] == -1 && in_degree_set[oht[iter->first][j]] <= threshold){
                                	  if (fht[oht[iter->first][j]] == -1){
                                          fht[oht[iter->first][j]] = min_proc;
                                          add_queue(oht[iter->first][j], min_proc, iter->first, countht[iter->first]);
                                      }
                                  }

                              }
                          }
                      }







                  }

                  for (target_hash_table_type::iterator iter = iht.begin(); iter != iht.end(); iter++){
                      if (flaght[iter->first] > 1){
                          logstream(LOG_INFO) << "flaght>2: " << iter->first << std::endl;
                      }

                  }
                  matrix_rpc.full_barrier();
                  if (l_procid == 0) {
                      logstream(LOG_INFO) << "matrix time: "
                                          << ti.current_time()
                                          << " secs"
                                          << std::endl;

                  }
                  iht.clear();
                  oht.clear();
                  //fht.clear();
                  // flaght.clear();
                  //countht.clear();
                  vht_exchange.flush();
                  vht_buffer_type vht_buffer;
                  proc = -1;
                  while (vht_exchange.recv(proc, vht_buffer)){
                      foreach(const hash_pair_type& pair, vht_buffer){
                          vht[pair.first] = pair.second;
                          //logstream(LOG_INFO) << "first = " << pair.first << " second = " << pair.second << std::endl;
                      }
                  }
                  vht_exchange.clear();
                  //fht.clear();
                  //vht_exchange.flush();                                   //
                  //vht_buffer_type vht_buffer;
                  //proc = -1;
                  // while (vht_exchange.recv(proc, vht_buffer)){
                  //	  foreach(const hash_pair_type& pair, vht_buffer){
                  //	  vht[pair.first] = pair.second;
                  //logstream(LOG_INFO) << "first = " << pair.first << " second = " << pair.second << std::endl;
                  //	  }
                  // }
                  //  vht_exchange.clear();
                  if (l_procid == 0) {
                      for (size_t ix = 0; ix < cnodes; ix++)
                          logstream(LOG_INFO) << " count[ " << ix << " ] = " << count[ix] << std::endl;
                  }
                  // re-send edges of high-degree vertices


                  for (size_t i = 0; i < matrix_edges.size(); i++){
                      procid_t owning_proc = 0;
                      edge_buffer_record& rec = matrix_edges[i];
                      if (dht[rec.target] <= threshold){
                          if (vht.find(rec.target) != vht.end() && vht[rec.target]!=-1){
                              owning_proc = vht[rec.target];
                          }
                          else{
                              owning_proc = rec.target%cnodes;
                          }
                          if (owning_proc != l_procid){
                              for (int j = 0; j < ceng; j++){
                                  matrix_edge_exchange.send(owning_proc%cnodes, rec);
                              }

                              // set re-sent edges as empty for skipping
                              matrix_edges[i] = edge_buffer_record();
                              --nedges;
                          }

                      }
                      else{
                    	  /*
                    	  if (vht.find(rec.source) != vht.end()&&vht[rec.source]!=-1){
                    	  	  owning_proc = vht[rec.source];
                    	  }*/

                    	  if (vht.find(rec.source) != vht.end()&&vht[rec.source]!=-1&&mht.find(rec.target)!=mht.end()){
                    		  owning_proc = mht[rec.target];
                    	  }

                          /*
                          else if(in_degree_set[rec.source] <= threshold){
                              owning_proc = rec.source%cnodes;
                          }
                           */
                          else {
                              int row, col;
                              col = (rec.source) % matrixSize;
                              //row = (rec.target) % matrixSize;
                              row = (rec.target) % matrixSize;
                              if (col == row) col = rand() % matrixSize;


                              owning_proc = myMatrix[row][col];
                          }



                          if (owning_proc != l_procid){

                              matrix_edge_exchange.send(owning_proc%cnodes, rec);


                              // set re-sent edges as empty for skipping
                              matrix_edges[i] = edge_buffer_record();
                              --nedges;
                          }
                      }
                  }


#ifdef TUNING
                  if(l_procid == 0) {
            logstream(LOG_INFO) << "resend edges time: "
                                << ti.current_time()
                                << " secs"
                                << std::endl;

          }
#endif

                  // receive edges of high-degree vertices
                  matrix_edge_exchange.flush();
#ifdef TUNING
                  logstream(LOG_INFO) << "receive high-degree edges: "
                              << matrix_edge_exchange.size() << std::endl;
#endif
                  proc = -1;
                  while(matrix_edge_exchange.recv(proc, edge_buffer)) {
                      foreach(const edge_buffer_record& rec, edge_buffer) {
                          matrix_edges.push_back(rec);
                          ++nedges;
                      }
                  }
                  matrix_edge_exchange.clear();

#ifdef TUNING
                  if(l_procid == 0) {
            logstream(LOG_INFO) << "receive resend edges time: "
                                << ti.current_time()
                                << " secs"
                                << std::endl;
          }
#endif
              }
          }

          if(l_procid == 0) {
              memory_info::log_usage("prepare matrix finalizing done.");
              logstream(LOG_EMPH) << "prepare matrix finalizing. ("
                                  << ti.current_time()
                                  << " secs)"
                                  << std::endl;
          }

          // connect to base finalize()
          modified_base_finalize(nedges);

          // set vertex degree type for topoX engine
          set_degree_type();

          if(l_procid == 0) {
              memory_info::log_usage("matrix finalizing graph done.");
              logstream(LOG_EMPH) << "matrix finalizing graph. ("
                                  << ti.current_time()
                                  << " secs)"
                                  << std::endl;
          }
      } // end of finalize

    void set_degree_type() {
      graphlab::timer ti;
      procid_t l_procid = matrix_rpc.procid();
      size_t high_master = 0, high_mirror = 0, low_master = 0, low_mirror = 0;

      for (size_t lvid = 0; lvid < graph.num_local_vertices(); lvid++) {
        vertex_record& vrec = graph.lvid2record[lvid];
        if (vrec.num_in_edges > threshold) {
          vrec.dtype = graph_type::HIGH;
          if (vrec.owner == l_procid) high_master ++;
          else high_mirror ++;
        } else {
          vrec.dtype = graph_type::LOW;
          if (vrec.owner == l_procid) low_master ++;
          else low_mirror ++;
        }
      }
      std::cout << "proc " << l_procid << ": " << high_master+high_mirror <<"|" <<low_master+low_mirror << std::endl;// << << << <<

#ifdef TUNING
      // Compute the total number of high-degree and low-degree vertices
      std::vector<size_t> swap_counts(matrix_rpc.numprocs());

      swap_counts[l_procid] = high_master;
      matrix_rpc.all_gather(swap_counts);
      high_master = 0;
      foreach(size_t count, swap_counts) high_master += count;

      swap_counts[l_procid] = high_mirror;
      matrix_rpc.all_gather(swap_counts);
      high_mirror = 0;
      foreach(size_t count, swap_counts) high_mirror += count;

      swap_counts[l_procid] = low_master;
      matrix_rpc.all_gather(swap_counts);
      low_master = 0;
      foreach(size_t count, swap_counts) low_master += count;

      swap_counts[l_procid] = low_mirror;
      matrix_rpc.all_gather(swap_counts);
      low_mirror = 0;
      foreach(size_t count, swap_counts) low_mirror += count;

      if(l_procid == 0) {
        logstream(LOG_EMPH) << "matrix info: master ["
                            << high_master << " "
                            << low_master << " "
                            << (float(high_master)/(high_master+low_master)) << "]"
                            << std::endl;
        if ((high_mirror + low_mirror) > 0)
        logstream(LOG_EMPH) << "matrix info: mirror ["
                            << high_mirror << " "
                            << low_mirror << " "
                            << (float(high_mirror)/(high_mirror+low_mirror)) << "]"
                            << std::endl;

        memory_info::log_usage("set vertex type done.");
        logstream(LOG_EMPH) << "set vertex type: "
                            << ti.current_time()
                            << " secs"
                            << std::endl;
      }
#endif
    }


    /*
     * do the same job as original base finalize except for
     * extracting edges from matrix_edges instead of original edge_buffer
     */
      /*
       * do the same job as original base finalize except for
       * extracting edges from matrix_edges instead of original edge_buffer
       */
      void modified_base_finalize(size_t nedges) {
          graphlab::timer ti;
          procid_t l_procid = matrix_rpc.procid();
          size_t nprocs = matrix_rpc.numprocs();
          size_t cnodes = nprocs / ceng;

          matrix_rpc.full_barrier();

          bool first_time_finalize = false;
          /**
           * Fast pass for first time finalization.
           */
          if (graph.is_dynamic()) {
              size_t nverts = graph.num_local_vertices();
              matrix_rpc.all_reduce(nverts);
              first_time_finalize = (nverts == 0);
          } else {
              first_time_finalize = false;
          }


          typedef typename hopscotch_map<vertex_id_type, lvid_type>::value_type
                  vid2lvid_pair_type;

          /**
           * \internal
           * Buffer storage for new vertices to the local graph.
           */
          typedef typename graph_type::hopscotch_map_type vid2lvid_map_type;
          vid2lvid_map_type vid2lvid_buffer;

          /**
           * \internal
           * The begining id assinged to the first new vertex.
           */
          const lvid_type lvid_start  = graph.vid2lvid.size();

          /**
           * \internal
           * Bit field incidate the vertex that is updated during the ingress.
           */
          dense_bitset updated_lvids(graph.vid2lvid.size());


          /**************************************************************************/
          /*                                                                        */
          /*                         Construct local graph                          */
          /*                                                                        */
          /**************************************************************************/
          { // Add all the edges to the local graph
              graph.local_graph.reserve_edge_space(nedges + 1);

              foreach(const edge_buffer_record& rec, matrix_edges) {
                  // skip re-sent edges
                  if (rec.source == vertex_id_type(-1)) continue;

                  // Get the source_vlid;
                  lvid_type source_lvid(-1);
                  if(graph.vid2lvid.find(rec.source) == graph.vid2lvid.end()) {
                      if (vid2lvid_buffer.find(rec.source) == vid2lvid_buffer.end()) {
                          source_lvid = lvid_start + vid2lvid_buffer.size();
                          vid2lvid_buffer[rec.source] = source_lvid;
                      } else {
                          source_lvid = vid2lvid_buffer[rec.source];
                      }
                  } else {
                      source_lvid = graph.vid2lvid[rec.source];
                      updated_lvids.set_bit(source_lvid);
                  }
                  // Get the target_lvid;
                  lvid_type target_lvid(-1);
                  if(graph.vid2lvid.find(rec.target) == graph.vid2lvid.end()) {
                      if (vid2lvid_buffer.find(rec.target) == vid2lvid_buffer.end()) {
                          target_lvid = lvid_start + vid2lvid_buffer.size();
                          vid2lvid_buffer[rec.target] = target_lvid;
                      } else {
                          target_lvid = vid2lvid_buffer[rec.target];
                      }
                  } else {
                      target_lvid = graph.vid2lvid[rec.target];
                      updated_lvids.set_bit(target_lvid);
                  }
                  graph.local_graph.add_edge(source_lvid, target_lvid, rec.edata);
              } // end for loop over buffers
              matrix_edges.clear();

              ASSERT_EQ(graph.vid2lvid.size() + vid2lvid_buffer.size(),
                        graph.local_graph.num_vertices());
#ifdef TUNING
              if(l_procid == 0)  {
          logstream(LOG_INFO) << "populating local graph: "
                              << ti.current_time()
                              << " secs"
                              << std::endl;
        }
#endif
              // Finalize local graph
              graph.local_graph.finalize();
#ifdef TUNING
              logstream(LOG_INFO) << "local graph info: " << std::endl
                            << "\t nverts: " << graph.local_graph.num_vertices()
                            << std::endl
                            << "\t nedges: " << graph.local_graph.num_edges()
                            << std::endl;

        if(l_procid == 0) {
          logstream(LOG_INFO) << "finalizing local graph: "
                              << ti.current_time()
                              << " secs"
                              << std::endl;
        }
#endif
          }


          /**************************************************************************/
          /*                                                                        */
          /*             Receive and add vertex data to masters                     */
          /*                                                                        */
          /**************************************************************************/
          // Setup the map containing all the vertices being negotiated by this machine
          {
              // receive any vertex data sent by other machines
              if (matrix_vertex_exchange.size() > 0) {
                  vertex_buffer_type vertex_buffer; procid_t sending_proc(-1);
                  while(matrix_vertex_exchange.recv(sending_proc, vertex_buffer)) {
                      foreach(const vertex_buffer_record& rec, vertex_buffer) {
                          lvid_type lvid(-1);
                          if (graph.vid2lvid.find(rec.vid) == graph.vid2lvid.end()) {
                              if (vid2lvid_buffer.find(rec.vid) == vid2lvid_buffer.end()) {
                                  lvid = lvid_start + vid2lvid_buffer.size();
                                  vid2lvid_buffer[rec.vid] = lvid;
                              } else {
                                  lvid = vid2lvid_buffer[rec.vid];
                              }
                          } else {
                              lvid = graph.vid2lvid[rec.vid];
                              updated_lvids.set_bit(lvid);
                          }
                          if (distributed_matrix_block_ingress::vertex_combine_strategy
                              && lvid < graph.num_local_vertices()) {
                              distributed_matrix_block_ingress::vertex_combine_strategy(
                                      graph.l_vertex(lvid).data(), rec.vdata);
                          } else {
                              graph.local_graph.add_vertex(lvid, rec.vdata);
                          }
                      }
                  }
                  matrix_vertex_exchange.clear();
#ifdef TUNING
                  logstream(LOG_INFO) << "base::#vert-msgs=" << matrix_vertex_exchange.size()
                              << std::endl;
          if(l_procid == 0) {
            logstream(LOG_INFO) << "adding vertex data: "
                                << ti.current_time()
                                << " secs"
                                << std::endl;
          }
#endif
              }
          } // end of loop to populate vrecmap


          /**************************************************************************/
          /*                                                                        */
          /*        Assign vertex data and allocate vertex (meta)data  space        */
          /*                                                                        */
          /**************************************************************************/
          {
              // determine masters for all negotiated vertices
              const size_t local_nverts = graph.vid2lvid.size() + vid2lvid_buffer.size();
              graph.lvid2record.reserve(local_nverts);
              graph.lvid2record.resize(local_nverts);
              graph.local_graph.resize(local_nverts);
              foreach(const vid2lvid_pair_type& pair, vid2lvid_buffer) {
                  vertex_record& vrec = graph.lvid2record[pair.second];
                  vrec.gvid = pair.first;
                  if (standalone)
                      vrec.owner = 0;
                  else{
                      if (vht.find(pair.first) != vht.end())
                          vrec.owner = vht[pair.first] ;
                      else if (mht.find(pair.first) != mht.end())
                          vrec.owner = mht[pair.first] ;
                      else vrec.owner = pair.first % cnodes + l_procid/cnodes*cnodes;
                  }
              }
              ASSERT_EQ(local_nverts, graph.local_graph.num_vertices());
              ASSERT_EQ(graph.lvid2record.size(), graph.local_graph.num_vertices());
#ifdef TUNING
              if(l_procid == 0) {
          logstream(LOG_INFO) << "allocating lvid2record: "
                              << ti.current_time()
                              << " secs"
                              << std::endl;
        }
#endif
          }


          /**************************************************************************/
          /*                                                                        */
          /*                          Master handshake                              */
          /*                                                                        */
          /**************************************************************************/
          if (!standalone) {
#ifdef _OPENMP
              buffered_exchange<vertex_id_type> vid_buffer(matrix_rpc.dc(), omp_get_max_threads());
#else
              buffered_exchange<vertex_id_type> vid_buffer(matrix_rpc.dc());
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
              // send not owned vids to their master
              for (lvid_type i = lvid_start; i < graph.lvid2record.size(); ++i) {
                  procid_t master = graph.lvid2record[i].owner;
                  if (master != l_procid)
#ifdef _OPENMP
                      vid_buffer.send(master, graph.lvid2record[i].gvid, omp_get_thread_num());
#else
                      vid_buffer.send(master, graph.lvid2record[i].gvid);
#endif
              }
              vid_buffer.flush();
              matrix_rpc.barrier();

              // receive all vids owned by me
              mutex flying_vids_lock;
              boost::unordered_map<vertex_id_type, mirror_type> flying_vids;
#ifdef _OPENMP
#pragma omp parallel
#endif
              {
                  typename buffered_exchange<vertex_id_type>::buffer_type buffer;
                  procid_t recvid = -1;
                  while(vid_buffer.recv(recvid, buffer)) {
                      foreach(const vertex_id_type vid, buffer) {
                          if (graph.vid2lvid.find(vid) == graph.vid2lvid.end()) {
                              if (vid2lvid_buffer.find(vid) == vid2lvid_buffer.end()) {
                                  flying_vids_lock.lock();
                                  mirror_type& mirrors = flying_vids[vid];
                                  mirrors.set_bit(recvid);
                                  flying_vids_lock.unlock();
                              } else {
                                  lvid_type lvid = vid2lvid_buffer[vid];
                                  graph.lvid2record[lvid]._mirrors.set_bit(recvid);
                              }
                          } else {
                              lvid_type lvid = graph.vid2lvid[vid];
                              graph.lvid2record[lvid]._mirrors.set_bit(recvid);
                              updated_lvids.set_bit(lvid);
                          }
                      }
                  }
              }
              vid_buffer.clear();

              if (!flying_vids.empty()) {
                  logstream(LOG_INFO) << "#flying-own-nverts="
                                      << flying_vids.size()
                                      << std::endl;

                  // reallocate spaces for the flying vertices.
                  size_t vsize_old = graph.lvid2record.size();
                  size_t vsize_new = vsize_old + flying_vids.size();
                  graph.lvid2record.resize(vsize_new);
                  graph.local_graph.resize(vsize_new);
                  for (typename boost::unordered_map<vertex_id_type, mirror_type>::iterator it = flying_vids.begin();
                       it != flying_vids.end(); ++it) {
                      lvid_type lvid = lvid_start + vid2lvid_buffer.size();
                      vertex_record& vrec = graph.lvid2record[lvid];
                      vertex_id_type gvid = it->first;
                      vrec.owner = l_procid;
                      vrec.gvid = gvid;
                      vrec._mirrors = it->second;
                      vid2lvid_buffer[gvid] = lvid;
                  }
              }
          } // end of master handshake

#ifdef TUNING
          if(l_procid == 0) {
        logstream(LOG_INFO) << "master handshake: "
                            << ti.current_time()
                            << " secs"
                            << std::endl;
      }
#endif


          /**************************************************************************/
          /*                                                                        */
          /*                        Merge in vid2lvid_buffer                        */
          /*                                                                        */
          /**************************************************************************/
          {
              if (graph.vid2lvid.size() == 0) {
                  graph.vid2lvid.swap(vid2lvid_buffer);
              } else {
                  graph.vid2lvid.rehash(graph.vid2lvid.size() + vid2lvid_buffer.size());
                  foreach (const typename vid2lvid_map_type::value_type& pair, vid2lvid_buffer) {
                      graph.vid2lvid.insert(pair);
                  }
                  vid2lvid_buffer.clear();
              }
          }


          /**************************************************************************/
          /*                                                                        */
          /*              Synchronize vertex data and meta information              */
          /*                                                                        */
          /**************************************************************************/
          // TODO:  optimization for standalone
          {
              // construct the vertex set of changed vertices

              // Fast pass for first time finalize;
              vertex_set changed_vset(true);

              // Compute the vertices that needs synchronization
              if (!first_time_finalize) {
                  vertex_set changed_vset = vertex_set(false);
                  changed_vset.make_explicit(graph);

                  updated_lvids.resize(graph.num_local_vertices());
                  for (lvid_type i = lvid_start; i <  graph.num_local_vertices(); ++i) {
                      updated_lvids.set_bit(i);
                  }
                  changed_vset.localvset = updated_lvids;
                  buffered_exchange<vertex_id_type> vset_exchange(matrix_rpc.dc());
                  // sync vset with all mirrors
                  changed_vset.synchronize_mirrors_to_master_or(graph, vset_exchange);
                  changed_vset.synchronize_master_to_mirrors(graph, vset_exchange);
              }

              graphlab::graph_gather_apply<graph_type, vertex_negotiator_record>
                      vrecord_sync_gas(graph,
                                       boost::bind(&distributed_matrix_block_ingress::finalize_gather, this, _1, _2),
                                       boost::bind(&distributed_matrix_block_ingress::finalize_apply, this, _1, _2, _3));
              vrecord_sync_gas.exec(changed_vset);

#ifdef TUNING
              if(l_procid == 0) {
          logstream(LOG_INFO) << "synchrionizing vertex (meta)data: "
                              << ti.current_time()
                              << " secs"
                              << std::endl;
        }
#endif
          }

          base_type::exchange_global_info(standalone);
#ifdef TUNING
          if(l_procid == 0) {
        logstream(LOG_INFO) << "exchange global info: "
                            << ti.current_time()
                            << " secs"
                            << std::endl;
      }
#endif

          if(l_procid == 0) {
              memory_info::log_usage("base finalizing done.");
              logstream(LOG_EMPH) << "base finalizing. ("
                                  << ti.current_time()
                                  << " secs)"
                                  << std::endl;
          }
      } // end of modified base finalize


  private:
    boost::function<void(vertex_data_type&, const vertex_data_type&)> vertex_combine_strategy;

    /**
     * \brief Gather the vertex distributed meta data.
     */
    vertex_negotiator_record finalize_gather(lvid_type& lvid, graph_type& graph) {
      vertex_negotiator_record accum;
      accum.num_in_edges = graph.local_graph.num_in_edges(lvid);
      accum.num_out_edges = graph.local_graph.num_out_edges(lvid);
      if (graph.l_is_master(lvid)) {
        accum.has_data = true;
        accum.vdata = graph.l_vertex(lvid).data();
        accum.mirrors = graph.lvid2record[lvid]._mirrors;
      }
      return accum;
    }

    /**
     * \brief Update the vertex data structures with the gathered vertex metadata.
     */
    void finalize_apply(lvid_type lvid, const vertex_negotiator_record& accum, graph_type& graph) {
      typename graph_type::vertex_record& vrec = graph.lvid2record[lvid];
      vrec.num_in_edges = accum.num_in_edges;
      vrec.num_out_edges = accum.num_out_edges;
      graph.l_vertex(lvid).data() = accum.vdata;
      vrec._mirrors = accum.mirrors;
    }
  };

}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>
#endif
