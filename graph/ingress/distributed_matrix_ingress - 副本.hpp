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

#ifndef GRAPHLAB_DISTRIBUTED_MATRIX_INGRESS_HPP
#define GRAPHLAB_DISTRIBUTED_MATRIX_INGRESS_HPP

#include <boost/functional/hash.hpp>

#include <graphlab/rpc/buffered_exchange.hpp>
#include <graphlab/graph/graph_basic_types.hpp>
#include <graphlab/graph/ingress/distributed_ingress_base.hpp>
#include <graphlab/graph/distributed_graph.hpp>

#include <graphlab/macros_def.hpp>
namespace graphlab {
  template<typename VertexData, typename EdgeData>
  class distributed_graph;

  /**
   * \brief Ingress object assigning edges using matrix.
   */
  template<typename VertexData, typename EdgeData>
  class distributed_matrix_ingress :
    public distributed_ingress_base<VertexData, EdgeData> {
  public:
    typedef distributed_graph<VertexData, EdgeData> graph_type;
    /// The type of the vertex data stored in the graph 
    typedef VertexData vertex_data_type;
    /// The type of the edge data stored in the graph 
    typedef EdgeData   edge_data_type;


    typedef distributed_ingress_base<VertexData, EdgeData> base_type;
   
  public:

    std::vector<std::vector<int> > myMatrix;
	int matrixSize, maxValue;

    distributed_matrix_ingress(distributed_control& dc, graph_type& graph) :
    base_type(dc, graph) {
    	std::fstream fin("test.txt");
    	fin >> matrixSize >> maxValue;
    	myMatrix.resize(matrixSize);
    	for(int i = 0;i<matrixSize;i++){
    		myMatrix.resize(matrixSize);
    		for(int j = 0; j< matrixSize; j++)
    			fin >> myMatrix[i][j];
    	}
    } // end of constructor

    ~distributed_matrix_ingress() { }

    /** \brief Add an edge to the ingress object. */
    virtual void add_edge(vertex_id_type source, vertex_id_type target,
                          const EdgeData& edata) {
        typedef typename base_type::edge_buffer_record edge_buffer_record;

      int row,col;

      int mixingPrime = 1125899906842597L;
      col = (source * mixingPrime) % matrixSize;
      row = (target * mixingPrime) % matrixSize;
      if(source == target) col = rand()%matrixSize;

      const procid_t owning_proc = myMatrix[row][col];
      const edge_buffer_record record(source, target, edata);

      if (mht.find(rec.source) == mht.end()){
          // update mht and nedges_incr
          for (procid_t p = 0; p < nprocs; ++p) {
            if (p != l_procid)
              mht_exchange.send(p, master_pair_type(source, source_owner_proc));
            else
              mht[source] = source_owner_proc;
          }
      }

#ifdef _OPENMP
      base_type::edge_exchange.send(owning_proc, record, omp_get_thread_num());
#else
      base_type::edge_exchange.send(owning_proc, record);
#endif
    } // end of add edge


  }; // end of distributed_random_ingress
}; // end of namespace graphlab
#include <graphlab/macros_undef.hpp>


#endif
