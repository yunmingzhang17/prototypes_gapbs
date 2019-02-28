#include <algorithm>
#include <cinttypes>

#include "platform_atomics.h"


/**
 * Phase-synchronous priority queue with dual representation
 * Representation 1: When using thread-local buckets, there is nothing stored in the data strucutre. It merely holds the current bucket index, next bucket index and other metadata. The real priority queue is distributed across threads. 
 * Representation 2: When using lazy buckets, the priority queue actually stores all the nodes with their buckets (the buckets are not distributed)
 **/
template<typename PriorityT_>
class PriorityQueue {

public:
  explicit PriorityQueue(bool use_lazy_bucket, PriorityT_* priorities) {
    priorities_ = priorities;
  }
  
  bool finished() {
    
  }

  void updatePriorityMin(NodeID dst, PriorityT_ new_p, PriorityT_ old_p){

  }

  PriorityT_* priorities_;
  const PriorityT_ kDistInf = numeric_limits<PriorityT_>::max()/2;
  const size_t kMaxBin = numeric_limits<size_t>::max()/2;

};
