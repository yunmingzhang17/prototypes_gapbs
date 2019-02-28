#include <algorithm>
#include <cinttypes>

#include "platform_atomics.h"


template<typename PriorityT_>
class PriorityQueue {

public:
  explicit PriorityQueue(bool use_lazy_bucket, PriorityT_* priorities) {
    priorities_ = priorities;
  }

  PriorityT_* priorities_;

};
