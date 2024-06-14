#include <cstdint>
#include <utility>

class ITernarizable {
  
  virtual short get_degree(uint32_t v) = 0;
  virtual std::pair<uint32_t, int> retrive_v_to_del(uint32_t v) = 0;

};
