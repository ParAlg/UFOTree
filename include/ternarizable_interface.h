#pragma once
#include <cstdint>
#include <utility>


namespace ufo {

template<typename aug_t>
class ITernarizable {
  virtual short get_degree(vertex_t v) = 0;
  virtual std::pair<vertex_t, aug_t> retrieve_v_to_del(vertex_t v) = 0;
};

}
