#include "vrt/vrt_read.h"
#include <iostream>

int main(int argc, char** argv) {
  std::cout<<"Tammo says hi\n";
  const uint32_t b[2] = {1, 2};
  uint64_t c = read_uint64(b);
  return 1;
}
