#include <cstddef>
// Get image from clipboard as an RGB image
//
// Outputs:
//   IM  h*w*c list of rgb pixel color values as uint8, running over height,
//     then width, then rgb channels
//   h  height
//   w  width
//   c  channels
bool paste(unsigned char *& IM,size_t & h, size_t & w, size_t & c);
