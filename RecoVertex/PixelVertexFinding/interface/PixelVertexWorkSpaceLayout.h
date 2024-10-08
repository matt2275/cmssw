#ifndef RecoVertex_PixelVertexFinding_interface_PixelVertexWorkSpaceLayout_h
#define RecoVertex_PixelVertexFinding_interface_PixelVertexWorkSpaceLayout_h

#include <alpaka/alpaka.hpp>

#include "DataFormats/SoATemplate/interface/SoALayout.h"

// Intermediate data used in the vertex reco algos
// For internal use only
namespace vertexFinder {

  GENERATE_SOA_LAYOUT(PixelVertexWSSoALayout,
                      SOA_COLUMN(uint16_t, itrk),            // index of original track
                      SOA_COLUMN(float, zt),                 // input track z at bs
                      SOA_COLUMN(float, ezt2),               // input error^2 on the above
                      SOA_COLUMN(float, ptt2),               // input pt^2 on the above
                      SOA_COLUMN(uint8_t, izt),              // interized z-position of input tracks
                      SOA_COLUMN(int32_t, iv),               // vertex index for each associated track
                      SOA_SCALAR(uint32_t, ntrks),           // number of "selected tracks"
                      SOA_SCALAR(uint32_t, nvIntermediate))  // the number of vertices after splitting pruning etc.

  using PixelVertexWorkSpaceSoALayout = PixelVertexWSSoALayout<>;
  using PixelVertexWorkSpaceSoAView = PixelVertexWSSoALayout<>::View;
  using PixelVertexWorkSpaceSoAConstView = PixelVertexWSSoALayout<>::ConstView;

  ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE void init(PixelVertexWorkSpaceSoAView& workspace_view) {
    workspace_view.ntrks() = 0;
    workspace_view.nvIntermediate() = 0;
  }

}  // namespace vertexFinder

#endif  // RecoVertex_PixelVertexFinding_interface_PixelVertexWorkSpaceLayout_h
