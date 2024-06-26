#ifndef IOPool_Streamer_EventMsgBuilder_h
#define IOPool_Streamer_EventMsgBuilder_h

#include "IOPool/Streamer/interface/MsgTools.h"

// ------------------ event message builder ----------------

namespace edm::streamer {
  class EventMsgBuilder {
  public:
    EventMsgBuilder(void* buf,
                    uint32 size,
                    uint32 run,
                    uint64 event,
                    uint32 lumi,
                    uint32 outModId,
                    uint32 droppedEventsCount,
                    std::vector<bool>& l1_bits,
                    uint8* hlt_bits,
                    uint32 hlt_bit_count,
                    uint32 adler32_chksum,
                    const char* host_name);

    void setOrigDataSize(uint32);
    uint8* startAddress() const { return buf_; }
    void setEventLength(uint32 len);
    void setBufAddr(uint8* buf_addr) { buf_ = buf_addr; }
    void setEventAddr(uint8* event_addr) { event_addr_ = event_addr; }
    uint8* eventAddr() const { return event_addr_; }
    uint32 headerSize() const { return event_addr_ - buf_; }
    uint32 size() const;
    uint32 bufferSize() const { return size_; }

    static uint32 computeHeaderSize(uint32 l1t_bit_count, uint32 hlt_bit_count);

  private:
    uint8* buf_;
    uint32 size_;
    uint8* event_addr_;
  };
}  // namespace edm::streamer
#endif
