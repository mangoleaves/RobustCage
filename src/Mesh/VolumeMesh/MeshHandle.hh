#pragma once

namespace Cage
{
namespace VolumeMesh
{
enum class HandleType
{
  Vertex,
  Edge,
  HalfEdge,
  Face,
  HalfFace,
  Cell
};

enum HalfOrientation : int
{
  Sequence = 0,
  Reverse = 1
};

template<HandleType HT>
class Handle
{
public:
  int m_idx;

  Handle() :m_idx(-1) {}
  Handle(int _idx) :m_idx(_idx) {}
  Handle(unsigned int _idx) :m_idx((int)_idx) {}
  Handle(size_t _idx) :m_idx((int)_idx) {}

  inline int idx() const { return m_idx; }
  inline bool isValid() const { return m_idx >= 0; }
  inline void setInvalid() { m_idx = -1; }

  // only compared with handles of same type.
  inline bool operator==(const Handle<HT>& h) const { return m_idx == h.m_idx; }
  inline bool operator!=(const Handle<HT>& h) const { return m_idx != h.m_idx; }
  inline bool operator<(const Handle<HT>& h) const { return m_idx < h.m_idx; }
  inline bool operator<=(const Handle<HT>& h) const { return m_idx <= h.m_idx; }
  inline bool operator>(const Handle<HT>& h) const { return m_idx > h.m_idx; }
  inline bool operator>=(const Handle<HT>& h) const { return m_idx >= h.m_idx; }

  static inline Handle<HT> invalidHandle() { return Handle<HT>(-1); }
};

typedef Handle<HandleType::Vertex> VertexHandle;
typedef Handle<HandleType::Edge> EdgeHandle;
typedef Handle<HandleType::HalfEdge> HalfEdgeHandle;
typedef Handle<HandleType::Face> FaceHandle;
typedef Handle<HandleType::HalfFace> HalfFaceHandle;
typedef Handle<HandleType::Cell> CellHandle;

/* half-elements orientation */
inline HalfOrientation orientation(HalfEdgeHandle heh)
{
  return HalfOrientation(heh.idx() & 0x1);
}
inline HalfOrientation orientation(HalfFaceHandle hfh)
{
  return HalfOrientation(hfh.idx() & 0x1);
}
/* convert between handles and half-handles */
inline HalfEdgeHandle reverse(HalfEdgeHandle heh)
{
  return HalfEdgeHandle(heh.idx() ^ 0x1);
}
inline HalfFaceHandle reverse(HalfFaceHandle hfh)
{
  return HalfFaceHandle(hfh.idx() ^ 0x1);
}
inline EdgeHandle edgeHdl(HalfEdgeHandle heh)
{
  return heh.idx() >> 1;  /* heh.idx / 2*/
}
inline FaceHandle faceHdl(HalfFaceHandle hfh)
{
  return hfh.idx() >> 1; /* hfh.idx / 2*/
}
inline HalfEdgeHandle halfEdgeHdl(EdgeHandle eh, HalfOrientation ori)
{
  return HalfEdgeHandle((eh.idx() << 1) + ori);
}
inline HalfFaceHandle halfFaceHdl(FaceHandle fh, HalfOrientation ori)
{
  return HalfFaceHandle((fh.idx() << 1) + ori);
}
}// namespace VolumeMesh
}// namespace Cage