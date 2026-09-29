//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef MOOSE_MFEM_ENABLED

#include "CellGridDataCollection.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <vector>

// File-local helpers live in a named namespace (not an anonymous one) so that they cannot clash
// with other files in a unity build. The using-directive is confined to this namespace.
namespace Moose::MFEM::CellGridDetail
{
using namespace mfem;

// ---------------------------------------------------------------------------
// Reference-element data of the VTK discontinuous-Galerkin cells.
//
// All coordinates are in the VTK reference element of each shape: [-1,1]^d
// for edges, quadrilaterals and hexahedra, the unit simplex for triangles and
// tetrahedra, the unit triangle times [-1,1] for wedges, and the pyramid with
// base [-1,1]^2 at t=0 and apex (0,0,1). The corner numbering of every VTK
// shape coincides with the MFEM reference vertex numbering.
//
// The HGRAD tables list the nodal point of each basis function of the VTK
// "C" bases (the bases are Lagrange interpolants). The HCURL/HDIV tables list,
// for each basis function of the VTK "I1" bases (Intrepid2 lowest-order
// Nedelec / Raviart-Thomas), the center of its edge/face, the edge vector
// (v1 - v0) or face normal (cross product of two face edges), and the value
// of the tangential/normal trace of that basis function there; the traces of
// all other basis functions vanish on that edge/face.
//
// The tables below were generated from the VTK 9.7.1 sources in
// Filters/CellGrid/Basis and checked by evaluating the VTK bases.
// ---------------------------------------------------------------------------

struct DofFunctional
{
  real_t x[3];  // edge/face center (VTK reference coordinates)
  real_t d[3];  // edge vector or face normal (VTK reference coordinates)
  real_t trace; // trace of the associated basis function at x along d
};

// clang-format off
static const real_t Edge_HGrad_C1[2][3] =
{
  {-1, 0, 0}, {1, 0, 0}
};
static const real_t Tri_HGrad_C1[3][3] =
{
  {0, 0, 0}, {1, 0, 0}, {0, 1, 0}
};
static const real_t Tri_HGrad_C2[6][3] =
{
  {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0.5, 0, 0},
  {0.5, 0.5, 0}, {0, 0.5, 0}
};
static const real_t Quad_HGrad_C1[4][3] =
{
  {-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}
};
static const real_t Quad_HGrad_C2[9][3] =
{
  {-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0},
  {0, -1, 0}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0},
  {0, 0, 0}
};
static const real_t Tet_HGrad_C1[4][3] =
{
  {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}
};
static const real_t Tet_HGrad_C2[10][3] =
{
  {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
  {0.5, 0, 0}, {0.5, 0.5, 0}, {0, 0.5, 0}, {0, 0, 0.5},
  {0.5, 0, 0.5}, {0, 0.5, 0.5}
};
static const real_t Hex_HGrad_C1[8][3] =
{
  {-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
  {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1}
};
static const real_t Hex_HGrad_C2[27][3] =
{
  {-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
  {-1, -1, 1}, {1, -1, 1}, {1, 1, 1}, {-1, 1, 1},
  {0, -1, -1}, {1, 0, -1}, {0, 1, -1}, {-1, 0, -1},
  {-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0},
  {0, -1, 1}, {1, 0, 1}, {0, 1, 1}, {-1, 0, 1},
  {0, 0, 0}, {0, 0, -1}, {0, 0, 1}, {-1, 0, 0},
  {1, 0, 0}, {0, -1, 0}, {0, 1, 0}
};
static const real_t Wdg_HGrad_C1[6][3] =
{
  {0, 0, -1}, {1, 0, -1}, {0, 1, -1}, {0, 0, 1},
  {1, 0, 1}, {0, 1, 1}
};
static const real_t Wdg_HGrad_C2[18][3] =
{
  {0, 0, -1}, {1, 0, -1}, {0, 1, -1}, {0, 0, 1},
  {1, 0, 1}, {0, 1, 1}, {0.5, 0, -1}, {0.5, 0.5, -1},
  {0, 0.5, -1}, {0, 0, 0}, {1, 0, 0}, {0, 1, 0},
  {0.5, 0, 1}, {0.5, 0.5, 1}, {0, 0.5, 1}, {0.5, 0, 0},
  {0.5, 0.5, 0}, {0, 0.5, 0}
};
static const real_t Pyr_HGrad_C1[5][3] =
{
  {-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0},
  {0, 0, 1}
};
static const DofFunctional Tet_HCurl_I1[6] =
{
  {{0.5, 0, 0}, {1, 0, 0}, 2},
  {{0.5, 0.5, 0}, {-1, 1, 0}, 2},
  {{0, 0.5, 0}, {0, -1, 0}, 2},
  {{0, 0, 0.5}, {0, 0, 1}, 2},
  {{0.5, 0, 0.5}, {-1, 0, 1}, 2},
  {{0, 0.5, 0.5}, {0, -1, 1}, 2}
};
static const DofFunctional Hex_HCurl_I1[12] =
{
  {{0, -1, -1}, {2, 0, 0}, 2},
  {{1, 0, -1}, {0, 2, 0}, 2},
  {{0, 1, -1}, {2, 0, 0}, -2},
  {{-1, 0, -1}, {0, 2, 0}, -2},
  {{0, -1, 1}, {2, 0, 0}, 2},
  {{1, 0, 1}, {0, 2, 0}, 2},
  {{0, 1, 1}, {2, 0, 0}, -2},
  {{-1, 0, 1}, {0, 2, 0}, -2},
  {{-1, -1, 0}, {0, 0, 2}, 2},
  {{1, -1, 0}, {0, 0, 2}, 2},
  {{1, 1, 0}, {0, 0, 2}, 2},
  {{-1, 1, 0}, {0, 0, 2}, 2}
};
static const DofFunctional Wdg_HCurl_I1[9] =
{
  {{0.5, 0, -1}, {1, 0, 0}, 2},
  {{0.5, 0.5, -1}, {-1, 1, 0}, 2},
  {{0, 0.5, -1}, {0, 1, 0}, -2},
  {{0.5, 0, 1}, {1, 0, 0}, 2},
  {{0.5, 0.5, 1}, {-1, 1, 0}, 2},
  {{0, 0.5, 1}, {0, -1, 0}, 2},
  {{0, 0, 0}, {0, 0, 2}, 2},
  {{1, 0, 0}, {0, 0, 2}, 2},
  {{0, 1, 0}, {0, 0, 2}, 2}
};
static const DofFunctional Tet_HDiv_I1[4] =
{
  {{1.0/3, 0, 1.0/3}, {0, -1, 0}, 1},
  {{1.0/3, 1.0/3, 1.0/3}, {1, 1, 1}, 1},
  {{0, 1.0/3, 1.0/3}, {-1, 0, 0}, 1},
  {{1.0/3, 1.0/3, 0}, {0, 0, -1}, 1}
};
static const DofFunctional Hex_HDiv_I1[6] =
{
  {{0, -1, 0}, {0, -4, 0}, 4},
  {{1, 0, 0}, {4, 0, 0}, 4},
  {{0, 1, 0}, {0, 4, 0}, 4},
  {{-1, 0, 0}, {-4, 0, 0}, 4},
  {{0, 0, -1}, {0, 0, -4}, 4},
  {{0, 0, 1}, {0, 0, 4}, 4}
};
static const DofFunctional Wdg_HDiv_I1[5] =
{
  {{0.5, 0, 0}, {0, -2, 0}, 4},
  {{0.5, 0.5, 0}, {2, 2, 0}, 4},
  {{0, 0.5, 0}, {-2, 0, 0}, 4},
  {{1.0/3, 1.0/3, -1}, {0, 0, -1}, 1},
  {{1.0/3, 1.0/3, 1}, {0, 0, 1}, 1}
};
// clang-format on

struct DGShape
{
  Geometry::Type geom;
  const char * type;  // VTK cell-metadata class name
  const char * shape; // VTK shape name
  int nv;             // number of corners
  int max_order;      // highest usable HGRAD order in VTK 9.7
  const real_t (*c1)[3];
  int n1;
  const real_t (*c2)[3];
  int n2;
  const DofFunctional * hcurl;
  int nhcurl;
  const DofFunctional * hdiv;
  int nhdiv;
  real_t A[3]; // diagonal of d(xi_mfem)/d(xi_vtk), if affine-diagonal
};

// Notes on max_order:
// - VTK 9.7 "EdgeC2" has a sign error in its first basis function, so edges
//   are limited to order 1 and subdivided instead.
// - The VTK quadratic pyramid space differs from the MFEM one, so pyramids
//   are limited to order 1 (MFEM's and VTK's linear pyramid spaces agree).
const DGShape *
GetShape(Geometry::Type geom)
{
  // clang-format off
  //  geometry, VTK type, VTK shape, corners, max order,
  //  HGRAD C1 nodes, HGRAD C2 nodes, HCURL I1 dofs, HDIV I1 dofs, d(xi_mfem)/d(xi_vtk)
  static const DGShape shapes[] = {
    {Geometry::SEGMENT, "vtkDGEdge", "edge", 2, 1,
     Edge_HGrad_C1, 2, nullptr, 0, nullptr, 0, nullptr, 0, {0.5, 0, 0}},
    {Geometry::TRIANGLE, "vtkDGTri", "triangle", 3, 2,
     Tri_HGrad_C1, 3, Tri_HGrad_C2, 6, nullptr, 0, nullptr, 0, {1, 1, 0}},
    {Geometry::SQUARE, "vtkDGQuad", "quadrilateral", 4, 2,
     Quad_HGrad_C1, 4, Quad_HGrad_C2, 9, nullptr, 0, nullptr, 0, {0.5, 0.5, 0}},
    {Geometry::TETRAHEDRON, "vtkDGTet", "tetrahedron", 4, 2,
     Tet_HGrad_C1, 4, Tet_HGrad_C2, 10, Tet_HCurl_I1, 6, Tet_HDiv_I1, 4, {1, 1, 1}},
    {Geometry::CUBE, "vtkDGHex", "hexahedron", 8, 2,
     Hex_HGrad_C1, 8, Hex_HGrad_C2, 27, Hex_HCurl_I1, 12, Hex_HDiv_I1, 6, {0.5, 0.5, 0.5}},
    {Geometry::PRISM, "vtkDGWdg", "wedge", 6, 2,
     Wdg_HGrad_C1, 6, Wdg_HGrad_C2, 18, Wdg_HCurl_I1, 9, Wdg_HDiv_I1, 5, {1, 1, 0.5}},
    {Geometry::PYRAMID, "vtkDGPyr", "pyramid", 5, 1,
     Pyr_HGrad_C1, 5, nullptr, 0, nullptr, 0, nullptr, 0, {0, 0, 0}},
  };
  // clang-format on
  for (const DGShape & s : shapes)
  {
    if (s.geom == geom)
    {
      return &s;
    }
  }
  return nullptr;
}

// Map VTK reference coordinates to MFEM reference coordinates.
void
VtkToMfemRef(Geometry::Type geom, const real_t xi[3], IntegrationPoint & ip)
{
  const real_t r = xi[0], s = xi[1], t = xi[2];
  switch (geom)
  {
    case Geometry::SEGMENT:
      ip.Set3((r + 1) / 2, 0, 0);
      break;
    case Geometry::TRIANGLE:
      ip.Set3(r, s, 0);
      break;
    case Geometry::SQUARE:
      ip.Set3((r + 1) / 2, (s + 1) / 2, 0);
      break;
    case Geometry::TETRAHEDRON:
      ip.Set3(r, s, t);
      break;
    case Geometry::CUBE:
      ip.Set3((r + 1) / 2, (s + 1) / 2, (t + 1) / 2);
      break;
    case Geometry::PRISM:
      ip.Set3(r, s, (t + 1) / 2);
      break;
    case Geometry::PYRAMID:
      ip.Set3((r + 1 - t) / 2, (s + 1 - t) / 2, t);
      break;
    default:
      MFEM_ABORT("unsupported geometry");
  }
}

// Maps a point of a parent reference element into one of its sub-cells.
class SubCellMap
{
  const FiniteElement * lin = nullptr;
  const IntegrationRule * pts = nullptr;
  const int * verts = nullptr;
  mutable Vector shape;

public:
  SubCellMap() = default;
  SubCellMap(const FiniteElement * lin_, const IntegrationRule * pts_, const int * verts_)
    : lin(lin_), pts(pts_), verts(verts_), shape(lin_->GetDof())
  {
  }

  void Map(IntegrationPoint & ip) const
  {
    if (!lin)
    {
      return;
    }
    lin->CalcShape(ip, shape);
    real_t x[3] = {0, 0, 0};
    for (int i = 0; i < lin->GetDof(); i++)
    {
      const IntegrationPoint & v = pts->IntPoint(verts[i]);
      x[0] += shape(i) * v.x;
      x[1] += shape(i) * v.y;
      x[2] += shape(i) * v.z;
    }
    ip.Set3(x[0], x[1], x[2]);
  }
};

// ---------------------------------------------------------------------------
// Minimal document model and JSON / MessagePack serialization.
// ---------------------------------------------------------------------------

class DocWriter
{
public:
  virtual ~DocWriter() = default;
  virtual void BeginObject(size_t n) = 0;
  virtual void EndObject() = 0;
  virtual void BeginArray(size_t n) = 0;
  virtual void EndArray() = 0;
  virtual void Key(const std::string & k) = 0;
  virtual void String(const std::string & s) = 0;
  virtual void Bool(bool b) = 0;
  virtual void Int(long long i) = 0;
  virtual void Real(double d) = 0;

  void Reals(const std::vector<double> & v)
  {
    BeginArray(v.size());
    for (double d : v)
    {
      Real(d);
    }
    EndArray();
  }
  void Ints(const std::vector<long long> & v)
  {
    BeginArray(v.size());
    for (long long i : v)
    {
      Int(i);
    }
    EndArray();
  }
};

class JSONWriter : public DocWriter
{
  std::ostream & os;
  std::vector<bool> first;
  bool after_key = false;
  int digits;
  char buf[40];

  void Sep()
  {
    if (after_key)
    {
      after_key = false;
      return;
    }
    if (!first.empty())
    {
      if (!first.back())
      {
        os << ',';
      }
      first.back() = false;
    }
  }

public:
  JSONWriter(std::ostream & os_, int digits_) : os(os_), digits(digits_) {}
  void BeginObject(size_t) override
  {
    Sep();
    os << '{';
    first.push_back(true);
  }
  void EndObject() override
  {
    first.pop_back();
    os << '}';
  }
  void BeginArray(size_t) override
  {
    Sep();
    os << '[';
    first.push_back(true);
  }
  void EndArray() override
  {
    first.pop_back();
    os << ']';
  }
  void Key(const std::string & k) override
  {
    String(k);
    os << ':';
    after_key = true;
  }
  void String(const std::string & s) override
  {
    Sep();
    os << '"';
    for (unsigned char c : s)
    {
      if (c == '"' || c == '\\')
      {
        os << '\\' << c;
      }
      else if (c < 0x20)
      {
        std::snprintf(buf, sizeof(buf), "\\u%04x", c);
        os << buf;
      }
      else
      {
        os << c;
      }
    }
    os << '"';
  }
  void Bool(bool b) override
  {
    Sep();
    os << (b ? "true" : "false");
  }
  void Int(long long i) override
  {
    Sep();
    os << i;
  }
  void Real(double d) override
  {
    Sep();
    if (!std::isfinite(d))
    {
      os << '0';
      return;
    } // JSON has no inf/nan
    std::snprintf(buf, sizeof(buf), "%.*g", digits, d);
    os << buf;
    // Keep the value a JSON floating-point number (e.g. "1" -> "1.0").
    bool has_dot = false;
    for (const char * p = buf; *p; p++)
    {
      if (*p == '.' || *p == 'e' || *p == 'n' || *p == 'i')
      {
        has_dot = true;
      }
    }
    if (!has_dot)
    {
      os << ".0";
    }
  }
};

class MsgPackWriter : public DocWriter
{
  std::ostream & os;
  // Number of items still expected in each open map/array. MessagePack
  // containers are length-prefixed, so a wrong count corrupts the file;
  // verify every count when the container is closed.
  std::vector<size_t> remaining;

  void Put(uint8_t b) { os.put(static_cast<char>(b)); }
  void PutBE(uint64_t v, int nbytes)
  {
    for (int i = nbytes - 1; i >= 0; i--)
    {
      Put(uint8_t(v >> (8 * i)));
    }
  }
  void Header(size_t n, uint8_t fix, size_t fixmax, uint8_t b16, uint8_t b32)
  {
    if (n <= fixmax)
    {
      Put(uint8_t(fix | n));
    }
    else if (n <= 0xffff)
    {
      Put(b16);
      PutBE(n, 2);
    }
    else
    {
      Put(b32);
      PutBE(n, 4);
    }
  }
  void Item()
  {
    if (remaining.empty())
    {
      return;
    }
    MFEM_VERIFY(remaining.back() > 0, "MessagePack container overflow");
    remaining.back()--;
  }
  void Close()
  {
    MFEM_VERIFY(!remaining.empty() && remaining.back() == 0, "MessagePack container size mismatch");
    remaining.pop_back();
  }

public:
  explicit MsgPackWriter(std::ostream & os_) : os(os_) {}
  void BeginObject(size_t n) override
  {
    Item();
    Header(n, 0x80, 15, 0xde, 0xdf);
    remaining.push_back(2 * n);
  }
  void EndObject() override { Close(); }
  void BeginArray(size_t n) override
  {
    Item();
    Header(n, 0x90, 15, 0xdc, 0xdd);
    remaining.push_back(n);
  }
  void EndArray() override { Close(); }
  void Key(const std::string & k) override { String(k); }
  void String(const std::string & s) override
  {
    Item();
    const size_t n = s.size();
    if (n <= 31)
    {
      Put(uint8_t(0xa0 | n));
    }
    else if (n <= 0xff)
    {
      Put(0xd9);
      PutBE(n, 1);
    }
    else if (n <= 0xffff)
    {
      Put(0xda);
      PutBE(n, 2);
    }
    else
    {
      Put(0xdb);
      PutBE(n, 4);
    }
    os.write(s.data(), static_cast<std::streamsize>(n));
  }
  void Bool(bool b) override
  {
    Item();
    Put(b ? 0xc3 : 0xc2);
  }
  void Int(long long i) override
  {
    Item();
    if (i >= 0 && i < 128)
    {
      Put(uint8_t(i));
    }
    else if (i >= 0 && i <= 0xffffffffLL)
    {
      Put(0xce);
      PutBE(uint64_t(i), 4);
    }
    else
    {
      Put(0xd3);
      PutBE(uint64_t(i), 8);
    }
  }
  void Real(double d) override
  {
    Item();
    uint64_t bits;
    std::memcpy(&bits, &d, sizeof(bits));
    Put(0xcb);
    PutBE(bits, 8);
  }
};

struct CGArray
{
  std::string name;
  bool integer = false;
  int comps = 1;
  bool default_scalars = false;
  std::vector<double> d;
  std::vector<long long> i;

  long long Tuples() const
  {
    const size_t n = integer ? i.size() : d.size();
    return comps > 0 ? static_cast<long long>(n / comps) : 0;
  }
};

struct CGGroup
{
  std::string name;
  std::vector<CGArray> arrays;
  CGArray & Add(const std::string & name_, bool integer, int comps)
  {
    arrays.emplace_back();
    CGArray & a = arrays.back();
    a.name = name_;
    a.integer = integer;
    a.comps = comps;
    return a;
  }
};

struct CGCellInfo
{
  std::string space, basis, sharing;
  int order = 1;
  // role -> (group, array)
  using ArrayRef = std::pair<std::string, std::string>; // (group, array)
  std::vector<std::pair<std::string, ArrayRef>> roles;
};

struct CGAttribute
{
  std::string name, space;
  int comps = 1;
  bool shape = false;
  std::vector<std::pair<std::string, CGCellInfo>> info; // cell type -> info
};

struct CGCellType
{
  std::string type, shape, conn_group, conn_array;
};

std::string
SpaceName(int ncomp)
{
  static const char * sup[] = {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};
  return ncomp < 10 ? std::string("ℝ") + sup[ncomp] : "ℝ^" + std::to_string(ncomp);
}

void
WriteArray(DocWriter & w, const CGArray & a)
{
  w.BeginObject(a.default_scalars ? 6 : 5);
  w.Key("name");
  w.String(a.name);
  w.Key("type");
  w.String(a.integer ? "vtktypeint64" : "double");
  w.Key("components");
  w.Int(a.comps);
  w.Key("tuples");
  w.Int(a.Tuples());
  w.Key("data");
  if (a.integer)
  {
    w.Ints(a.i);
  }
  else
  {
    w.Reals(a.d);
  }
  if (a.default_scalars)
  {
    w.Key("default_scalars");
    w.Bool(true);
  }
  w.EndObject();
}

void
WriteDocument(DocWriter & w,
              const std::vector<CGGroup> & groups,
              const std::vector<CGAttribute> & atts,
              const std::vector<CGCellType> & types,
              int cycle,
              real_t time)
{
  w.BeginObject(10);
  w.Key("data-type");
  w.String("cell-grid");
  w.Key("format-version");
  w.Int(1);
  w.Key("schema-name");
  w.String("dg leaf");
  w.Key("schema-version");
  w.Int(0);
  w.Key("content-version");
  w.Int(0);
  w.Key("mfem-cycle");
  w.Int(cycle);
  w.Key("mfem-time");
  w.Real(time);

  w.Key("arrays");
  w.BeginObject(groups.size());
  for (const CGGroup & g : groups)
  {
    w.Key(g.name);
    w.BeginArray(g.arrays.size());
    for (const CGArray & a : g.arrays)
    {
      WriteArray(w, a);
    }
    w.EndArray();
  }
  w.EndObject();

  w.Key("attributes");
  w.BeginArray(atts.size());
  for (const CGAttribute & a : atts)
  {
    w.BeginObject(a.shape ? 5 : 4);
    w.Key("name");
    w.String(a.name);
    w.Key("space");
    w.String(a.space);
    w.Key("components");
    w.Int(a.comps);
    if (a.shape)
    {
      w.Key("shape");
      w.Bool(true);
    }
    w.Key("cell-info");
    w.BeginObject(a.info.size());
    for (const auto & ci : a.info)
    {
      const CGCellInfo & c = ci.second;
      const bool shared = !c.sharing.empty();
      w.Key(ci.first);
      w.BeginObject(4 + (shared ? 1 : 0));
      w.Key("function-space");
      w.String(c.space);
      w.Key("basis");
      w.String(c.basis);
      w.Key("order");
      w.Int(c.order);
      if (shared)
      {
        w.Key("dof-sharing");
        w.String(c.sharing);
      }
      w.Key("arrays");
      w.BeginObject(c.roles.size());
      for (const auto & r : c.roles)
      {
        w.Key(r.first);
        w.BeginArray(2);
        w.String(r.second.first);
        w.String(r.second.second);
        w.EndArray();
      }
      w.EndObject();
      w.EndObject();
    }
    w.EndObject();
    w.EndObject();
  }
  w.EndArray();

  w.Key("cell-types");
  w.BeginArray(types.size());
  for (const CGCellType & t : types)
  {
    w.BeginObject(2);
    w.Key("type");
    w.String(t.type);
    w.Key("cell-spec");
    w.BeginObject(2);
    w.Key("connectivity");
    w.BeginArray(2);
    w.String(t.conn_group);
    w.String(t.conn_array);
    w.EndArray();
    w.Key("shape");
    w.String(t.shape);
    w.EndObject();
    w.EndObject();
  }
  w.EndArray();
  w.EndObject();
}

void
WriteHeader(std::ostream & out, CellGridEncoding enc)
{
  if (enc == CellGridEncoding::MESSAGEPACK)
  {
    out << "vtkCellGrid MessagePack v1\n";
  }
}

std::unique_ptr<DocWriter>
MakeWriter(std::ostream & out, CellGridEncoding enc, int digits)
{
  if (enc == CellGridEncoding::MESSAGEPACK)
  {
    return std::unique_ptr<DocWriter>(new MsgPackWriter(out));
  }
  return std::unique_ptr<DocWriter>(new JSONWriter(out, digits));
}

// How one field is represented on one cell type.
enum class FieldKind
{
  HGRAD,
  CONSTANT,
  HCURL,
  HDIV
};

struct FieldPlan
{
  FieldKind kind = FieldKind::HGRAD;
  int order = 1;
};

// Settings passed from CellGridDataCollection to WriteCellGrid.
struct WriteOptions
{
  int levels_of_detail = 0;
  bool native_vector_spaces = true;
  bool serial = true;
  int myid = 0;
  int precision = 17;
  int cycle = -1;
  real_t time = 0.0;
  CellGridEncoding encoding = CellGridEncoding::MESSAGEPACK;
};

// Write the cell-grid document for the local part of the mesh.
void
WriteCellGrid(std::ostream & out,
              mfem::Mesh * mesh,
              const std::map<std::string, GridFunction *> & field_map,
              const WriteOptions & opt)
{
  const int levels_of_detail = opt.levels_of_detail;
  const bool native_vector_spaces = opt.native_vector_spaces;
  const bool serial = opt.serial;
  const int myid = opt.myid;
  const int precision = opt.precision;
  const int cycle = opt.cycle;
  const real_t time = opt.time;
  const CellGridEncoding encoding = opt.encoding;
  const int dim = mesh->Dimension();
  const int sdim = mesh->SpaceDimension();
  const GridFunction * nodes = mesh->GetNodes();
  const int mesh_order = nodes ? nodes->FESpace()->GetMaxElementOrder() : 1;

  // Elements grouped by geometry, in order of first appearance.
  std::vector<Geometry::Type> geoms;
  std::map<int, std::vector<int>> elems;
  for (int e = 0; e < mesh->GetNE(); e++)
  {
    const Geometry::Type g = mesh->GetElementGeometry(e);
    MFEM_VERIFY(GetShape(g), "unsupported element geometry " << g);
    if (elems[g].empty())
    {
      geoms.push_back(g);
    }
    elems[g].push_back(e);
  }

  std::vector<std::pair<std::string, const GridFunction *>> fields;
  for (const auto & f : field_map)
  {
    fields.emplace_back(f.first, f.second);
  }

  // Highest polynomial order needed and the resulting subdivision.
  int needed = mesh_order;
  for (const auto & f : fields)
  {
    needed = std::max(needed, f.second->FESpace()->GetMaxElementOrder());
  }
  int R = levels_of_detail;
  if (R == 0)
  {
    R = 1;
    for (Geometry::Type g : geoms)
    {
      const DGShape * s = GetShape(g);
      if (g == Geometry::PYRAMID)
      {
        continue;
      }
      R = std::max(R, (needed + s->max_order - 1) / s->max_order);
    }
  }
  if (elems.count(Geometry::PYRAMID) && needed > 1)
  {
    MFEM_WARNING("CellGridDataCollection: pyramids are written with linear"
                 " bases; data of order "
                 << needed << " is interpolated.");
  }

  const bool shared_vertices = (nodes == nullptr && R == 1);

  // Pointers into 'groups' and into each group's arrays are kept while the
  // data is filled in, so reserve all capacity up front.
  std::vector<CGGroup> groups;
  groups.reserve(1 + geoms.size());
  std::vector<CGAttribute> atts;
  std::vector<CGCellType> types;

  // Point group holding the shape DOFs.
  groups.emplace_back();
  groups[0].name = "points";
  CGArray * coords = &groups[0].Add("coords", false, 3);
  coords->default_scalars = true;
  if (shared_vertices)
  {
    coords->d.reserve(3 * mesh->GetNV());
    for (int v = 0; v < mesh->GetNV(); v++)
    {
      const real_t * x = mesh->GetVertex(v);
      for (int c = 0; c < 3; c++)
      {
        coords->d.push_back(c < sdim ? x[c] : 0.0);
      }
    }
  }

  // Attribute skeletons.
  CGAttribute shape_att;
  shape_att.name = "shape";
  shape_att.space = SpaceName(3);
  shape_att.comps = 3;
  shape_att.shape = true;
  std::vector<CGAttribute> field_atts(fields.size());
  for (size_t f = 0; f < fields.size(); f++)
  {
    const GridFunction & gf = *fields[f].second;
    const int mt = gf.FESpace()->GetTypicalFE()->GetMapType();
    int nc = gf.VectorDim();
    if (mt == FiniteElement::H_CURL || mt == FiniteElement::H_DIV)
    {
      nc = std::max(nc, sdim);
    }
    const int out_nc = (nc == 2) ? 3 : nc;
    field_atts[f].name = fields[f].first;
    field_atts[f].comps = out_nc;
    field_atts[f].space = SpaceName(out_nc);
  }
  CGAttribute attr_att, rank_att;
  attr_att.name = "attribute";
  attr_att.space = SpaceName(1);
  rank_att.name = "rank";
  rank_att.space = SpaceName(1);

  mfem::GeometryRefiner refiner(Quadrature1D::ClosedUniform);
  mfem::LinearFECollection lin_fec;
  mfem::Vector val, pval;
  mfem::DenseMatrix Jv(3), adjT(3);

  for (Geometry::Type g : geoms)
  {
    const DGShape & S = *GetShape(g);
    const std::vector<int> & elist = elems[g];

    // Subdivision of this geometry.
    const int Rg = (g == Geometry::PYRAMID) ? 1 : R;
    RefinedGeometry * RG = nullptr;
    int nsub = 1;
    if (Rg > 1)
    {
      RG = refiner.Refine(g, Rg);
      nsub = RG->RefGeoms.Size() / S.nv;
    }
    const long long ncells = (long long)elist.size() * nsub;

    groups.emplace_back();
    CGGroup & grp = groups.back();
    grp.name = S.type;
    grp.arrays.reserve(3 + fields.size());

    // Shape attribute.
    const int shape_order = shared_vertices ? 1 : std::min(mesh_order, S.max_order);
    const real_t(*shape_nodes)[3] = (shape_order == 2) ? S.c2 : S.c1;
    const int n_shape = (shape_order == 2) ? S.n2 : S.n1;

    // Cell-spec connectivity. VTK reads the DOFs of every attribute with
    // shared DOFs (here: the shape) through this array, so it lists all
    // shape nodes of a cell; the corners come first in every VTK basis.
    CGArray & conn = grp.Add("_mfem_conn", true, n_shape);
    conn.i.reserve(ncells * n_shape);

    // Plans for the fields on this geometry.
    std::vector<FieldPlan> plans(fields.size());
    std::vector<CGArray *> farrays(fields.size());
    for (size_t f = 0; f < fields.size(); f++)
    {
      const FiniteElementSpace & fes = *fields[f].second->FESpace();
      int order = 0, mt = FiniteElement::VALUE;
      for (int e : elist)
      {
        const FiniteElement * fe = fes.GetFE(e);
        order = std::max(order, fe->GetOrder());
        mt = fe->GetMapType();
      }
      FieldPlan & p = plans[f];
      const bool vec_space = (mt == FiniteElement::H_CURL || mt == FiniteElement::H_DIV);
      if (vec_space && native_vector_spaces && Rg == 1 && order == 1 && dim == 3 && sdim == 3 &&
          S.hcurl && mesh_order <= S.max_order)
      {
        p.kind = (mt == FiniteElement::H_CURL) ? FieldKind::HCURL : FieldKind::HDIV;
        p.order = 1;
      }
      else if (!vec_space && order == 0)
      {
        p.kind = FieldKind::CONSTANT;
        p.order = 0;
      }
      else
      {
        p.kind = FieldKind::HGRAD;
        p.order = std::max(1, std::min(order, S.max_order));
      }
    }
    for (size_t f = 0; f < fields.size(); f++)
    {
      const FieldPlan & p = plans[f];
      const int out_nc = field_atts[f].comps;
      int comps;
      switch (p.kind)
      {
        case FieldKind::HCURL:
          comps = S.nhcurl;
          break;
        case FieldKind::HDIV:
          comps = S.nhdiv;
          break;
        case FieldKind::CONSTANT:
          comps = out_nc;
          break;
        default:
          comps = out_nc * ((p.order == 2) ? S.n2 : S.n1);
      }
      farrays[f] = &grp.Add(fields[f].first, false, comps);
      farrays[f]->d.reserve(ncells * comps);
    }
    CGArray & attr_arr = grp.Add("_mfem_attribute", false, 1);
    CGArray * rank_arr = serial ? nullptr : &grp.Add("_mfem_rank", false, 1);

    // Loop over elements and their sub-cells.
    IntegrationPoint ip;
    for (int e : elist)
    {
      ElementTransformation * T = mesh->GetElementTransformation(e);
      for (int sc = 0; sc < nsub; sc++)
      {
        SubCellMap sub;
        if (RG)
        {
          sub = SubCellMap(
              lin_fec.FiniteElementForGeometry(g), &RG->RefPts, RG->RefGeoms.GetData() + sc * S.nv);
        }
        const long long first_node = (long long)coords->d.size() / 3;

        // Shape and connectivity.
        if (shared_vertices)
        {
          Array<int> v;
          mesh->GetElementVertices(e, v);
          for (int i = 0; i < S.nv; i++)
          {
            conn.i.push_back(v[i]);
          }
        }
        else
        {
          Vector x(3);
          for (int n = 0; n < n_shape; n++)
          {
            VtkToMfemRef(g, shape_nodes[n], ip);
            sub.Map(ip);
            T->SetIntPoint(&ip);
            x.SetSize(sdim);
            T->Transform(ip, x);
            for (int c = 0; c < 3; c++)
            {
              coords->d.push_back(c < sdim ? x(c) : 0.0);
            }
            conn.i.push_back(first_node + n);
          }
        }

        // Fields.
        for (size_t f = 0; f < fields.size(); f++)
        {
          const GridFunction & gf = *fields[f].second;
          const FieldPlan & p = plans[f];
          const int out_nc = field_atts[f].comps;
          std::vector<double> & data = farrays[f]->d;
          auto push_value = [&](const IntegrationPoint & pt)
          {
            T->SetIntPoint(&pt);
            gf.GetVectorValue(*T, pt, val);
            for (int c = 0; c < out_nc; c++)
            {
              data.push_back(c < val.Size() ? val(c) : 0.0);
            }
          };
          if (p.kind == FieldKind::CONSTANT)
          {
            ip = Geometries.GetCenter(g);
            sub.Map(ip);
            push_value(ip);
          }
          else if (p.kind == FieldKind::HGRAD)
          {
            const real_t(*pts)[3] = (p.order == 2) ? S.c2 : S.c1;
            const int np = (p.order == 2) ? S.n2 : S.n1;
            for (int n = 0; n < np; n++)
            {
              VtkToMfemRef(g, pts[n], ip);
              sub.Map(ip);
              push_value(ip);
            }
          }
          else
          {
            // Lowest-order Nedelec / Raviart-Thomas: apply the DOF
            // functionals dual to the VTK basis. With J the Jacobian
            // w.r.t. VTK reference coordinates, the HCURL coefficient
            // is u.(J d)/trace and the HDIV one u.(adj(J)^T n)/trace.
            const bool curl = (p.kind == FieldKind::HCURL);
            const DofFunctional * dofs = curl ? S.hcurl : S.hdiv;
            const int nd = curl ? S.nhcurl : S.nhdiv;
            for (int k = 0; k < nd; k++)
            {
              VtkToMfemRef(g, dofs[k].x, ip);
              T->SetIntPoint(&ip);
              const mfem::DenseMatrix & J = T->Jacobian();
              for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                {
                  Jv(i, j) = J(i, j) * S.A[j];
                }
              gf.GetVectorValue(*T, ip, val);
              Vector d(const_cast<real_t *>(dofs[k].d), 3), w(3);
              if (curl)
              {
                Jv.Mult(d, w);
              }
              else
              {
                CalcAdjugateTranspose(Jv, adjT);
                adjT.Mult(d, w);
              }
              data.push_back((val * w) / dofs[k].trace);
            }
          }
        }

        attr_arr.d.push_back(mesh->GetAttribute(e));
        if (rank_arr)
        {
          rank_arr->d.push_back(myid);
        }
      }
    }

    // Cell type and attribute metadata.
    types.push_back({S.type, S.shape, S.type, "_mfem_conn"});

    CGCellInfo si;
    si.space = "HGRAD";
    si.basis = "C";
    si.order = shape_order;
    si.sharing = "points";
    si.roles = {{"connectivity", {S.type, conn.name}}, {"values", {"points", "coords"}}};
    shape_att.info.emplace_back(S.type, si);

    for (size_t f = 0; f < fields.size(); f++)
    {
      CGCellInfo ci;
      switch (plans[f].kind)
      {
        case FieldKind::HGRAD:
          ci.space = "HGRAD";
          ci.basis = "C";
          break;
        case FieldKind::CONSTANT:
          ci.space = "constant";
          ci.basis = "C";
          break;
        case FieldKind::HCURL:
          ci.space = "HCURL";
          ci.basis = "I";
          break;
        case FieldKind::HDIV:
          ci.space = "HDIV";
          ci.basis = "I";
          break;
      }
      ci.order = plans[f].order;
      ci.roles = {{"values", {S.type, farrays[f]->name}}};
      field_atts[f].info.emplace_back(S.type, ci);
    }
    CGCellInfo ai;
    ai.space = "constant";
    ai.basis = "C";
    ai.order = 0;
    ai.roles = {{"values", {S.type, "_mfem_attribute"}}};
    attr_att.info.emplace_back(S.type, ai);
    if (rank_arr)
    {
      ai.roles = {{"values", {S.type, "_mfem_rank"}}};
      rank_att.info.emplace_back(S.type, ai);
    }
  }

  atts.push_back(shape_att);
  for (CGAttribute & a : field_atts)
  {
    atts.push_back(a);
  }
  atts.push_back(attr_att);
  if (!serial)
  {
    atts.push_back(rank_att);
  }

  WriteHeader(out, encoding);
  auto w = MakeWriter(out, encoding, precision);
  WriteDocument(*w, groups, atts, types, cycle, time);
}

} // namespace Moose::MFEM::CellGridDetail

namespace Moose::MFEM
{

CellGridDataCollection::CellGridDataCollection(const std::string & collection_name,
                                               mfem::Mesh * mesh_)
  : mfem::DataCollection(collection_name, mesh_)
{
  precision = 17;
}

void
CellGridDataCollection::setLevelsOfDetail(int levels_of_detail)
{
  MFEM_VERIFY(levels_of_detail >= 0, "invalid levels of detail");
  _levels_of_detail = levels_of_detail;
}

const char *
CellGridDataCollection::cellTypeName(mfem::Geometry::Type geom)
{
  const CellGridDetail::DGShape * s = CellGridDetail::GetShape(geom);
  return s ? s->type : nullptr;
}

std::string
CellGridDataCollection::GenerateCollectionPath() const
{
  return prefix_path + GetCollectionName();
}

std::string
CellGridDataCollection::GenerateFileName() const
{
  std::string fname = GetCollectionName();
  if (cycle >= 0)
  {
    fname += "_" + mfem::to_padded_string(cycle, pad_digits_cycle);
  }
  return fname + ".dg";
}

std::string
CellGridDataCollection::GeneratePiecePath() const
{
  return cycle >= 0 ? "Cycle" + mfem::to_padded_string(cycle, pad_digits_cycle)
                    : GetCollectionName() + "_pieces";
}

std::string
CellGridDataCollection::GeneratePieceFileName(int rank) const
{
  return "proc" + mfem::to_padded_string(rank, pad_digits_rank) + ".dg";
}

void
CellGridDataCollection::Load(int)
{
  MFEM_ABORT("CellGridDataCollection does not support Load().");
}

void
CellGridDataCollection::Save()
{
  MFEM_VERIFY(mesh != nullptr, "no mesh set");
  MFEM_VERIFY(mesh->NURBSext == nullptr, "CellGridDataCollection does not support NURBS meshes");
  error = NO_ERROR;
  for (const auto & qf : GetQFieldMap())
  {
    MFEM_WARNING("CellGridDataCollection: skipping QuadratureFunction '" << qf.first << "'");
  }

  // The top-level file is what ParaView opens. Its reader
  // (vtkCompositeCellGridReader, VTK 9.7) parses that file as JSON only, but
  // reads the pieces it lists with vtkCellGridReader, which also accepts
  // MessagePack. So MessagePack data and parallel data are written as pieces
  // with a JSON composite index on top; serial JSON data is a single file.
  const bool use_index = !serial || _encoding == CellGridEncoding::MESSAGEPACK;

  const std::string col_path = GenerateCollectionPath();
  const std::string dir = use_index ? col_path + "/" + GeneratePiecePath() : col_path;
  if (create_directory(dir, mesh, myid) != 0)
  {
    error = WRITE_ERROR;
    MFEM_WARNING("Error creating directory: " << dir);
    return;
  }

  const std::string fname =
      use_index ? dir + "/" + GeneratePieceFileName(myid) : col_path + "/" + GenerateFileName();
  {
    const std::ios::openmode mode = (_encoding == CellGridEncoding::MESSAGEPACK)
                                        ? std::ios::out | std::ios::binary
                                        : std::ios::out;
    std::ofstream out(fname, mode);
    if (!out)
    {
      error = WRITE_ERROR;
      MFEM_WARNING("cannot open " << fname);
    }
    else
    {
      CellGridDetail::WriteOptions opt;
      opt.levels_of_detail = _levels_of_detail;
      opt.native_vector_spaces = _native_vector_spaces;
      opt.serial = serial;
      opt.myid = myid;
      opt.precision = precision;
      opt.cycle = cycle;
      opt.time = time;
      opt.encoding = _encoding;
      CellGridDetail::WriteCellGrid(out, mesh, GetFieldMap(), opt);
    }
    if (!out)
    {
      error = WRITE_ERROR;
    }
  }
  if (!use_index)
  {
    return;
  }

  // Rank 0 writes the composite index listing the non-empty pieces.
  std::vector<int> all_ne(num_procs, 0);
  int ne = mesh->GetNE();
#ifdef MFEM_USE_MPI
  if (!serial)
  {
    MPI_Gather(&ne, 1, MPI_INT, all_ne.data(), 1, MPI_INT, 0, m_comm);
  }
  else
#endif
  {
    all_ne[0] = ne;
  }
  if (myid == 0)
  {
    std::vector<std::string> files;
    for (int r = 0; r < num_procs; r++)
    {
      if (all_ne[r] > 0)
      {
        files.push_back(GeneratePiecePath() + "/" + GeneratePieceFileName(r));
      }
    }
    const std::string iname = col_path + "/" + GenerateFileName();
    std::ofstream out(iname);
    if (!out)
    {
      error = WRITE_ERROR;
      MFEM_WARNING("cannot open " << iname);
      return;
    }
    CellGridDetail::JSONWriter w(out, precision);
    w.BeginObject(4);
    w.Key("data-type");
    w.String("composite");
    w.Key("mfem-cycle");
    w.Int(cycle);
    w.Key("mfem-time");
    w.Real(time);
    w.Key("group");
    w.BeginObject(2);
    w.Key("group-type");
    w.String("collection");
    w.Key("files");
    w.BeginArray(files.size());
    for (const std::string & f : files)
    {
      w.String(f);
    }
    w.EndArray();
    w.EndObject();
    w.EndObject();
    out << '\n';
    if (!out)
    {
      error = WRITE_ERROR;
    }
  }
}

} // namespace Moose::MFEM

#endif
