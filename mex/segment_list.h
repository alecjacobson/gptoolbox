// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2019 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COPYLEFT_CGAL_SEGMENT_LIST_H
#define IGL_COPYLEFT_CGAL_SEGMENT_LIST_H
#include <Eigen/Core>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      // Convert a edge-complex (V,E) to a list of CGAL segments
      //
      // Templates:
      //   Kernal  CGAL computation and construction kernel (e.g.
      //     CGAL::Exact_predicates_exact_constructions_kernel)
      // Inputs:
      //   V  #V by 2 list of vertex positions
      //   E  #E by 2 list of triangle indices
      // Outputs:
      //   S  #E list of CGAL segments
      template <
        typename DerivedV,
        typename DerivedE,
        typename Kernel>
      inline void segment_list(
        const Eigen::MatrixBase<DerivedV> & V,
        const Eigen::MatrixBase<DerivedE> & E,
        std::vector<CGAL::Segment_2<Kernel> > & S);
    }
  }
}

// Implementation
template <
  typename DerivedV,
  typename DerivedE,
  typename Kernel>
inline void igl::copyleft::cgal::segment_list(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedE> & E,
  std::vector<CGAL::Segment_2<Kernel> > & S)
{
  S.reserve(E.rows());
  for(int i = 0;i<E.rows();i++)
  {
    S.emplace_back(
      CGAL::Point_2<Kernel>(V(E(i,1),0),V(E(i,1),1)),
      CGAL::Point_2<Kernel>(V(E(i,0),0),V(E(i,0),1)));
  }
}

#endif

