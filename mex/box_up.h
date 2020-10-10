// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2019 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_COPYLEFT_CGAL_BOX_UP_H
#define IGL_COPYLEFT_CGAL_BOX_UP_H
#include <Eigen/Core>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/box_intersection_d.h>

namespace igl
{
  namespace copyleft
  {
    namespace cgal
    {
      // Box up segments
      //
      // Inputs:
      //   S  #S list of CGAL segments (non-const for iterator)
      // Outputs:
      //   B  #S list of CGAL Boxes
      template <
        typename Kernel>
      inline void box_up(
        /*const*/ std::vector<CGAL::Segment_2<Kernel> > & S,
        std::vector<
        CGAL::Box_intersection_d::Box_with_handle_d<
          double,2,typename std::vector<CGAL::Segment_2<Kernel> >::iterator > > & B);
    }
  }
}

// Implementation

template <
  typename Kernel>
inline void igl::copyleft::cgal::box_up(
  std::vector<CGAL::Segment_2<Kernel> > & S,
  std::vector<
  CGAL::Box_intersection_d::Box_with_handle_d<
    double,2,typename std::vector<CGAL::Segment_2<Kernel> >::iterator > > & B)
{
  typedef typename std::vector<CGAL::Segment_2<Kernel> > Segments;
  typedef typename Segments::iterator SegmentsIterator;
  typedef typename Segments::const_iterator SegmentsConstIterator;
  typedef 
    CGAL::Box_intersection_d::Box_with_handle_d<double,2,SegmentsIterator> 
    Box;
  B.reserve(S.size());
  for ( 
    SegmentsIterator sit = S.begin(); 
    sit != S.end(); 
    ++sit)
  {
    B.push_back(Box(sit->bbox(), sit));
  }
}

#endif


