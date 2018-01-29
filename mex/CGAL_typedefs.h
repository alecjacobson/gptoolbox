#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

// Precision Kernel
typedef CGAL::Simple_cartesian<double> K;

typedef K::FT FT;
typedef K::Ray_3      Ray3;
typedef K::Line_3     Line3;
typedef K::Point_3    Point3;
typedef K::Vector_3   Vector3;
typedef K::Triangle_3 Triangle3;
typedef K::Segment_2 Segment2;

typedef std::list<Triangle3> TriangleList;
typedef TriangleList::iterator TriangleListIterator;
typedef TriangleList::const_iterator TriangleListConstIterator;
typedef CGAL::AABB_triangle_primitive<K,TriangleListIterator> Primitive;
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> TriTree;
typedef TriTree::Object_and_primitive_id Object_and_primitive_id;
typedef TriTree::Primitive_id Primitive_id;

typedef std::list<Segment2> SegmentList;
typedef SegmentList::iterator SegmentListIterator;
typedef SegmentList::const_iterator SegmentListConstIterator;
//typedef CGAL::AABB_segment_primitive<K,SegmentListIterator> SegPrimitive;
//typedef CGAL::AABB_traits<K, SegPrimitive> AABB_segment_traits;
//typedef CGAL::AABB_tree<AABB_segment_traits> SegTree;
//typedef SegTree::Object_and_primitive_id Object_and_primitive_id;
//typedef SegTree::Primitive_id Primitive_id;

#endif
