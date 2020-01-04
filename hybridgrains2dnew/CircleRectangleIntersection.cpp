//
//  CircleRectangleIntersection.cpp
//  CircleRasterizer
//
//  Created by Yonghao Yue on 7/13/16.
//  Copyright © 2016 Yonghao Yue. All rights reserved.
//

#include "CircleRectangleIntersection.h"
#include <assert.h>
#include <iostream>
#include <math.h>

// http://stackoverflow.com/questions/622287/area-of-intersection-between-circle-and-rectangle

// returns the positive root of intersection of line y = h with circle centered
// at the origin and radius r
static double section(const double &h, const double &r = 1) {
  // assume r is positive, leads to some simplifications in the formula below
  // (can factor out r from the square root)
  assert(r >= 0);
  // http://www.wolframalpha.com/input/?i=r+*+sin%28acos%28x+%2F+r%29%29+%3D+h
  return (h < r) ? sqrt(r * r - h * h) : 0;
}

// indefinite integral of circle segment
static double g(const double &x, const double &h, const double &r = 1) {
  // http://www.wolframalpha.com/input/?i=r+*+sin%28acos%28x+%2F+r%29%29+-+h
  return .5 *
         (sqrt(1 - x * x / (r * r)) * x * r + r * r * asin(x / r) - 2 * h * x);
}

// area of intersection of an infinitely tall box with left edge at x0, right
// edge at x1, bottom edge at h and top edge at infinity, with circle centered
// at the origin with radius r
static double area(double x0, double x1, const double &h, const double &r) {
  if (x0 > x1)
    std::swap(x0, x1); // this must be sorted otherwise we get negative area
  double s = section(h, r);
  return g(std::max<double>(-s, std::min<double>(s, x1)), h, r) -
         g(std::max<double>(-s, std::min<double>(s, x0)), h,
           r); // integrate the area
}

// area of the intersection of a finite box with a circle centered at the origin
// with radius r
static double area(double x0, double x1, double y0, double y1,
                   const double &r) {
  if (y0 > y1)
    std::swap(y0, y1); // this will simplify the reasoning
  if (y0 < 0) {
    if (y1 < 0)
      return area(
          x0, x1, -y0, -y1,
          r); // the box is completely under, just flip it above and try again
    else
      return area(x0, x1, 0, -y0, r) +
             area(x0, x1, 0, y1, r); // the box is both above and below, divide
                                     // it to two boxes and go again
  } else {
    assert(y1 >= 0); // y0 >= 0, which means that y1 >= 0 also (y1 >= y0)
                     // because of the swap at the beginning
    return area(x0, x1, y0, r) -
           area(x0, x1, y1,
                r); // area of the lower box minus area of the higher box
  }
}

double circleRectangleIntersectionArea(
    double x0, double x1, double y0, double y1, const double &cx,
    const double &cy,
    const double
        &r) // area of the intersection of a general box with a general circle
{
  x0 -= cx;
  x1 -= cx;
  y0 -= cy;
  y1 -= cy;
  // get rid of the circle center

  return area(x0, x1, y0, y1, r);
}
