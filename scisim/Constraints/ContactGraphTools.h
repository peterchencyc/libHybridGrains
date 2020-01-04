// ContactGraphTools.h
//
// Breannan Smith
// Last updated: 10/27/2015

#ifndef CONTACT_GRAPH_TOOLS_H
#define CONTACT_GRAPH_TOOLS_H

#include <memory>
#include <vector>

class Constraint;

namespace ContactGraphTools {
void reduceBandwidth(const unsigned long num_bodies,
                     std::vector<std::unique_ptr<Constraint>> &constraints);
}

#endif
