#pragma once

#include "defs.h"
#include "parser.h"
#include "shortest_certificate.h"
#include "geometry_basics.h"
#include "certificate.h"
#include "frechet_light.h"

#include <list>
#include <string>
#include <array>
#include <vector>
#include <set>
#include <tuple>
#include <algorithm>
#include <queue>
#include <map>
#include <iterator>


class ShortestCertificate {

public:
struct Freespace_Node {
    
private:
    std::vector<Freespace_Node*> m_successor;
    CInterval m_interval;
           
public:
    Freespace_Node* parent;
    PointID c1;
    PointID c2;            
    u_int32_t steps;
    bool visited;  
    bool direction; // true->horizontal, false->vertical

public:
    Freespace_Node()
        : m_interval(CInterval(1000000, 0.0, 1000000, 0.42)),
          m_successor(std::vector<Freespace_Node*>()),
          c1(1000000),
          c2(1000000),
          visited(false),
          parent(nullptr),
          steps(0),
          direction(true) {}

    Freespace_Node(CInterval interval, PointID c1_point, PointID c2_point, bool dir) 
        : m_interval(interval), 
          m_successor(std::vector<Freespace_Node*>()),
          c1(c1_point),
          c2(c2_point),
          visited(false),
          steps(0),
          direction(dir) {}

    void add_edge(Freespace_Node* successor) {
        m_successor.push_back(successor);
    }

    std::vector<Freespace_Node*> get_edges() {
        return m_successor;
    }

    std::vector<Freespace_Node*> get_successor_list() {
        return m_successor;
    }

    size_t num_edges() {
        return m_successor.size();
    }

    PointID get_pointid() {
        return m_interval.begin.getPoint();
    }

    CInterval get_interval() {
        return m_interval;
    }

    void setEndFraction(distance_t frac) {
        m_interval.end.setFraction(frac);
    }

    void add_parent(Freespace_Node* parent_node) {
        parent = parent_node;
    }

    bool is_visited() {
        return visited;
    }

};
public:
struct Box_Set {
    
public:
    std::vector<Freespace_Node> horizontal;
    std::vector<Freespace_Node> vertical;

    Box_Set() : horizontal(std::vector<Freespace_Node>()), vertical(std::vector<Freespace_Node>()) {};
};

public:
    Box_Set** matrix;
	Curve curve1;
    Curve curve2;   
    distance_t delta;

public:
    ShortestCertificate(std::string c_file1, std::string c_file2, std::string dist_delta) :
    delta(std::stod(dist_delta))
    {
        curve1 = parser::readCurve(c_file1);
	    curve2 = parser::readCurve(c_file2);

        const size_t CURVE1_SIZE = curve1.size();
        const size_t CURVE2_SIZE = curve2.size();

        matrix = new Box_Set*[CURVE2_SIZE];
        for(size_t i = 0; i < CURVE2_SIZE; i++) matrix[i] = new Box_Set[CURVE1_SIZE];
    }

    ~ShortestCertificate() {
        for(size_t i = 0; i < curve2.size()-1; i++) delete[] matrix[i];
        delete[] matrix;
    }

public:
    void printUsage();
    Certificate no_certificate(Curve& curve1, Curve& curve2, Box_Set** matrix, double delta);
    Certificate yes_certificate(Curve& curve1, Curve& curve2, Box_Set** matrix, double delta);
};
