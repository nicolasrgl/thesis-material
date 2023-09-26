#include "shortest_certificate.h"


inline bool inFreeSpace(Curve& curve1, Curve& curve2, CPoint p1, CPoint p2, distance_t delta) {

    distance_t point_distance = std::sqrt(std::pow(curve1.interpolate_at(p1).y - curve2.interpolate_at(p2).y, 2) + std::pow(curve1.interpolate_at(p1).x - curve2.interpolate_at(p2).x, 2));
    return (point_distance <= delta);
}

inline void splitCurve(Curve& curve, std::vector<std::vector<distance_t>>& approximations, std::vector<PointID>& curve_points, FrechetLight& frechet, PointID index) {

    if(curve.size() <= 2) return;
    distance_t costs[curve.size()];

    for(PointID i = 2; i < curve.size()-1; ++i) {

        Curve constant_curve_left({curve[0], curve[i]});                   
        Curve partial_curve_left = Curve(std::vector<Point>(curve.begin(), curve.begin() + i));
        approximations[index][index+i] = frechet.calcDistance(constant_curve_left, partial_curve_left);
        
        Curve constant_curve_right({curve[i], curve[curve.size()-1]});                   
        Curve partial_curve_right = Curve(std::vector<Point>(curve.begin() + i, curve.end()));
        approximations[index+i][index+curve.size()-1] = frechet.calcDistance(constant_curve_right, partial_curve_right);
        
        costs[i] = std::max(approximations[index][index+i], approximations[index+i][index+ curve.size()-1]);
    }

    PointID min = 2;
    for(size_t j = 3; j < curve.size()-1; ++j) {
        if(costs[min] > costs[j]) min = j;
    }

    curve_points.push_back(min + index);
    Curve left = Curve(std::vector<Point>(curve.begin(), curve.begin() + min));
    Curve right = Curve(std::vector<Point>(curve.begin() + min, curve.end()));

    splitCurve(left, approximations, curve_points, frechet, index);
    splitCurve(right, approximations, curve_points, frechet, min+index);
}

using SuccessorList = std::list<std::pair<PointID, PointID>>;
inline std::vector<std::vector<SuccessorList>> find_jumps(Curve& curve1, Curve& curve2, double delta) {

    const size_t N = curve1.size();
    const size_t M = curve2.size();
    FrechetLight frechet;
    std::vector<std::vector<SuccessorList>> adj_matrix(M, std::vector<SuccessorList>(N));

    std::vector<std::vector<distance_t>> approx_x(N, std::vector<distance_t>(N));
    std::vector<std::vector<distance_t>> approx_y(M, std::vector<distance_t>(M));

    std::vector<PointID> curve_points_1;
    std::vector<PointID> curve_points_2;

    splitCurve(curve1, approx_x, curve_points_1, frechet, 0);
    splitCurve(curve2, approx_y, curve_points_2, frechet, 0);

    curve_points_1.push_back(0);
    curve_points_1.push_back(curve1.size()-1);
    approx_x[0][curve1.size()-1] = delta;
    curve_points_2.push_back(0);
    curve_points_2.push_back(curve2.size()-1);
    approx_y[0][curve2.size()-1] = delta;

    std::sort(curve_points_1.begin(), curve_points_1.end());
    std::sort(curve_points_2.begin(), curve_points_2.end());

    for(size_t i = 0; i < curve_points_1.size(); i++) {
        for(size_t iter = i+2; iter < curve_points_1.size(); iter++) {

            Curve constant_curve({curve1[curve_points_1[i]], curve1[curve_points_1[iter]]});                   
            Curve partial_curve = Curve(std::vector<Point>(curve1.begin() + curve_points_1[i], curve1.begin() + curve_points_1[iter]));
            approx_x[curve_points_1[i]][curve_points_1[iter]] = frechet.calcDistance(constant_curve, partial_curve);
        }
    }

    for(size_t j = 0; j < curve_points_2.size(); j++) {
        for(size_t iter = j+2; iter < curve_points_2.size(); iter++) {

            Curve constant_curve({curve2[curve_points_2[j]], curve2[curve_points_2[j]]});                   
            Curve partial_curve = Curve(std::vector<Point>(curve2.begin() + curve_points_2[j], curve2.begin() + curve_points_2[iter]));
            approx_y[curve_points_2[j]][curve_points_2[iter]] = frechet.calcDistance(constant_curve, partial_curve);
        }
    }

    for(size_t j = 0; j < curve_points_2.size(); ++j) {
        for(size_t i = 0; i < curve_points_1.size(); ++i) {

            for(size_t current_y = j; current_y < curve_points_2.size(); ++current_y) { //iterate for each point O(n^2)
                for(size_t current_x = i; current_x < curve_points_1.size(); ++current_x) {

                    //compute distance between length 2 curves O(1)
                    distance_t end_approx = std::sqrt(std::pow(curve2[curve_points_2[current_y]].y - curve1[curve_points_1[current_x]].y, 2) + std::pow(curve2[curve_points_2[current_y]].x - curve1[curve_points_1[current_x]].x, 2));
                    distance_t start_approx = std::sqrt(std::pow(curve2[curve_points_2[j]].y - curve1[curve_points_1[i]].y, 2) + std::pow(curve2[curve_points_2[j]].x - curve1[curve_points_1[i]].x, 2));
                    distance_t line_approx = std::max(start_approx, end_approx);

                    if(approx_y[curve_points_2[j]][curve_points_2[current_y]] + approx_x[curve_points_1[i]][curve_points_1[current_x]] + line_approx <= delta) {               
                        adj_matrix[curve_points_2[j]][curve_points_1[i]].push_back(std::make_pair(curve_points_2[current_y], curve_points_1[current_x]));
                    }
                }
            }
        }
    }

    return adj_matrix;
}

inline Certificate pruning_diagonals(Certificate c, Curve curve1, Curve curve2, distance_t delta) {

    std::vector<std::pair<CPoint, CPoint>> traversal;
    FrechetLight frechet;

    // get all traversal PointIDs
    for(CPosition point : c.getTraversal()) {

       traversal.push_back(std::make_pair(point[0], point[1]));
    }
    
    std::map<size_t, std::set<size_t>> shortest_traversal;

    for(size_t i = 0; i < traversal.size()-1; i++) {
                shortest_traversal[i].insert(i+1);
            }

    for(size_t j = 0; j < traversal.size(); ++j) {
        for(size_t i = j+2; i < traversal.size(); ++i) {

            auto start = traversal[j];
            auto end = traversal[i];

            distance_t x1 = (double) traversal[j].first.getPoint() + traversal[j].first.getFraction();
            distance_t y1 = (double) traversal[j].second.getPoint() + traversal[j].second.getFraction();

            distance_t x2 = (double) traversal[i].first.getPoint() + traversal[i].first.getFraction();
            distance_t y2 = (double) traversal[i].second.getPoint() + traversal[i].second.getFraction();

            //compute gradient for curve pieces
            double segment_size_1 = x2 - x1;
            double segment_size_2 = y2 - y1;

            //create line function y = m * x + t
            double m = (segment_size_1 < segment_size_2) ? segment_size_2 / segment_size_1 : (segment_size_2 == 0) ? 0 : segment_size_1 / segment_size_2;
            double t = (segment_size_1 < segment_size_2) ? y2 - m * x2 : x2 - m * y2;
            bool valid = true;

            for(double i2 = std::ceil(x1); i2 <= std::floor(x2); i2++) {
                double j2 = m * i2 + t;
                double floor_j = std::floor(j2);

                CPoint p1((uint32_t) i2, 0.0);
                CPoint p2((uint32_t) floor_j, j2 - floor_j);

                if(j2 >= curve2.size()) {
                    valid = false;
                    break;
                }
                if(!inFreeSpace(curve1, curve2, p1, p2, delta)) {
                    valid = false;
                }
            }
            if(m != 0) m = 1/m;
            t = (segment_size_1 < segment_size_2) ? x2 - m * y2 : y2 - m * x2;

            for(double i2 = std::ceil(y1); i2 <= std::floor(y2); i2++) {
                double j2 = m * i2 + t;
                double floor_j = std::floor(j2);

                CPoint p2((uint32_t) i2, 0.0);
                CPoint p1((uint32_t) floor_j, j2 - floor_j);

                if(j2 >= curve1.size()) {
                    valid = false;
                    break;
                }
                if(!inFreeSpace(curve1, curve2, p1, p2, delta)) {
                    valid = false;
                }
            }
            if(valid) shortest_traversal[j].insert(i);
        }
    }

    size_t parent[traversal.size()];
    bool visited[traversal.size()];
    memset(visited, 0, sizeof(visited));

    parent[0] = UINT32_MAX;
    std::queue<size_t> q;
    q.push(0);
    visited[0] = 1;
    long counter = 0;

    while(!q.empty()) {
        auto index = q.front();
        q.pop(); 
        counter++;
        for(auto succ : shortest_traversal[index]) {
            if(!visited[succ]) {

                visited[succ] = 1;
                q.push(succ);
                parent[succ] = index;
            }
        }
    }

    std::vector<size_t> indices;
    size_t iter = traversal.size()-1;
    indices.push_back(iter);

    while(parent[iter] != UINT32_MAX) {
        indices.push_back(parent[iter]);
        iter = parent[iter];
    }

    std::reverse(indices.begin(), indices.end());

    Certificate c_ret;
    c_ret.setCurves(&curve1, &curve2);
    c_ret.setDistance(delta);
    for(auto index : indices) c_ret.addPoint({traversal[index].first, traversal[index].second});
    c_ret.validate();
    c_ret.setAnswer(1);

    return c_ret;
}

Certificate ShortestCertificate::no_certificate(Curve& curve1, Curve& curve2, Box_Set** matrix, double delta) {

    size_t n = curve1.size(), m = curve2.size();
    distance_t x = curve2[0].x - curve1[0].x, y = curve2[0].y - curve1[0].y;

    if(std::sqrt(x*x + y*y) > delta) {
        Certificate c;
        c.setCurves(&curve1, &curve2);
        c.setDistance(delta);
        c.addPoint({CPoint(0, 0.0), CPoint(0, 0.0)});
        c.validate();
        c.setAnswer(false);
        c.dump_certificate();
        return c;
    }
    x = curve2[curve2.size()-1].x - curve1[curve1.size()-1].x;
    y = curve2[curve2.size()-1].y - curve1[curve1.size()-1].y;
    if(std::sqrt(x*x + y*y) > delta) {
        Certificate c;
        c.setCurves(&curve1, &curve2);
        c.setDistance(delta);
        c.addPoint({CPoint(curve1.size()-1, 0.0), CPoint(curve2.size()-1, 0.0)});
        c.validate();
        c.setAnswer(false);
        c.dump_certificate();
        return c;
    }

    for(PointID j = 0; j < m; ++j) {   
        for(PointID i = n-1; i > 0; --i) {

            Interval interval_horizontal = IntersectionAlgorithm::intersection_interval(curve2[j], delta, curve1[i-1], curve1[i]);

            if(interval_horizontal.is_empty()) {
                // add Node to grid with interval from [0, 1]
                matrix[j][i-1].horizontal.push_back(Freespace_Node(CInterval(i-1, 0.0, i-1, 1.0), i-1, j, true));
            }
            else if(interval_horizontal.begin != 0.0 && interval_horizontal.end == 1.0) {
                // add Node to grid with interval from [0, begin)
                matrix[j][i-1].horizontal.push_back(Freespace_Node(CInterval(i-1, 0.0, i-1, interval_horizontal .begin), i-1, j, true));
            }
            else if(interval_horizontal.begin == 0.0 && interval_horizontal.end != 1.0) {
                // add Node to grid with interval from (end, 1]
                matrix[j][i-1].horizontal.push_back(Freespace_Node(CInterval(i-1, interval_horizontal.end, i-1, 1.0), i-1, j, true));            
            }
            else if(interval_horizontal.begin != 0.0 && interval_horizontal.end != 1.0 && !interval_horizontal.is_empty()) {
                // add 2 (!) nodes to grid with separated intervals
                matrix[j][i-1].horizontal.push_back(Freespace_Node(CInterval(i-1, 0.0, i-1, interval_horizontal.begin), i-1, j, true));
                matrix[j][i-1].horizontal.push_back(Freespace_Node(CInterval(i-1, interval_horizontal.end, i-1, 1.0), i-1, j, true));
            }
        }
    }

    for(PointID i = n; i > 0; --i) {  
        for(PointID j = 0; j < m-1; ++j) { 

            Interval interval_vertical = IntersectionAlgorithm::intersection_interval(curve1[i-1], delta, curve2[j], curve2[j+1]);

            if(interval_vertical.is_empty()) {
                // add Node to grid with interval from [0, 1]
                matrix[j][i-1].vertical.push_back(Freespace_Node(CInterval(i, 0.0, i, 1.0), i-1, j, false));
            }
            else if(interval_vertical.begin != 0.0 && interval_vertical .end == 1.0) {
                // add Node to grid with interval from [0, begin)
                matrix[j][i-1].vertical.push_back(Freespace_Node(CInterval(i, 0.0, i, interval_vertical.begin), i-1, j, false));
            }
            else if(interval_vertical.begin == 0.0 && interval_vertical .end != 1.0) {
                // add Node to grid with interval from (end, 1]
                matrix[j][i-1].vertical.push_back(Freespace_Node(CInterval(i, interval_vertical.end, i-1, 1.0), i-1, j, false));
            }
            else if(interval_vertical .begin != 0.0 && interval_vertical.end != 1.0 && !interval_vertical.is_empty()) {
                // add 2 (!) nodes to grid with separated intervals
                matrix[j][i-1].vertical.push_back(Freespace_Node(CInterval(i, 0.0, i, interval_vertical.begin), i-1, j, false));
                matrix[j][i-1].vertical.push_back(Freespace_Node(CInterval(i, interval_vertical.end, i, 1.0), i-1, j, false));
            }
        }
    }

    for(PointID j = 0; j < m; ++j) {
        for(PointID i = n; i > 0; --i) {

            //HORIZONTAL
            if(!matrix[j][i-1].horizontal.empty()) {
            
                Freespace_Node* currentNode_horizontal = &(matrix[j][i-1].horizontal[0]);
                
                if(!matrix[j][i-1].vertical.empty()) {
                    Freespace_Node* currentNode_vertical = &(matrix[j][i-1].vertical[0]);

                    //single check neighbouring vertical interval if connected
                    if(currentNode_horizontal->get_interval().begin.getFraction() == 0.0 && currentNode_vertical->get_interval().begin.getFraction() == 0.0) {
                        currentNode_horizontal->add_edge(currentNode_vertical);
                    }
                }

                //iterate left until 0 and add potential edges, break when empty or green gap                
                for(PointID row = i; row > 0; --row) {

                    if(matrix[j][row-1].horizontal.empty()) {
                        break;
                    }

                    distance_t begin_fraction = matrix[j][row-1].horizontal[0].get_interval().begin.getFraction();
                    distance_t end_fraction = matrix[j][row-1].horizontal[0].get_interval().end.getFraction();                
                    
                    if(begin_fraction == 0.0 && end_fraction == 0.0) {
                        currentNode_horizontal->add_edge(&(matrix[j][row-1].horizontal[0]));
                    }
                    else if(matrix[j][row-1].horizontal.size() == 2) {
                        currentNode_horizontal->add_edge(&(matrix[j][row-1].horizontal[1]));
                        break;
                    }
                    else if(begin_fraction != 0.0 && end_fraction == 0.0) {              
                        currentNode_horizontal->add_edge(&(matrix[j][row-1].horizontal[0]));
                        break;
                    }
                    else if(begin_fraction == 0.0 && end_fraction != 0.0) {
                        break;
                    }                   
                }          
            }

            //VERTICAL
            if(!matrix[j][i-1].vertical.empty()) {
            
                size_t oneOrTwo = matrix[j][i-1].vertical.size() == 2 ? 1 : 0;
                Freespace_Node* currentNode_vertical = &(matrix[j][i-1].vertical[oneOrTwo]);
                
                if(j+1 < m && i-2 != std::numeric_limits<uint32_t>::max() && !matrix[j+1][i-2].horizontal.empty()) {
                    Freespace_Node* currentNode_horizontal = &(matrix[j+1][i-2].horizontal[0]);

                    //single check neighbouring vertical interval if connected
                    if(currentNode_vertical->get_interval().end.getFraction() == 0.0 && currentNode_horizontal->get_interval().end.getFraction() == 0.0) {
                        currentNode_vertical->add_edge(currentNode_horizontal);
                    }
                }
            
                //iterate up until m-1 and add potential edges, break when empty or green gap
                for(PointID column = j; column < m-1; ++column) {
                    
                    if(matrix[column][i-1].vertical.empty()) break;
                                       
                    distance_t end_fraction = matrix[column][i-1].vertical[0].get_interval().end.getFraction();
                    
                    if(end_fraction != 0.0) {
                        currentNode_vertical->add_edge(&(matrix[column][i-1].vertical[0]));
                        break;
                    }               
                    if(matrix[column+1][i-1].vertical.empty()) {
                        break;
                    }

                    distance_t begin_fraction = matrix[column+1][i-1].vertical[0].get_interval().begin.getFraction();

                    if(begin_fraction == 0.0 && end_fraction == 0.0) {
                        currentNode_vertical->add_edge(&(matrix[column+1][i-1].vertical[0]));
                    }  
                    else break;           
                }

                //Check lower_right intervals 
                for(PointID row = i; row < n-1; ++row) {
                    if(matrix[j][row].vertical.size() > 0) {
                        size_t oneOrTwo = matrix[j][row].vertical.size() == 2 ? 1 : 0;
                        distance_t correct_fraction = currentNode_vertical->get_interval().end.getFraction() == 0.0 ? 1.0 : currentNode_vertical->get_interval().end.getFraction();
                        if(correct_fraction > matrix[j][row].vertical[oneOrTwo].get_interval().begin.getFraction()) {
                            currentNode_vertical->add_edge(&(matrix[j][row].vertical[oneOrTwo]));
                        }
                    }
                }
            }         
        }
    }

    //ADD dummy start and end node and connect to possible points on lower and right (upper, left)
    Freespace_Node s, e;

    for(PointID i = 0; i < n; ++i) {
        if(matrix[0][i].vertical.size() > 0) {
            if(matrix[0][i].vertical[0].get_interval().begin.getFraction() == 0.0)
                s.add_edge(&(matrix[0][i].vertical[0]));  
        }

        if(matrix[m-2][i].vertical.size() > 0) {
            size_t oneOrTwo = matrix[m-2][i].vertical.size() == 2 ? 1 : 0;
            if(matrix[m-2][i].vertical[oneOrTwo].get_interval().end.getFraction() == 0.0)
                matrix[m-2][i].vertical[oneOrTwo].add_edge(&e);
        }
    }

    for(PointID j = 0; j < m; ++j) {
        if(matrix[j][0].horizontal.size() > 0) {
            if(matrix[j][0].horizontal[0].get_interval().begin.getFraction() == 0.0)
                matrix[j][0].horizontal[0].add_edge(&e);
        }

        if(matrix[j][n-2].horizontal.size() > 0) {
            size_t oneOrTwo = matrix[j][n-2].horizontal.size() == 2 ? 1 : 0;
            if(matrix[j][n-2].horizontal[oneOrTwo].get_interval().end.getFraction() == 0.0)
                s.add_edge(&(matrix[j][n-2].horizontal[oneOrTwo]));
        }
    }

    std::queue<Freespace_Node*> q;
    q.push(&s);
    s.visited = true;
    int queue_count = 0;

    while(!q.empty() && !e.is_visited()) {
        queue_count++;
        Freespace_Node* current_node = q.front();
        q.pop(); 
        for(auto succ : current_node->get_successor_list()) {
            if(!succ->is_visited()) {

                succ->visited = true;
                q.push(succ);
                succ->steps = current_node->steps + 1;
                succ->parent = current_node;
            }
        }
    }

    std::vector<Freespace_Node*> path;

    if(!e.visited) {
        std::cout << "YES instance." << std::endl;
        return Certificate();
    }
    else {
        Freespace_Node* temp = &e;
        while(temp != nullptr) {
            path.push_back(temp);
            temp = temp->parent;
        }
    
        //reverse order and post process certificate
        std::reverse(path.begin(), path.end());
        std::vector<Freespace_Node*> final_path;
        final_path.push_back(path[1]);

        //due to using edges as nodes, it is importent to merge two nodes into the respective point that connects
        //two edges if they are connected
        Freespace_Node copy_node;
        Freespace_Node corner_node;
        for(size_t i = 1; i < path.size()-1; i++) {
            
            if(path[i-1]->c1 == path[i]->c1) {
                if(path[i]->c1 == path[i+1]->c1 + 1 && path[i]->c2 == path[i+1]->c2 - 1) {
                    final_path[final_path.size()-1]->c2 = final_path[final_path.size()-1]->c2 + 1;
                }
            }

            //2 edge run around the start node
            if(path[i]->c2 == 0) {
                if(path[i]->c1 == path[i+1]->c1 + 1 && path[i]->c2 == path[i+1]->c2 - 1) {
                    copy_node = *path[i];
                    copy_node.c2 += 1;
                    final_path.push_back(&copy_node);
                }
            }

            //2 edge run around the end node
            if(path[i]->c1 == curve1.size()-2 && path[i]->c2 == curve2.size()-2) {
                if(path[i]->c1 == path[i+1]->c1 && path[i]->c2 == path[i+1]->c2) {
                    copy_node = *path[i];
                    final_path.push_back(&copy_node);
                    path[i+1]->c2 += 1;
                }
            }
            final_path.push_back(path[i+1]);
        }

        //remove dummy end
        final_path.pop_back();

        //adjust coordinate for endpoint, since intervals end at point -1
        Freespace_Node final_adjustment_node;
        if(final_path[final_path.size()-2]->c2 == curve2.size()-2 && final_path[final_path.size()-1]->c2 == curve2.size()-2) {
            final_adjustment_node = *final_path[final_path.size()-1];
            final_path.push_back(&final_adjustment_node);
        }

        if(final_path[0]->c1 == curve1.size()-2 && !(final_path[0]->c2 == 0)) {
            final_path[0]->c1 = final_path[0]->c1 + 1;
        }
        if(final_path[final_path.size()-1]->c2 == curve2.size()-2 && !(final_path[final_path.size()-1]->c1 == 0)) {
            final_path[final_path.size()-1]->c2 = final_path[final_path.size()-1]->c2 + 1;
        }

        std::set<Freespace_Node*> duplicates;
        std::vector<Freespace_Node*> final_path_final;

        for(size_t i = 0; i < final_path.size(); i++) {
            duplicates.insert(final_path[i]);
        }
        for(size_t i = 0; i < final_path.size(); i++) {
            if(duplicates.find(final_path[i]) != duplicates.end()) {
                final_path_final.push_back(final_path[i]);
                duplicates.erase(final_path[i]);
            }
        }
        final_path = final_path_final;
        
        Certificate c;
        c.setCurves(&curve1, &curve2);
        c.setDistance(delta);

        for(size_t i = 0; i < final_path.size(); i++) {
            
            bool doublejump_flag = false;
            //final adjustment to interval borders and lower right jumps
            if(final_path[i]->get_interval().end.getFraction() != 0.0) final_path[i]->setEndFraction(final_path[i]->get_interval().end.getFraction()-0.00001);
            
            if(final_path[i]->get_interval().end.getFraction() == 0.0 && !final_path[i]->direction && (i != 0 && i != final_path.size()-1) && final_path[i+1]->c1 > final_path[i]->c1) {
                
                if(final_path[i]->c1 > final_path[i-1]->c1 && final_path[i]->c2 == final_path[i-1]->c2) {
                    CPosition cp = { CPoint(final_path[i]->c1, 0.0), CPoint(final_path[i]->c2, final_path[i]->get_interval().begin.getFraction()+0.0000001)};
                    c.addPoint(cp);
                   final_path[i]->c2 += 1; 
                   doublejump_flag = true;
                }
                else {
                    final_path[i]->c2 += 1; 
                }           
            } 

            auto frac = final_path[i]->get_interval().begin.getFraction() != 0.0 && !doublejump_flag ? final_path[i]->get_interval().begin.getFraction()+0.00000001 : final_path[i]->get_interval().end.getFraction();

            if(final_path[i]->direction) {
                CPosition cp = {CPoint(final_path[i]->c1, frac), CPoint(final_path[i]->c2, 0.0)};
                c.addPoint(cp);
            }
            else {
                CPosition cp = { CPoint(final_path[i]->c1, 0.0), CPoint(final_path[i]->c2, frac)};
                c.addPoint(cp);
            }
        }
        c.validate();
        c.setAnswer(false);
        c.dump_certificate();
        return c;
    }
}

Certificate ShortestCertificate::yes_certificate(Curve& curve1, Curve& curve2, Box_Set** matrix, double delta) {

    size_t n = curve1.size(), m = curve2.size();
    distance_t x = curve2[0].x - curve1[0].x, y = curve2[0].y - curve1[0].y;
    
    if(std::sqrt(x*x + y*y) > delta) {
        std::cout << "NO instance." << std::endl;
        return Certificate();
    }
    x = curve2[curve2.size()-1].x - curve1[curve1.size()-1].x;
    y = curve2[curve2.size()-1].y - curve1[curve1.size()-1].y;
    if(std::sqrt(x*x + y*y) > delta) {
        std::cout << "NO instance." << std::endl;
        return Certificate();
    }

    for(PointID j = 0; j < m-1; ++j) {  //Iterate each row

        std::set<distance_t> fractions;
        fractions.insert(0.0);
        PointID start = 0;
        PointID i = 0;

        while(i <= n) { //Iterate each element in row

            Interval interval_vertical = IntersectionAlgorithm::intersection_interval(curve1[i], delta, curve2[j], curve2[j+1]);

            if(!interval_vertical.is_empty() && !(i == n)) { //While not empty collect intervals
                fractions.insert(interval_vertical.begin);
                fractions.insert(interval_vertical.end);
            }
            else {
                fractions.insert(1.0);
                for(PointID section = start; section < i; ++section) { //if empty was found iterate from current start until i and create nodes

                    Interval interval_current = IntersectionAlgorithm::intersection_interval(curve1[section], delta, curve2[j], curve2[j+1]);
                    interval_current.begin += 0.00000001;
                    interval_current.end -= 0.00000001;

                    for(auto it = fractions.begin(); it != std::prev(fractions.end()); ++it) { //iterate over set and set fraction intervals
                       
                        Interval interval_section = Interval(*it, *std::next(it));                 
                        
                        if(interval_current.intersects(interval_section)) {
                            
                            CInterval c = CInterval(j, interval_section.begin, j, interval_section.end);
                            matrix[j][section].vertical.push_back(Freespace_Node(c, section, j, false));
                        }
                        else {
                            matrix[j][section].vertical.push_back(Freespace_Node());
                        } 
                    }
                }
                start = i+1;                
                fractions.clear();
            }
        ++i;      
        }
    }   

    for(PointID i = 0; i < n-1; ++i) {  //Iterate each column

        std::set<distance_t> fractions;
        fractions.insert(0.0);
        PointID start = 0;
        PointID j = 0;

        while(j <= m) { //Iterate each element in column

            Interval interval_horizontal = IntersectionAlgorithm::intersection_interval(curve2[j], delta, curve1[i], curve1[i+1]);

            if(!interval_horizontal.is_empty() && !(j == m)) { //While not empty collect intervals
                fractions.insert(interval_horizontal.begin);
                fractions.insert(interval_horizontal.end);
            }
            else {
                fractions.insert(1.0);
                for(PointID section = start; section < j; ++section) { //if empty was found iterate from current start until i and create nodes

                    Interval interval_current = IntersectionAlgorithm::intersection_interval(curve2[section], delta, curve1[i], curve1[i+1]);
                    interval_current.begin += 0.00000001;
                    interval_current.end -= 0.00000001;

                    for(auto it = fractions.begin(); it != std::prev(fractions.end()); ++it) { //iterate over set and set fraction intervals
                        
                        Interval interval_section = Interval(*it, *std::next(it));              
                        
                        if(interval_current.intersects(interval_section)) {
                            
                            CInterval c = CInterval(i, interval_section.begin, i, interval_section.end);
                            matrix[section][i].horizontal.push_back(Freespace_Node(c, i, section, true));
                        }
                        else {
                            matrix[section][i].horizontal.push_back(Freespace_Node());
                        } 
                    }
                }
                start = j+1;
                fractions.clear();
            }
        ++j;      
        }
    }
    matrix[m-1][n-1].horizontal.push_back(Freespace_Node(CInterval(curve2.size()-1, 0.0, curve2.size()-1, 0.0), curve1.size()-1, curve2.size()-1, 1));
    matrix[m-1][n-1].vertical.push_back(Freespace_Node(CInterval(curve1.size()-1, 0.0, curve1.size()-1, 0.0), curve1.size()-1, curve2.size()-1, 0));

    matrix[0][0].horizontal[0].visited = true;
    matrix[0][0].vertical[0].visited = true;

    for(PointID j = 0; j < m; ++j) {
        for(PointID i = 0; i < n; ++i) {
            //vertical connection to all possible intervals above

            //current Node in box
            for(size_t elem = 0; elem < matrix[j][i].vertical.size(); elem++) {
                
                Freespace_Node* current_node = &(matrix[j][i].vertical[elem]);
                if(!current_node->is_visited()) continue;

                size_t node_counter = elem + 1;
                PointID position_c1 = i;
                PointID position_c2 = j;

                //add top bound of box to edge set
                if(!(j == m-1) && !(i == n-1)) {
                    for(size_t elem_hor = 0; elem_hor < matrix[j+1][i].horizontal.size(); elem_hor++) {
                        if(!(matrix[j+1][i].horizontal[elem_hor].c1 == 1000000)) {
                            current_node->add_edge(&(matrix[j+1][i].horizontal[elem_hor]));
                            matrix[j+1][i].horizontal[elem_hor].visited = true;
                        }
                    }

                    for(size_t elem_ver = elem; elem_ver < matrix[j][i+1].vertical.size(); elem_ver++) {
                        if(!(matrix[j][i+1].vertical[elem_ver].c1 == 1000000)) {
                            current_node->add_edge(&(matrix[j][i+1].vertical[elem_ver]));
                            matrix[j][i+1].vertical[elem_ver].visited = true;
                        }
                    }
                }

                //add vertical nodes in x direction while proper
                while(position_c1 <= n) {

                    if(position_c1 == n || matrix[j][position_c1].vertical.empty() || matrix[j][position_c1].vertical[elem].c2 == 1000000) break;

                    current_node->add_edge(&(matrix[j][position_c1].vertical[elem]));
                    matrix[j][position_c1].vertical[elem].visited = true;
                    ++position_c1;
                }

                //add vertical nodes in y direction while proper
                while(position_c2 <= m) {
                    
                    if(node_counter == matrix[position_c2][i].vertical.size()) {
                        ++position_c2;
                        node_counter = 0;
                    }
                    if(position_c2 == m || matrix[position_c2][i].vertical.empty() || matrix[position_c2][i].vertical[node_counter].c2 == 1000000) break;

                    current_node->add_edge(&(matrix[position_c2][i].vertical[node_counter]));
                    matrix[position_c2][i].vertical[node_counter].visited = true;
                    ++node_counter;
                }
            }
        }
    }

    //compute for horizontal intervals
    for(PointID j = 0; j < m; ++j) {
        for(PointID i = 0; i < n; ++i) {
            //vertical connection to all possible intervals above

            //current Node in box
            for(size_t elem = 0; elem < matrix[j][i].horizontal.size(); elem++) {
                
                Freespace_Node* current_node = &(matrix[j][i].horizontal[elem]);
                if(!current_node->is_visited()) continue;

                size_t node_counter = elem + 1;
                PointID position_c1 = i;
                PointID position_c2 = j;

                //add right and top bound of box to edge set
                if(!(i == n-1) && !(j == m-1)) {
                    for(size_t elem_ver = 0; elem_ver < matrix[j][i+1].vertical.size(); elem_ver++) {
                        if(!(matrix[j][i+1].vertical[elem_ver].c1 == 1000000)) {
                            current_node->add_edge(&(matrix[j][i+1].vertical[elem_ver]));
                            matrix[j][i+1].vertical[elem_ver].visited = true;
                        }
                    }

                    for(size_t elem_hor = elem; elem_hor < matrix[j+1][i].horizontal.size(); elem_hor++) {
                        if(!(matrix[j+1][i].horizontal[elem_hor].c1 == 1000000)) {
                            current_node->add_edge(&(matrix[j+1][i].horizontal[elem_hor]));
                            matrix[j+1][i].horizontal[elem_hor].visited = true;
                        }
                    }
                }

                //add horizontal nodes in y direction while proper
                while(position_c2 <= m) {
                    if(position_c2 == m || matrix[position_c2][i].horizontal.empty() || matrix[position_c2][i].horizontal[elem].c2 == 1000000) break;
                    current_node->add_edge(&(matrix[position_c2][i].horizontal[elem]));
                    matrix[position_c2][i].horizontal[elem].visited = true;
                    ++position_c2;
                }

                //add horizontal nodes in x direction while proper
                while(position_c1 <= n) {
                    
                    if(node_counter == matrix[j][position_c1].horizontal.size()) {
                        ++position_c1;
                        node_counter = 0;
                    }
                    if(position_c1 == n || matrix[j][position_c1].horizontal.empty() || matrix[j][position_c1].horizontal[node_counter].c2 == 1000000) break;

                    current_node->add_edge(&(matrix[j][position_c1].horizontal[node_counter]));
                    matrix[j][position_c1].horizontal[node_counter].visited = true;
                    ++node_counter;
                }
            }
        }
    }

    //HEURISTIC ADDITION 
    /*auto adj_matrix = find_jumps(curve1, curve2, delta);

    for(PointID j = 0; j < m-1; ++j) {
        for(PointID i = 0; i < n-1; ++i) {
            for(auto elem : adj_matrix[j][i]) {
                //std::cout << j << " " << i << " | " << elem.first << " " << elem.second << std::endl;
                if(!(elem.first == m-1) && !(elem.second == n-1)) matrix[j][i].horizontal[0].add_edge(&matrix[elem.first][elem.second].horizontal[0]);
                if(!(elem.first == m-1) && !(elem.second == n-1)) matrix[j][i].vertical[0].add_edge(&matrix[elem.first][elem.second].vertical[0]);
            }
        }
    }*/

    Freespace_Node s, e;
    e.visited = true;

    s.add_edge(&(matrix[0][0].vertical[0]));
    s.add_edge(&(matrix[0][0].horizontal[0]));
    if(!matrix[m-1][n-1].horizontal.empty()) matrix[m-1][n-1].horizontal[0].add_edge(&e);
    if(!matrix[m-1][n-1].vertical.empty()) matrix[m-1][n-1].vertical[0].add_edge(&e);

    std::queue<Freespace_Node*> q;
    q.push(&s);
    s.visited = false;
    int queue_count = 0;

    while(!q.empty() && e.is_visited()) {
        queue_count++;
        Freespace_Node* current_node = q.front();
        q.pop();             

        for(auto succ : current_node->get_successor_list()) {
            if(succ->is_visited()) {

                succ->visited = false;
                q.push(succ);
                succ->steps = current_node->steps + 1;
                succ->parent = current_node;
            }
        }
    }
    std::vector<Freespace_Node*> path;

    if(e.visited) {
        std::cout << "NO instance." << std::endl;
        return Certificate();
    }
    else {
        Freespace_Node* temp = &e;
        while(temp != nullptr) {
            path.push_back(temp);
            temp = temp->parent;
        }

        path.pop_back();
        std::reverse(path.begin(), path.end());
        path.pop_back();

        path[path.size()-1]->get_interval().begin.setFraction(0.0);
        path[path.size()-1]->get_interval().end.setFraction(0.0);
    
        Certificate c;
        c.setCurves(&curve1, &curve2);
        c.setDistance(delta);

        for(size_t i = 0; i < path.size()-1; i++) {
            
            if(path[i]->direction) {
                CPosition cp = {CPoint(path[i]->c1, path[i]->get_interval().begin.getFraction()), CPoint(path[i]->c2, 0.0)};
                c.addPoint(cp);
            }
            else {
                CPosition cp = {CPoint(path[i]->c1, 0.0), CPoint(path[i]->c2, path[i]->get_interval().begin.getFraction())};
                c.addPoint(cp);
            }
        }

        c.addPoint({CPoint(path[path.size()-1]->c1, 0.0), CPoint(path[path.size()-1]->c2, 0.0)});
    
        c.validate();
        c.setAnswer(true);
        c.dump_certificate();

        //DIAGONAL ADDITION
        /*Certificate alt = pruning_diagonals(c, curve1, curve2, delta);
        std::cout << "\nHEURISTIC" << std::endl;
        alt.dump_certificate();*/

        return c;
    }
    return Certificate();
}
