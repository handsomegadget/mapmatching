#include <bits/stdc++.h>
using namespace std;
typedef long double ld;
typedef std::pair<int, int> pii;
const ld c1 = 111194.92664; //上海地区单位度数纬度对应的米
const ld c2 = 95112.16616; //上海地区单位经度对应的米
ld MIN_x = 130, MAX_x = 0, MIN_y = 40, MAX_y = 0.0;

struct pt {
    ld x, y;
};
int total_roads;

class Road {
public:
    Road() {}

    std::vector<pt> RoadPosition;
    int ID;
    int start_pt;
    int terminal_pt;
    int RoadLevel;
    int pt_num = 0;
    ld RoadLength = 0;
};

class Edges {
public:
    void AddIn(Road &e) { edges.push_back(e); }

    std::vector<Road> edges;
};

std::set<int> node[60010]; //记录每个节点出发的路的编号
class Grid {
public:
    Grid() {}
    std::map<pii, std::set<int>> Grids;
};

class Candidate_Points_Set {
public:
    int quantity = 0;
    std::vector<pt> candidate_points;
    std::vector<ld> distance_to_begin;
    std::vector<int> RoadID;
    std::vector<ld> OBP; //观测概率
    std::vector<std::vector<ld>> TRP; //传递概率
};

std::vector<Candidate_Points_Set> candidates;
struct trajectory_element {
    int timestamp;
    pt point;
} trajectory[2300];
Edges E;
Grid Cells;
ld smaller(ld a, ld b){
    if(a>b)return b;
    else return a;
}
ld distance_of_points(pt x0, pt y0) {
    return sqrt((x0.x - y0.x) * (x0.x - y0.x) + (x0.y - y0.y) * (x0.y - y0.y));
}

ld distance_of_point_line(pt p, pt j, pt k) {
    ld x1 = j.x;
    ld x2 = k.x;
    ld y1 = j.y;
    ld y2 = k.y;
    ld x = p.x;
    ld y = p.y;
    ld slope = (y1 - y2) / (x1 - x2);
    ld intercept = (x2 * y1 - x1 * y2) / (x2 - x1);
    if (((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1)) >= 0 && ((x1 - x2) * (x - x2) + (y1 - y2) * (y - y2)) >= 0)
        return (abs(-slope * x - intercept + y)) / sqrt(1 + (slope) * (slope));
    else if ((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1) < 0)
        return sqrt((x - x1) * (x - x1) + (y - y1) * (y - y1));
    else return sqrt((x - x2) * (x - x2) + (y - y2) * (y - y2));
}

ld CalCulate_Road_Length(Road &r) {
    ld res = 0;
    for (int i = 0; i < r.pt_num - 1; i++) {
        res += distance_of_points(r.RoadPosition[i], r.RoadPosition[i + 1]);
    }
    return res;
}

ld calculate_RoadStart_to_Candidate(Road &r, int min_pos, pt candi) {
    ld res = 0;
    int i;
    for (i = 0; i < min_pos; i++) {
        res += distance_of_points(r.RoadPosition[i], r.RoadPosition[i + 1]);
    }
    res += distance_of_points(r.RoadPosition[i], candi);
    return res;
}

void FindNearestPoint(pt p, pt j, pt k, pt &ans) {
    ld x1 = j.x;
    ld x2 = k.x;
    ld y1 = j.y;
    ld y2 = k.y;
    ld x = p.x;
    ld y = p.y;
    if (((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1)) >= 0 && ((x1 - x2) * (x - x2) + (y1 - y2) * (y - y2)) >= 0) {
        ld A1 = y2 - y1;
        ld B1 = x1 - x2;
        ld C1 = y1 * (x1 - x2) - x1 * (y1 - y2);
        ld A2 = x1 - x2;
        ld B2 = y1 - y2;
        ld C2 = (x1 - x2) * x + (y1 - y2) * y;
        ans.x = (C1 * B2 - C2 * B1) / (A1 * B2 - A2 * B1);
        ans.y = (A1 * C2 - A2 * C1) / (A1 * B2 - A2 * B1);
    } else if ((x2 - x1) * (x - x1) + (y2 - y1) * (y - y1) < 0) {
        ans.x = x1;
        ans.y = y1;
    } else {
        ans.x = x2;
        ans.y = y2;
    }
}

void FindCandidatesNearby(int seq, Candidate_Points_Set &x) {
    pt p = trajectory[seq].point;
    int m = (int) floorl((MAX_y - p.y) / 50);
    int n = (int) floorl((p.x - MIN_x) / 50);
    pii location = std::make_pair(m, n);
    std::set<int> road_overlap;
    std::set<int> node_overlap;
    for (int i = location.first - 1; i <= location.first + 1; i++) {
        for (int j = location.second - 1; j <= location.second + 1; j++) {
            auto it = Cells.Grids.find(std::make_pair(i, j));
            if (it != Cells.Grids.end()) {
                for (auto k = it->second.begin(); k != it->second.end(); k++) {
                    Road temp = E.edges[*k];
                    if (road_overlap.find(temp.ID) != road_overlap.end()) continue;
                    road_overlap.insert(temp.ID);
                    pt res_candi;
                    ld min_dis = 50;
                    int min_pos = -1;
                    for (int d = 0; d < temp.pt_num - 1; d++) {
                        ld dis = distance_of_point_line(p, temp.RoadPosition[d], temp.RoadPosition[d + 1]);
                        if (dis < min_dis) {
                            min_pos = d;
                            min_dis = dis;
                        }
                    }
                    if (min_dis < 49) {
                        FindNearestPoint(p, temp.RoadPosition[min_pos], temp.RoadPosition[min_pos + 1], res_candi);
                        if (distance_of_points(res_candi, temp.RoadPosition[0]) < 1e-10) {
                            int l = temp.start_pt;
                            if (node_overlap.find(l) != node_overlap.end()) continue;
                            node_overlap.insert(l);
                        }
                        if (distance_of_points(res_candi, temp.RoadPosition[temp.pt_num - 1]) < 1e-10) {
                            int l = temp.terminal_pt;
                            if (node_overlap.find(l) != node_overlap.end()) continue;
                            node_overlap.insert(l);
                        }
                        x.quantity++;
                        x.candidate_points.push_back(res_candi);
                        x.distance_to_begin.push_back(calculate_RoadStart_to_Candidate(temp, min_pos, res_candi));
                        x.RoadID.push_back(temp.ID);
                        x.OBP.push_back(pow(2, -0.00125 * min_dis * min_dis));
                    }
                }
            }
        }
    }
    if (x.quantity == 0) {
        road_overlap.clear();
        node_overlap.clear();
        for (int i = location.first - 2; i <= location.first + 2; i++) {
            for (int j = location.second - 2; j <= location.second + 2; j++) {
                auto it = Cells.Grids.find(std::make_pair(i, j));
                if (it != Cells.Grids.end()) {
                    for (auto k = it->second.begin(); k != it->second.end(); k++) {
                        Road temp = E.edges[*k];
                        if (road_overlap.find(temp.ID) != road_overlap.end()) continue;
                        road_overlap.insert(temp.ID);
                        pt res_candi;
                        ld min_dis = 75;
                        int min_pos = -1;
                        for (int d = 0; d < temp.pt_num - 1; d++) {
                            ld dis = distance_of_point_line(p, temp.RoadPosition[d], temp.RoadPosition[d + 1]);
                            if (dis < min_dis) {
                                min_pos = d;
                                min_dis = dis;
                            }
                        }
                        if (min_dis < 74) {
                            FindNearestPoint(p, temp.RoadPosition[min_pos], temp.RoadPosition[min_pos + 1], res_candi);
                            if (distance_of_points(res_candi, temp.RoadPosition[0]) < 1e-10) {
                                int l = temp.start_pt;
                                if (node_overlap.find(l) != node_overlap.end()) continue;
                                node_overlap.insert(l);
                            }
                            if (distance_of_points(res_candi, temp.RoadPosition[temp.pt_num - 1]) < 1e-10) {
                                int l = temp.terminal_pt;
                                if (node_overlap.find(l) != node_overlap.end()) continue;
                                node_overlap.insert(l);
                            }
                            x.quantity++;
                            x.candidate_points.push_back(res_candi);
                            x.distance_to_begin.push_back(calculate_RoadStart_to_Candidate(temp, min_pos, res_candi));
                            x.RoadID.push_back(temp.ID);
                            x.OBP.push_back(pow(2, -0.00125 * min_dis * min_dis));
                        }
                    }
                }
            }
        }
    }
}

void FindFirstCandidate(Candidate_Points_Set &x) {
    pt p = trajectory[0].point;
    int m = (int) floorl((MAX_y - p.y) / 50);
    int n = (int) floorl((p.x - MIN_x) / 50);
    pii location = std::make_pair(m, n);
    std::set<int> road_overlap;
    std::set<int> node_overlap;
    for (int i = location.first - 1; i <= location.first + 1; i++) {
        for (int j = location.second - 1; j <= location.second + 1; j++) {
            auto it = Cells.Grids.find(std::make_pair(i, j));
            if (it != Cells.Grids.end()) {
                for (auto k = it->second.begin(); k != it->second.end(); k++) {
                    Road temp = E.edges[*k];
                    if (road_overlap.find(temp.ID) != road_overlap.end()) continue;
                    road_overlap.insert(temp.ID);
                    pt res_candi;
                    ld min_dis = 50;
                    int min_pos = -1;
                    for (int d = 0; d < temp.pt_num - 1; d++) {
                        ld dis = distance_of_point_line(p, temp.RoadPosition[d], temp.RoadPosition[d + 1]);
                        if (dis < min_dis) {
                            min_pos = d;
                            min_dis = dis;
                        }
                    }
                    if (min_dis < 49) {
                        FindNearestPoint(p, temp.RoadPosition[min_pos], temp.RoadPosition[min_pos + 1], res_candi);
                        if (distance_of_points(res_candi, temp.RoadPosition[0]) < 1e-10) {
                            int l = temp.start_pt;
                            if (node_overlap.find(l) != node_overlap.end()) continue;
                            node_overlap.insert(l);
                        }
                        if (distance_of_points(res_candi, temp.RoadPosition[temp.pt_num - 1]) < 1e-10) {
                            int l = temp.terminal_pt;
                            if (node_overlap.find(l) != node_overlap.end()) continue;
                            node_overlap.insert(l);
                        }
                        x.quantity++;
                        x.candidate_points.push_back(res_candi);
                        x.distance_to_begin.push_back(calculate_RoadStart_to_Candidate(temp, min_pos, res_candi));
                        x.RoadID.push_back(temp.ID);
                        x.OBP.push_back(pow(2, -0.00125 * min_dis * min_dis));
                    }
                }
            }
        }
    }
    if (x.quantity == 0) {
        road_overlap.clear();
        node_overlap.clear();
        for (int i = location.first - 2; i <= location.first + 2; i++) {
            for (int j = location.second - 2; j <= location.second + 2; j++) {
                auto it = Cells.Grids.find(std::make_pair(i, j));
                if (it != Cells.Grids.end()) {
                    for (auto k = it->second.begin(); k != it->second.end(); k++) {
                        Road temp = E.edges[*k];
                        if (road_overlap.find(temp.ID) != road_overlap.end()) continue;
                        road_overlap.insert(temp.ID);
                        pt res_candi;
                        ld min_dis = 75;
                        int min_pos = -1;
                        for (int d = 0; d < temp.pt_num - 1; d++) {
                            ld dis = distance_of_point_line(p, temp.RoadPosition[d], temp.RoadPosition[d + 1]);
                            if (dis < min_dis) {
                                min_pos = d;
                                min_dis = dis;
                            }
                        }
                        if (min_dis < 74) {
                            FindNearestPoint(p, temp.RoadPosition[min_pos], temp.RoadPosition[min_pos + 1], res_candi);
                            if (distance_of_points(res_candi, temp.RoadPosition[0]) < 1e-10) {
                                int l = temp.start_pt;
                                if (node_overlap.find(l) != node_overlap.end()) continue;
                                node_overlap.insert(l);
                            }
                            if (distance_of_points(res_candi, temp.RoadPosition[temp.pt_num - 1]) < 1e-10) {
                                int l = temp.terminal_pt;
                                if (node_overlap.find(l) != node_overlap.end()) continue;
                                node_overlap.insert(l);
                            }
                            x.quantity++;
                            x.candidate_points.push_back(res_candi);
                            x.distance_to_begin.push_back(calculate_RoadStart_to_Candidate(temp, min_pos, res_candi));
                            x.RoadID.push_back(temp.ID);
                            x.OBP.push_back(pow(2, -0.00125 * min_dis * min_dis));
                        }
                    }
                }
            }
        }
    }
    assert(x.quantity > 0);
}

void CalculateTRP(int num) {
    for (int cnt = 1; cnt < num; cnt++) {
        for (int i = 0; i < candidates[cnt].quantity; i++) {
            std::vector<ld> trp;
            std::vector<ld> tmp;
            for (int j = 0; j < candidates[cnt - 1].quantity; j++) {
                ld DisOfGPSPoint = distance_of_points(candidates[cnt].candidate_points[i],
                                                      candidates[cnt - 1].candidate_points[j]);
                ld ShortestRoute = 99999;
                int _start = candidates[cnt - 1].RoadID[j];
                int _end = candidates[cnt].RoadID[i];
                if (_start == _end) {
                    ShortestRoute = abs(candidates[cnt].distance_to_begin[i] - candidates[cnt - 1].distance_to_begin[j]);

                } else {
                    for (auto i1 = node[E.edges[_start].terminal_pt].begin();
                         i1 != node[E.edges[_start].terminal_pt].end(); i1++) {

                        if (*i1 == _end) {
                            ld temp = E.edges[_start].RoadLength + candidates[cnt].distance_to_begin[i] -
                                      candidates[cnt - 1].distance_to_begin[j];
                            ShortestRoute = smaller(abs(temp), ShortestRoute);
                        } else {
                            for (auto i2 = node[E.edges[*i1].terminal_pt].begin();
                                 i2 != node[E.edges[*i1].terminal_pt].end(); i2++) {
                                if (*i2 == _end) {
                                    ld temp = E.edges[_start].RoadLength + E.edges[*i1].RoadLength +
                                              candidates[cnt].distance_to_begin[i] -
                                              candidates[cnt - 1].distance_to_begin[j];
                                    ShortestRoute = smaller(abs(temp), ShortestRoute);
                                } else {
                                    for (auto i3 = node[E.edges[*i2].terminal_pt].begin();
                                         i3 != node[E.edges[*i2].terminal_pt].end(); i3++) {
                                        if (*i3 == _end) {
                                            ld temp = E.edges[_start].RoadLength + E.edges[*i2].RoadLength +
                                                      E.edges[*i1].RoadLength + candidates[cnt].distance_to_begin[i] -
                                                      candidates[cnt - 1].distance_to_begin[j];
                                            ShortestRoute = smaller(abs(temp), ShortestRoute);
                                        }

                                    }
                                }
                            }
                        }
                    }
                }
                ld c, d;
                ld AverageSpeed, ExpectedSpeed;
                if(DisOfGPSPoint == 0){
                    c = 1;
                    //ExpectedSpeed = speed[E.edges[_end].RoadLevel];
                    //AverageSpeed = ExpectedSpeed;
                }
                else {
                    //ExpectedSpeed = (speed[E.edges[_start].RoadLevel]+speed[E.edges[_end].RoadLevel])/2.0;
                    //AverageSpeed = ShortestRoute/(trajectory[cnt].timestamp - trajectory[cnt-1].timestamp);
                    c = DisOfGPSPoint/ShortestRoute;}
                //d = pow(2, -0.05*(AverageSpeed-ExpectedSpeed)*(AverageSpeed-ExpectedSpeed));
                trp.push_back(c);
                //tmp.push_back(d);
            }
            candidates[cnt].TRP.push_back(trp);
            //candidates[cnt].TMP.push_back(tmp);
        }
    }
}

void ReadInRoad() {
    scanf("%d", &total_roads);
    for (int i = 0; i < total_roads; i++) {
        Road r;
        char Good_For_Nothing[100];
        scanf("%d", &r.ID);
        scanf("%d", &r.start_pt);
        scanf("%d", &r.terminal_pt);
        cin >> Good_For_Nothing;
        scanf("%d", &r.RoadLevel);
        scanf("%d", &r.pt_num);
        for (int j = 0; j < r.pt_num; j++) {
            pt p;
            scanf("%Lf%Lf", &p.y, &p.x);
            MIN_x = std::min(MIN_x, p.x);
            MAX_y = std::max(MAX_y, p.y);
            p.y *= c1;
            p.x *= c2;
            r.RoadPosition.push_back(p);
        }
        r.RoadLength = CalCulate_Road_Length(r);
        E.AddIn(r);
        node[r.start_pt].insert(r.ID);
    }
}

int readTrajectory() {
    memset(trajectory, 0, 2300 * sizeof(trajectory_element));
    int temp;
    int i = 0;
    while (true) {
        scanf("%d", &temp);
        if (temp < 10000) break;
        trajectory[i].timestamp = temp;
        ld y, t;
        scanf("%Lf%Lf", &y, &t);
        trajectory[i].point.x = t * c2;
        trajectory[i].point.y = y * c1;
        ++i;
    }
    return i;
}

void InitializeGrid() {
    MIN_x *= c2;
    MAX_y *= c1;
    MIN_x -= 60;
    MAX_y += 60;
    for (int i = 0; i < total_roads; i++) {
        for (int j = 0; j < E.edges[i].pt_num; j++) {
            int m = (int) floorl((MAX_y - E.edges[i].RoadPosition[j].y) / 50);
            int n = (int) floorl((E.edges[i].RoadPosition[j].x - MIN_x) / 50);
            Cells.Grids[std::make_pair(m, n)].insert(i);
        }
        for (int j = 0; j < E.edges[i].pt_num - 1; j++) {
            int m = (int) floorl(
                    (MAX_y - (E.edges[i].RoadPosition[j].y + 3 * E.edges[i].RoadPosition[j + 1].y) / 4) / 50);
            int n = (int) floorl(
                    ((E.edges[i].RoadPosition[j].x + 3 * E.edges[i].RoadPosition[j + 1].x) / 4 - MIN_x) / 50);
            Cells.Grids[std::make_pair(m, n)].insert(i);
        }
        for (int j = 0; j < E.edges[i].pt_num - 1; j++) {
            int m = (int) floorl((MAX_y - (E.edges[i].RoadPosition[j].y + E.edges[i].RoadPosition[j + 1].y) / 2) / 50);
            int n = (int) floorl(((E.edges[i].RoadPosition[j].x + E.edges[i].RoadPosition[j + 1].x) / 2 - MIN_x) / 50);
            Cells.Grids[std::make_pair(m, n)].insert(i);
        }
        for (int j = 0; j < E.edges[i].pt_num - 1; j++) {
            int m = (int) floorl(
                    (MAX_y - (3 * E.edges[i].RoadPosition[j].y + E.edges[i].RoadPosition[j + 1].y) / 4) / 50);
            int n = (int) floorl(
                    ((3 * E.edges[i].RoadPosition[j].x + E.edges[i].RoadPosition[j + 1].x) / 4 - MIN_x) / 50);
            Cells.Grids[std::make_pair(m, n)].insert(i);
        }
    }
}

void HMM_Matching(int num){
    ld f[2300][28] = {{0}};
    int pre[2300][28] = {{-1}};
    ld max = -1;
    for(int p = 0; p < candidates[0].quantity; p++)
        f[0][p] = candidates[0].OBP[p];
    for(int i = 1; i < num; i++){
        for(int j = 0; j < candidates[i].quantity; j++){
            max = -9999;
            for(int k = 0; k < candidates[i-1].quantity; k++){
                ld alt = f[i-1][k]*candidates[i].OBP[j]*candidates[i].TRP[j][k];
                if(alt > max){
                    max = alt;
                    pre[i][j] = k;
                }
                f[i][j] = max;
            }
        }
    }
    std::stack<int> matched_points;
    ld MAX = -1;
    int c = -1;
    for(int i = 0; i < candidates[num-1].quantity; i++){
        if(f[num-1][i] > MAX){
            MAX = f[num-1][i];
            c = i;
        }
    }
    for(int i = num-1; i > 0; i--){
        matched_points.push(c);
        c = pre[i][c];
    }
    matched_points.push(c);
    for(int i = 0; i < num; i++){
        int res = matched_points.top();
        printf("%d ", candidates[i].RoadID[res]);
        //cout << candidates[i].RoadID[res] << ' ';
        matched_points.pop();
    }
    cout << '\n';
}

int main() {
    ReadInRoad();
    InitializeGrid();
    int test_rounds;
    cin >> test_rounds;
    cout << test_rounds << std::endl;
    while (test_rounds--) {
        int length_of_trajectory = readTrajectory();
        Candidate_Points_Set candi;
        FindFirstCandidate(candi);
        candidates.push_back(candi);
        for (int i = 1; i < length_of_trajectory; i++) {
            Candidate_Points_Set _candi;
            FindCandidatesNearby(i, _candi);
            candidates.push_back(_candi);
        }
        CalculateTRP(length_of_trajectory);
        HMM_Matching(length_of_trajectory);
        std::vector<Candidate_Points_Set>().swap(candidates);
    }
    return 0;
}



