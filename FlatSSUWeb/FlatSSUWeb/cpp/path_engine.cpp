// ============================================================
//  path_engine.cpp  —  FlatSSU Path Engine (annotated)
//
//  캠퍼스 내 지점을 노드로, 노드 간 연결 간선을 정의하여
//  다익스트라 알고리즘을 통해 최단/편한 경로를 계산한다.
//  WebAssembly (emscripten) 및 Windows CLI 양쪽에서 동작하도록 설계되었습니다.
// ============================================================


#include <cstdio>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <queue>
#include <algorithm>
#include <limits>
#include <cctype>       // isdigit
#include <filesystem>   // 디버깅용 현재 작업 디렉토리 출력

#include "json.hpp"     // nlohmann::json 헤더

// ─── emscripten 헤더 (웹 빌드 시에만) ──────────────────────────
#ifdef __EMSCRIPTEN__
#include <emscripten/bind.h>
namespace em = emscripten;           // 접두어 짧게 별칭

#endif
using namespace std;
using json = nlohmann::json;

// ─────────────────────────────────────────────────────────────
//  자료 구조
// ─────────────────────────────────────────────────────────────
// Node: 각 지점(건물 입구 등)의 이름과 위도/경도 저장
struct Node { string name; double lat, lng; };
struct Edge { int to, dist, conv; };          // conv 0-600
using Graph = vector<vector<Edge>>;

// 전역 변수: 노드 목록, 이름->인덱스 맵, 그래프
static vector<Node>             nodes;
static unordered_map<string, int> name2id;
static Graph                    graph;

// ─────────────────────────────────────────────────────────────
//  비용 함수 파라미터
// ─────────────────────────────────────────────────────────────
enum Mode { SHORTEST, CONVENIENT };
static constexpr int MAX_CONV = 600;   // conv 상한
static constexpr int BETA = 5;    // 패널티 가중치

inline int edgeCost(const Edge& e, Mode m) {
    if (m == SHORTEST) {
        // 최단 경로: 단순 거리 사용
        return e.dist;
    }
    else {
        // 편한 경로: 거리 + (MAX_CONV - conv) * BETA
        // conv가 높을수록 패널티가 작아짐
        return e.dist + (MAX_CONV - e.conv) * BETA;
    }
}
// ─────────────────────────────────────────────────────────────
//  다익스트라 + 경로 복원
// ─────────────────────────────────────────────────────────────
vector<int> dijkstra(const Graph& G, int S, int T, Mode m)
{
    const int INF = numeric_limits<int>::max();
    int N = static_cast<int>(G.size());
    vector<int> dist(N, INF), prev(N, -1);
    using P = pair<int, int>; priority_queue<P, vector<P>, greater<P>> pq;
    dist[S] = 0; pq.push({ 0,S });

    while (!pq.empty()) {
        auto [cd, u] = pq.top(); pq.pop();
        if (cd > dist[u]) continue;
        if (u == T) break;
        for (const auto& e : G[u]) {
            int nd = cd + edgeCost(e, m);
            if (nd < dist[e.to]) { dist[e.to] = nd; prev[e.to] = u; pq.push({ nd,e.to }); }
        }
    }
    vector<int> path;
    if (dist[T] < INF) { for (int v = T; v != -1; v = prev[v]) path.push_back(v); reverse(path.begin(), path.end()); }
    return path;
}

// ─────────────────────────────────────────────────────────────
//  문자열 유틸 : trim & 숫자 판별
// ─────────────────────────────────────────────────────────────
static inline void trim(string& s)
{
    const char* ws = " \t\r\n";
    s.erase(0, s.find_first_not_of(ws));
    s.erase(s.find_last_not_of(ws) + 1);
}
static inline bool isNumber(const string& s)
{
    return !s.empty() &&
        all_of(s.begin(), s.end(),
            [](unsigned char c) { return isdigit(c); });
}

// ─────────────────────────────────────────────────────────────
//  그래프 초기화
// ─────────────────────────────────────────────────────────────
void initGraph()
{
    cout << "[DEBUG] CWD = " << filesystem::current_path() << '\n';

    // 1) 노드
    ifstream fn("flatssu_nodes.json");
    if (!fn) { cerr << "[ERR] flatssu_nodes.json 열기 실패\n"; return; }
    json jn; fn >> jn; fn.close();
    nodes.clear(); name2id.clear();
    for (auto& it : jn) {
        string nm = it["name"].get<string>();
        double la = it["lat"].get<double>();
        double ln = it["lng"].get<double>();
        name2id[nm] = static_cast<int>(nodes.size());
        nodes.push_back({ nm,la,ln });
    }

    // 2) 간선
    ifstream fe("node_connect.csv");
    if (!fe) { cerr << "[ERR] node_connect.csv 열기 실패\n"; return; }
    string line; getline(fe, line);       // 헤더 skip
    graph.assign(nodes.size(), {});
    int row = 1, bad = 0;
    while (getline(fe, line)) {
        ++row;
        string a, b, sd, sc; stringstream ss(line);
        getline(ss, a, ','); getline(ss, b, ','); getline(ss, sd, ','); getline(ss, sc, ',');
        trim(a); trim(b); trim(sd); trim(sc);
        if (a.empty() || b.empty() || !isNumber(sd) || !isNumber(sc)) {
            cerr << "[WARN] bad row " << row << " : " << line << '\n'; ++bad; continue;
        }
        int u = name2id[a], v = name2id[b];
        int d = stoi(sd), cv = stoi(sc);
        cv = min(max(cv, 0), MAX_CONV);        // 범위 클램핑
        graph[u].push_back({ v,d,cv });
        graph[v].push_back({ u,d,cv });
    }
    size_t ec = 0; for (auto& vec : graph) ec += vec.size();
    cout << "[INFO] initGraph : 노드 " << nodes.size() << ", 간선 " << ec / 2
        << "개, 무시 " << bad << "행\n";
}

// ─────────────────────────────────────────────────────────────
//  WebAssembly 바인딩
// ─────────────────────────────────────────────────────────────
vector<int> findShortest(int s, int t) { return dijkstra(graph, s, t, SHORTEST); }
vector<int> findConvenient(int s, int t) { return dijkstra(graph, s, t, CONVENIENT); }

#ifdef __EMSCRIPTEN__
EMSCRIPTEN_BINDINGS(FlatSSUPath) {
    em::function("initGraph", &initGraph);
    em::function("findShortest", &findShortest);
    em::function("findConvenient", &findConvenient);
    em::register_vector<int>("VectorInt");
}
#endif

// ─────────────────────────────────────────────────────────────
//  CLI 테스트 (Windows 실행)
// ─────────────────────────────────────────────────────────────
#ifndef __EMSCRIPTEN__
int main()
{
    initGraph();
    if (graph.empty()) return 1;
    int s, t; cout << "start idx : "; cin >> s; cout << "end idx : "; cin >> t;
    auto p = dijkstra(graph, s, t, CONVENIENT);
    for (size_t i = 0; i < p.size(); ++i) {
        cout << nodes[p[i]].name << (i + 1 < p.size() ? " -> " : "\n");
    }
}
#endif
