#pragma once

#include <iostream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <limits.h>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cstdint>
#include <functional>

// ─────────────────────── БАЗОВЫЕ СТРУКТУРЫ ───────────────────────
struct Edge  { int target, weight; };
struct Node  { int v, g, f, t;  bool operator>(const Node& o) const { return f > o.f; } };
struct Agent {
    int id;
    std::vector<int> route;          // шахта‑терминал‑…
    int idx = 0;
    bool finished = false;
    int start() const { return route[idx]; }
    int goal()  const { return route[(idx+1) % route.size()]; }
    void advance(){ idx = (idx+1) % route.size(); }
};

std::vector<int> dijkstra(const std::vector<std::vector<Edge>>& g, int s)
{
    int n = g.size();
    std::vector<int> d(n, INT_MAX); d[s] = 0;
    std::priority_queue<std::pair<int,int>,
                        std::vector<std::pair<int,int>>,
                        std::greater<>> pq;
    pq.push({0, s});
    while (!pq.empty()) {
        auto [dist,u] = pq.top(); pq.pop();
        if (dist>d[u]) continue;
        for (auto& e: g[u])
            if (d[u]+e.weight < d[e.target]) {
                d[e.target]=d[u]+e.weight;
                pq.push({d[e.target], e.target});
            }
    }
    return d;
}

using Reserve = std::unordered_map<int, std::unordered_set<int>>;  // time ➜ {vertices}

std::vector<int> a_star_with_reserve(const std::vector<std::vector<Edge>>& g,
                                     int start, int goal,
                                     const std::vector<int>& h,
                                     const Reserve& res,     // global reserve
                                     int tOffset)
{
    struct OpenNode {
        int v;           // вершина
        int g;           // фактическая стоимость от старта
        int t;           // время (tick) прихода в v
        int f;           // g + h
        bool operator>(const OpenNode& o) const { return f > o.f; }
    };

    std::priority_queue<OpenNode,
                        std::vector<OpenNode>,
                        std::greater<OpenNode>> open;
    open.push({start, 0, 0, h[start]});

    // key = (vertex << 20) ^ time   (20 бит достаточно до 1 M tick’ов)
    auto key = [](int v, int t){
        return (static_cast<std::uint64_t>(v) << 20) ^ t;
    };
    std::unordered_map<std::uint64_t, int> bestG;   // key ➜ g‑cost
    bestG[key(start,0)] = 0;

    // проверка: занята ли вершина v в глобальный момент (t + tOffset)
    auto busy_vertex = [&](int t, int v) -> bool
    {
        auto it = res.find(t + tOffset);          // переводим во «внешние» тики
        return it != res.end() && it->second.count(v);
    };

    // для восстановления пути
    struct Back { int prev_v, prev_t; };
    std::unordered_map<std::uint64_t, Back> back;

    while(!open.empty()){
        auto cur = open.top(); open.pop();

        // если пришли в цель и в резерве свободно
        if(cur.v == goal){
            std::vector<int> path;
            int cv = cur.v, ct = cur.t;
            while(!(cv == start && ct == 0)){
                path.push_back(cv);
                auto b = back[key(cv,ct)];
                cv = b.prev_v;
                ct = b.prev_t;
            }
            path.push_back(start);
            std::reverse(path.begin(), path.end());
            return path;
        }

        auto try_relax = [&](int from, int to, int w) {
            int nxt_g = cur.g + w;
            int nxt_t = cur.t + w;
            int nxt_f = nxt_g + h[to];

            // 1) вершина занята
            if (busy_vertex(nxt_t,  to)) return;
            // 2) head‑on swap
            if (busy_vertex(nxt_t,  from) &&
                busy_vertex(cur.t,  to))   return;

            uint64_t k = key(to, nxt_t);
            if (!bestG.count(k) || nxt_g < bestG[k]) {
                bestG[k] = nxt_g;
                back[k]  = {from, cur.t};
                open.push({to, nxt_g, nxt_t, nxt_f});
            }
        };

        /* 1) реальные рёбра */
        for(const Edge& e : g[cur.v]){
            try_relax(cur.v, e.target, e.weight);
        }
        /* 2) действие «подождать 1 tick» */
        try_relax(cur.v, cur.v, 1);
    }
    return {};           // пустой – пути нет
}

// ───────────── ЧТЕНИЕ ФАЙЛОВ ─────────────

std::vector<std::vector<Edge>> read_graph(const std::string& fname)
{
    std::ifstream f(fname);
    if (!f) { std::cerr << "Graph.txt?\n"; std::exit(1); }

    /* ── читаем ПЕРВУЮ строку целиком ───────────────────────── */
    std::string line;
    if (!std::getline(f, line)) { std::cerr << "пустой файл\n"; std::exit(1); }

    std::istringstream first(line);
    std::vector<int> row;
    for (int w; first >> w;) row.push_back(w);

    int n = row.size();
    std::vector<std::vector<Edge>> g(n);

    /* заполняем первую строку */
    for (int j = 0; j < n; ++j)
        if (row[j] > 0) g[0].push_back({j, row[j]});

    /* читаем ещё (n‑1) строк, игнорируя всё, что идёт ниже */
    for (int i = 1; i < n; ++i)
    {
        std::getline(f, line);
        if (line.empty()) { --i; continue; }          // пропускаем пустые

        std::istringstream is(line);
        for (int j = 0, w; j < n && is >> w; ++j)
            if (w > 0) g[i].push_back({j, w});
    }
    return g;
}


// ───── Scenario + чтение ─────
struct Scenario {
    int total, nMines, nTerms, nAgents;
    std::vector<int> ids;    // список id длиной P
    std::vector<int> mines;  // шахты   длиной P
    std::vector<int> terms;  // терминалы длиной P
};

Scenario read_scenario(const std::string& fname)
{
    std::ifstream f(fname);
    if(!f){ std::cerr<<"terminal.txt?\n"; std::exit(1); }

    Scenario sc;
    if(!(f >> sc.total >> sc.nMines >> sc.nTerms >> sc.nAgents)){
        std::cerr<<"bad header\n"; std::exit(1);
    }

    auto read_line = [](std::ifstream& fs){
        std::vector<int> v;
        std::string line;          std::getline(fs >> std::ws, line);
        std::istringstream is(line);
        for(int x; is >> x;) v.push_back(x);
        return v;
    };

    sc.ids   = read_line(f);   // строка‑2
    sc.mines = read_line(f);   // строка‑3
    sc.terms = read_line(f);   // строка‑4

    size_t P = sc.ids.size();
    if(P == 0 || sc.mines.size()!=P || sc.terms.size()!=P){
        std::cerr<<"rows 2‑4 have different lengths\n"; std::exit(1);
    }
    return sc;
}

// ───── make_agents ─────
std::vector<Agent> make_agents(const Scenario& sc)
{
    std::vector<Agent> res(sc.nAgents);

    for(size_t i = 0; i < sc.ids.size(); ++i){
        int aid = sc.ids[i];
        if(aid >= sc.nAgents){
            std::cerr<<"id "<<aid<<" > nAgents\n"; std::exit(1);
        }
        res[aid].id = aid;
        res[aid].route.push_back(sc.mines[i]);
        res[aid].route.push_back(sc.terms[i]);
    }

    for(auto& ag : res){
        if(ag.route.empty()){
          ag.finished = true;
          continue;
        }
        ag.route.push_back(ag.route.front());   // замыкаем цикл
    }
    return res;
}

// ──────────────────  SIMULATE  ───────────────────
void simulate(std::vector<std::vector<Edge>>& g, std::vector<Agent>& agents)
{
    //std::ofstream log("log.txt");     ⬅️log⬅️

    std::vector<int> lastStanding(agents.size(), 0);

    Reserve globalRes;
    std::vector<int> agentClock(agents.size(), 1);

    bool allDone = false;
    while(!allDone)
    {
        allDone = true;
        for (auto& ag : agents) if (!ag.finished) { allDone = false; break; }
        if (allDone) break;

        // 1. планируем только “живых” агентов
        std::vector<std::vector<int>> paths(agents.size());
        for (size_t k=0; k<agents.size(); ++k)
        {
            auto& ag = agents[k];
            if (ag.finished) continue;

            auto heur = dijkstra(g, ag.goal());
            paths[k] = a_star_with_reserve(g,
                                           ag.start(), ag.goal(),
                                           heur,
                                           globalRes,
                                           agentClock[k]);

            if (paths[k].empty()) {
                //log << "Agent "<<ag.id<<" : no path\n";       ⬅️log⬅️
                ag.finished = true;          // чтобы не зациклиться
                continue;
            }

            /* резервируем путь */
            int t = 0;
            for (size_t s=0; s<paths[k].size(); ++s)
            {
                int v = paths[k][s];
                int absT = agentClock[k]+t;
                globalRes[absT].insert(v);
                if (s+1<paths[k].size())         // идём к следующей
                {
                    int u = paths[k][s+1];

                    if (u == v) continue;      // ⬅️ ожидание, ничего не прибавляем
                    int w = 1;
                    for (auto& e:g[v]) if (e.target==u){ w=e.weight; break; }
                    t += w;
                }
            }
            /* навечно резервируем финал */
            int arrival = agentClock[k]+t;
            globalRes[arrival    ].insert(ag.goal());   // момент прихода
            globalRes[arrival + 1].insert(ag.goal());   // вечная пауза
        }

        // 2. проигрываем и отмечаем завершившихся
        for (size_t k=0; k<agents.size(); ++k)
        {
            auto& ag = agents[k];
            if (ag.finished || paths[k].empty()) continue;

            int local=0;
            for (size_t s=0; s<paths[k].size(); ++s)
            {
                int v = paths[k][s];
                int gT = agentClock[k]+local;
                bool isMine     = (v == ag.route[ag.idx]);           // текущая шахта
                bool isTerminal = (v == ag.route[ag.idx+1]);         // текущий терминал

                /*      ⬅️log⬅️
                log << "Time " << gT << " : Agent " << ag.id
                    << " reached " << v;
                if (isTerminal) log << " (terminal)";
                else if (isMine) log << " (mine)";
                log << '\n';
                */

                if (s+1<paths[k].size()){
                    int u=paths[k][s+1], w=1;
                    for(auto& e:g[v]) if(e.target==u){w=e.weight;break;}
                    local+=w;
                }
            }
            /* финальная секунда стояния */
            
            /*      ⬅️log⬅️
            log<<"Time "<<agentClock[k]+local<<" : Agent "<<ag.id
               <<" standing at "<<ag.goal()<<'\n';
            */
            lastStanding[ag.id] = std::max(lastStanding[ag.id], agentClock[k] + local);
            agentClock[k]+=local+1;   // +1 за паузу
            ag.advance();
            if (ag.idx==ag.route.size()-2) ag.finished=true;
        }
    }

    {
      std::ofstream pso("pso.txt");
      for (size_t i = 0; i < lastStanding.size(); ++i) {
          // +1, потому что по условию нужно «последний tick + 1»
          pso << lastStanding[i] + (lastStanding[i] ? 1 : 0);
          if (i + 1 < lastStanding.size()) pso << ' ';
      }
      pso << '\n';
    }

    //log.close();
}

void run_simulation(const std::string& graph_file, const std::string& scenario_file) {
    auto graph = read_graph(graph_file);
    auto sc = read_scenario(scenario_file);
    auto agents = make_agents(sc);
    simulate(graph, agents);
}

// ──────────────────────────  main()  ───────────────────────────
//int main()
//{
//    run_simulation("Graph1.txt", "terminal1.txt");
//
//    std::cout << "Готово, смотрите pso.txt\n";
//    return 0;
//}
