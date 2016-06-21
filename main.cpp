#include "vp.h"
#include <cmath>
#include <iostream>
#include <numeric>
#include <cassert>
#include <omp.h>
#include <fstream>
#include <map>


size_t NO_PLACE = -1;


struct TPoint {
    double x;
    double y;
    size_t place_id;
    int minute;
    int wday;
    int week;
    int time;
    double accuracy;


    TPoint(double x, double y, size_t place, int time, double accuracy) :
        x(x), y(y), place_id(place), minute(time % 1440), wday((time / 1440) % 7), week((time / 1440 /  7)), time(time), accuracy(accuracy)
    {}

    TPoint(double x, double y, int time, double accuracy) :
        x(x), y(y), place_id(NO_PLACE), minute(time % 1440), wday((time / 1440) % 7), week((time / 1440 /  7)), time(time), accuracy(accuracy)
    {}

    TPoint(const TPoint& other) : 
        x(other.x), y(other.y), place_id(other.place_id), minute(other.minute), wday(other.wday), week(other.week), time(other.time), accuracy(other.accuracy)
    {}

    double cycle_div(int a, int b, int period) const {
        return std::min(std::min(a, b) + period - std::max(a, b), std::max(a, b) - std::min(a, b));
    }

    TPoint operator - (const TPoint& other) const {
        auto ans = TPoint(*this);
        ans.x = std::abs(x - other.x);
        ans.y = std::abs(y - other.y);
        ans.place_id = NO_PLACE;
        ans.minute =  cycle_div(minute, other.minute, 1440);
        ans.accuracy = std::min(accuracy, other.accuracy);
        ans.wday = cycle_div(wday, other.wday, 7);
        ans.time = std::abs(time - other.time);
        ans.week = std::abs(week - other.week);
        return ans;
    }

    double abs() const {
        return sqrt(x * x + 4 * y * y + double(minute * minute) / (60 * 60 * 125 * 125) + double(wday * wday) / (500 * 500) + double(time * time) / (size_t(60) * 24 * 7 * 500 * 60 * 24 * 7 * 500));
    }
};


inline double dist(const TPoint &p1, const TPoint &p2) {
    return (p1 - p2).abs();
}


inline double mapn(size_t targ, const std::vector<size_t>& targs, size_t n) {
    double map = 0;
    bool found = false;
    n = std::min(n, targs.size());
    for (size_t i = 0; i < n; ++i) {
        found |= targ == targs[i]; 
        map += double(found) / n;
    }

    return map;
}

inline double proxmity_score(size_t targ, const std::vector<TPoint>& targs, size_t n) {
    double map = 0;
    double found = 0.0;
    n = std::min(n, targs.size());
    for (size_t i = 0; i < n; ++i) {
        found += (targ == targs[i].place_id); 
        map += double(found) / n;
    }

    return map;
}

std::map<size_t, std::vector<double>> get_time_coefs(const std::vector<TPoint>& pts) {
    // std::vector<size_t> total_counter(24, 0);
    std::map<size_t, std::vector<double>> result;
    std::map<size_t, size_t> total;

    size_t shift = std::min(pts.size(), size_t(3000000));
    for (size_t i = (pts.size() - shift); i < pts.size(); ++ i) {
        int hour = (pts[i].minute / 60) % 24;
        // total_counter[hour] += 1;
        size_t place_id = pts[i].place_id;
        
        if (result.find(place_id) == result.end()) {
            result[place_id] = std::vector<double>(24, 1);
            total[place_id] = 24;
        }

        total[place_id] += 1;
        result[place_id][hour] += 1;
    }
    for (const auto& a : result) {
        for (size_t i = 0; i < 24; ++i) {
            result[a.first][i] /= total[a.first];
        }
    }
    return result;
}

inline std::vector<size_t> rank(
    const std::vector<TPoint>& targs,
    const std::map<size_t, std::vector<double>>& time_coefs,
    int minute) {

    // std::vector<size_t> places(targs.size(), size_t(0)); 
    // for (const auto& a : targs) places.emplace_back(a.place_id);

    std::map<size_t, double> ll;
    // for (const auto& a: targs) ll[a.place_id] = mapn(a.place_id, places, 20);
    for (const auto& a: targs) ll[a.place_id] = proxmity_score(a.place_id, targs, 100);

    std::vector<std::pair<size_t, double>> pairs;
    for (const auto& a : ll) {
        auto it = time_coefs.find(a.first);
        if (it != time_coefs.end()) {
            pairs.push_back(std::make_pair(a.first, a.second * it->second[(minute / 60) % 24])); // 
        } else {
            pairs.push_back(std::make_pair(a.first, a.second / 24)); 
        }
    }

    sort(pairs.rbegin(), pairs.rend(), [=](std::pair<size_t, double>& a, std::pair<size_t, double>& b) {
       return a.second < b.second;
    });

    std::vector<size_t> ans;
    for (const auto& a : pairs)
        ans.push_back(a.first);
    
    return ans;
}



std::vector<TPoint> read_file(char *fname, bool no_ans){
    double start = omp_get_wtime();

    std::cout << "reading : " << fname << std::endl;

    std::vector<TPoint> data;

    std::string buf;
    double x, y, time, acc;
    size_t place_id = 0, cnt = 0;

    std::ifstream fin(fname);
    getline (fin, buf);

    while (getline (fin, buf,',')) {
        cnt += 1;
        if (!(cnt % 2000000))
            std::cout << '\t' << cnt / 1000000 << "M records done\n";

        getline (fin, buf, ',' );
        x = std::stod(buf);
        getline (fin, buf, ',' );
        y = std::stod(buf);
        getline (fin, buf, ',' );
        acc = std::stod(buf);
        if (no_ans) {
            getline (fin, buf);
            time = std::stol(buf);
            place_id = NO_PLACE;
        } else {
            getline (fin, buf, ',' );
            time = std::stol(buf);
            getline (fin, buf);
            place_id = std::stoll(buf);
        }

        data.emplace_back(x, y, place_id, time, acc);
    }

    double end = omp_get_wtime();
    std::cout << "creating particles took: " << end - start << "sec" << std::endl;

    return data;
}

int main(int argc, char ** argv) {
    std::vector<TPoint> test  = read_file(argv[2], argc > 4);
    std::cout << test.size() << " lines in test" << std::endl;
    std::vector<TPoint> train = read_file(argv[1], false);
    std::cout << train.size() << " lines in train" << std::endl;
    std::map<size_t, std::vector<double>> time_coefs = get_time_coefs(train);

    VpTree<TPoint, dist> tree;

    double start = omp_get_wtime();
    tree.create(train);
    train.clear();
    double end = omp_get_wtime();
    std::cout << "creating tree took: " << end - start << "sec" << std::endl;

    int k = 100;

    start = omp_get_wtime();

    double precision = 0.0;

    std::ofstream fout(argv[3]);
    fout << "row_id,place_id" << std::endl;
    for (size_t i = 0; i < test.size(); i++) {
        std::vector<double> distances;
        std::vector<TPoint> neighbors;
        tree.search(test[i], k, &neighbors, &distances);

        fout << i << ",";
        std::vector<size_t> ans = rank(neighbors, time_coefs, test[i].minute);
        for (size_t j = 0; j < std::min(ans.size(), size_t(3)); ++j) 
            fout << ans[j] << " ";

        fout << std::endl;

        precision += mapn(test[i].place_id, ans, 3);

        if (!(i % 40000)) {
            std::cout << i << " samples, precision: " <<  precision / i << std::endl;
            // for (const auto& a: ans) std::cout << a << " ";
            // std::cout << std::endl;
        }
    }

    end = omp_get_wtime();
    std::cout << "scoring took: " << end - start << " seconds" << std::endl;

    return 0;
}