#include "tsp.hpp"

std::function<bool(double)> isinf_tsp = [](double elem) { return std::isnan(elem); };

std::vector<int> tsp(std::vector<std::vector<double>> cost_matrix){
    double LB=0;
    ivec result;
    ipair_list_t idx_list;
    ipair_dict_t pair;
    ipair_t step;
    for(int loop=0; loop<cost_matrix.size()-2; loop++){
        reduce_min_rows(cost_matrix, LB);
        reduce_min_cols(cost_matrix, LB);
        pair.clear();
        find_pairs(cost_matrix, pair, idx_list);
        step = find_best_pair(pair);
        idx_list.push_back(step);
        make_next_matrix(cost_matrix, step);
    }
    solve_matrix_2x2(cost_matrix, LB, idx_list);
    if(LB==0) {
        solve_zeros_matrix(cost_matrix, result, idx_list);
        return result;
    }
    result = sort_result(idx_list);
    return result;
}

double get_forbidden_cost(){return INF;}

void reduce_min_rows(matrix_t &cost_matrix, double &LB){
    for (int i = 0; i < cost_matrix.size(); i++) {
        double min = find_min_row(cost_matrix,i);
        if(isinf_tsp(min)){min=0;}
        LB += min;
        for (int idx = 0; idx < cost_matrix.size(); idx++) {
            cost_matrix[i][idx] -= min;
        }
    }
}

void reduce_min_cols(matrix_t &cost_matrix, double &LB){
    for (int j = 0; j < cost_matrix.size(); j++) {
        double min = find_min_col(cost_matrix,j);
        if(isinf_tsp(min)){min=0;}
        LB += min;
        for (int idx = 0; idx < cost_matrix.size(); idx++) {
            cost_matrix[idx][j] -= min;
        }
    }
}

double find_min_row(matrix_t &cost_matrix, int row_idx){
    double min = cost_matrix[row_idx][0];
    if (isinf_tsp(min)) {
        int idx = 0;
        while (isinf_tsp(min) and idx < cost_matrix.size()-1) {
            idx++;
            min = cost_matrix[row_idx][idx];
        }
    }
    for (int j = 0; j < cost_matrix.size(); j++) {
        if (cost_matrix[row_idx][j] < min) {
            min = cost_matrix[row_idx][j];
        }
    }
    return min;
}
double find_min_col(matrix_t &cost_matrix, int col_idx){
    double min = cost_matrix[0][col_idx];
    if (isinf_tsp(min)){
        int idx=0;
        while(isinf_tsp(min) and idx < cost_matrix.size()-1){
            idx++;
            min = cost_matrix[idx][col_idx];
        }
    }
    for (int i = 0; i < cost_matrix.size(); i++) {
        if (cost_matrix[i][col_idx] < min) {
            min = cost_matrix[i][col_idx];
        }
    }
    return min;
}


void find_pairs(matrix_t &cost_matrix, ipair_dict_t &pair, ipair_list_t &idx_list){
    double arg=0;
    for (int i = 0; i < cost_matrix.size(); i++) {
        for (int j = 0; j < cost_matrix.size(); j++) {
            if (cost_matrix[i][j] == 0) {
                cost_matrix[i][j] = INF;
                double min = find_min_row(cost_matrix, i);
                arg += min;
                min = find_min_col(cost_matrix, j);
                arg += min;
                cost_matrix[i][j] = 0;
                if(!isinf_tsp(arg)){pair[std::make_pair(i+1, j+1)] = arg;}
                arg=0;
            }
        }
    }
}

ipair_t find_best_pair(ipair_dict_t &pair) {
    double max = -1;
    ipair_t best;
    for (auto &next : pair) {
        if (pair.size() == 1) {
            best = next.first;
        } else {
            if (max < next.second) {
                max = next.second;
                best = next.first;
            }
        }
    }
    return best;
}

void make_next_matrix(matrix_t &cost_matrix, ipair_t step){
    for (int i = 0; i < cost_matrix.size(); i++) {
        cost_matrix[i][step.second-1] = INF;
        cost_matrix[step.first-1][i] = INF;
    }
    cost_matrix[step.second-1][step.first-1]=INF;
}

void solve_matrix_2x2(matrix_t &cost_matrix, double &LB, ipair_list_t &idx_list) {
    reduce_min_rows(cost_matrix, LB);
    reduce_min_cols(cost_matrix, LB);
    std::vector<ipair_t> last;
    for (int i = 0; i < cost_matrix.size(); i++) {
        for (int j = 0; j < cost_matrix.size(); j++) {
            if (!isinf_tsp(cost_matrix[i][j])) {
                last.push_back(std::make_pair(i + 1, j + 1));
            }
        }
    }
    auto pair_1 = last.begin();
    auto pair_2 = last.begin() + 1;
    auto pair_3 = last.begin() + 2;
    if (pair_1->first==pair_2->first){
        if (pair_1->second==pair_3->second){
            idx_list.push_back(std::make_pair(pair_2->first, pair_2->second));
            idx_list.push_back(std::make_pair(pair_3->first, pair_3->second));
        }else{
            idx_list.push_back(std::make_pair(pair_1->first, pair_1->second));
            idx_list.push_back(std::make_pair(pair_3->first, pair_3->second));
        }
    }
    if (pair_1->first==pair_3->first){
        if (pair_1->second==pair_2->second){
            idx_list.push_back(std::make_pair(pair_2->first, pair_2->second));
            idx_list.push_back(std::make_pair(pair_3->first, pair_3->second));
        }else{
            idx_list.push_back(std::make_pair(pair_1->first, pair_1->second));
            idx_list.push_back(std::make_pair(pair_2->first, pair_2->second));
        }
    }
    if (pair_2->first==pair_3->first){
        if (pair_2->second==pair_1->second){
            idx_list.push_back(std::make_pair(pair_1->first, pair_1->second));
            idx_list.push_back(std::make_pair(pair_3->first, pair_3->second));
        }else{
            idx_list.push_back(std::make_pair(pair_1->first, pair_1->second));
            idx_list.push_back(std::make_pair(pair_2->first, pair_2->second));
        }
    }
}

ivec sort_result(ipair_list_t&idx_list){
    ivec result;
    auto first = idx_list.begin();
    result.push_back(first->first);
    result.push_back(first->second);
    while(result.size() != idx_list.size()+1){
        for(auto &elem : idx_list){
            if (elem.first == result.back()){
                result.push_back(elem.second);
            }
        }
    }
    return result;
}

void solve_zeros_matrix(matrix_t &cost_matrix, ivec &result, ipair_list_t &idx_list){
    result = sort_result(idx_list);
    result.erase(result.end()-1);
    for (int z = 0; z!=cost_matrix.size(); z++) {
        if(std::find(result.begin(), result.end(), z+1) == result.end()){
            result.push_back(z+1);
        }
    }
    result.push_back(result[0]);
}