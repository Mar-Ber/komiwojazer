#ifndef KOMIWOJAZER_TSP_HPP
#define KOMIWOJAZER_TSP_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <map>
#include <list>
#include <algorithm>

#define INF (NAN)
using matrix_t = std::vector<std::vector<double>>;
using ipair_t = std::pair<int, int>;
using ipair_dict_t = std::map<std::pair<int, int>,double>;
using ipair_list_t = std::list<std::pair<int, int>>;
using ivec = std::vector<int>;

std::vector<int> tsp(std::vector<std::vector<double>> cost_matrix);
double get_forbidden_cost();

void reduce_min_rows(matrix_t &cost_matrix, double &LB);
void reduce_min_cols(matrix_t &cost_matrix, double &LB);
double find_min_row(matrix_t &cost_matrix, int row_idx);
double find_min_col(matrix_t &cost_matrix, int col_idx);
void find_pairs(matrix_t &cost_matrix, ipair_dict_t &pair, ipair_list_t &idx_list);
ipair_t find_best_pair(ipair_dict_t &pair);
void make_next_matrix(matrix_t &cost_matrix, ipair_t step);
void solve_matrix_2x2(matrix_t &cost_matrix, double &LB, ipair_list_t &idx_list);
ivec sort_result(ipair_list_t&idx_list);
void solve_zeros_matrix(matrix_t &cost_matrix, ivec &result, ipair_list_t &idx_list);

#endif //KOMIWOJAZER_TSP_HPP
