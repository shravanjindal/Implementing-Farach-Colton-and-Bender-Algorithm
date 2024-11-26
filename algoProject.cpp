#include<bits/stdc++.h>
using namespace std;

// Tree Representation
int n; // Number of nodes
vector<vector<int>> tree; // Adjacency list

// Variables for LCA and shortest distance
int block_size, block_count;
vector<int> first_visit;
vector<int> euler_tour;
vector<int> depth; // Stores depth of each node
vector<int> log_values;
vector<vector<int>> sparse_table;
vector<vector<vector<int>>> precomputed_blocks;
vector<int> block_masks;

// Depth-First Search (DFS) to generate Euler Tour and depths
void perform_dfs(int current, int parent, int current_depth) {
    first_visit[current] = euler_tour.size(); // Record first visit index
    euler_tour.push_back(current); // Add node to Euler Tour
    depth[current] = current_depth;

    for (int neighbor : tree[current]) {
        if (neighbor == parent)
            continue;
        perform_dfs(neighbor, current, current_depth + 1); // Visit child
        euler_tour.push_back(current); // Add current node again after returning
    }
}

// Helper function to compare nodes by depth in Euler Tour
int get_min_by_depth(int index1, int index2) {
    return depth[euler_tour[index1]] < depth[euler_tour[index2]] ? index1 : index2;
}

// Preprocessing to compute LCA structure
void initialize_lca(int root) {
    // Step 1: Initialize data structures
    first_visit.assign(n, -1);
    depth.assign(n, 0);
    euler_tour.reserve(2 * n); // Max size of Euler Tour

    // Step 2: Generate Euler Tour and depths via DFS
    perform_dfs(root, -1, 0);

    // Step 3: Precompute log values
    int euler_size = euler_tour.size();
    log_values.resize(euler_size + 1, -1);
    for (int i = 1; i <= euler_size; i++)
        log_values[i] = log_values[i / 2] + 1;

    // Step 4: Block decomposition
    block_size = max(1, log_values[euler_size] / 2); // Block size is log/2
    block_count = (euler_size + block_size - 1) / block_size; // Number of blocks

    // Step 5: Precompute minimum in each block and sparse table
    sparse_table.assign(block_count, vector<int>(log_values[block_count] + 1));
    for (int i = 0, block_index = 0, block_offset = 0; i < euler_size; i++, block_offset++) {
        if (block_offset == block_size)
            block_offset = 0, block_index++;
        if (block_offset == 0 || get_min_by_depth(i, sparse_table[block_index][0]) == i)
            sparse_table[block_index][0] = i;
    }
    for (int level = 1; level <= log_values[block_count]; level++) {
        for (int i = 0; i < block_count; i++) {
            int next_block = i + (1 << (level - 1));
            if (next_block >= block_count)
                sparse_table[i][level] = sparse_table[i][level - 1];
            else
                sparse_table[i][level] = get_min_by_depth(sparse_table[i][level - 1], sparse_table[next_block][level - 1]);
        }
    }

    // Step 6: Mask generation for blocks
    block_masks.assign(block_count, 0);
    for (int i = 0, block_index = 0, block_offset = 0; i < euler_size; i++, block_offset++) {
        if (block_offset == block_size)
            block_offset = 0, block_index++;
        if (block_offset > 0 && (i >= euler_size || get_min_by_depth(i - 1, i) == i - 1))
            block_masks[block_index] |= 1 << (block_offset - 1);
    }

    // Step 7: Precompute RMQs for all unique blocks
    int possibilities = 1 << (block_size - 1);
    precomputed_blocks.resize(possibilities);
    for (int block_index = 0; block_index < block_count; block_index++) {
        int mask = block_masks[block_index];
        if (!precomputed_blocks[mask].empty())
            continue;
        precomputed_blocks[mask].assign(block_size, vector<int>(block_size));
        for (int l = 0; l < block_size; l++) {
            precomputed_blocks[mask][l][l] = l;
            for (int r = l + 1; r < block_size; r++) {
                precomputed_blocks[mask][l][r] = precomputed_blocks[mask][l][r - 1];
                if (block_index * block_size + r < euler_size)
                    precomputed_blocks[mask][l][r] = get_min_by_depth(
                        block_index * block_size + precomputed_blocks[mask][l][r],
                        block_index * block_size + r
                    ) - block_index * block_size;
            }
        }
    }
}

// Retrieve LCA from within a block
int lca_within_block(int block, int left, int right) {
    return precomputed_blocks[block_masks[block]][left][right] + block * block_size;
}

// Find the Lowest Common Ancestor (LCA) of two nodes
int find_lca(int node1, int node2) {
    int left = first_visit[node1];
    int right = first_visit[node2];
    if (left > right) swap(left, right);

    int block_left = left / block_size;
    int block_right = right / block_size;

    if (block_left == block_right)
        return euler_tour[lca_within_block(block_left, left % block_size, right % block_size)];

    int left_min = lca_within_block(block_left, left % block_size, block_size - 1);
    int right_min = lca_within_block(block_right, 0, right % block_size);
    int result = get_min_by_depth(left_min, right_min);

    if (block_left + 1 < block_right) {
        int range_log = log_values[block_right - block_left - 1];
        int middle_min1 = sparse_table[block_left + 1][range_log];
        int middle_min2 = sparse_table[block_right - (1 << range_log)][range_log];
        result = get_min_by_depth(result, get_min_by_depth(middle_min1, middle_min2));
    }
    return euler_tour[result];
}

// Utility to compute the shortest distance between two nodes
int compute_shortest_distance(int node1, int node2) {
    int lca_node = find_lca(node1, node2);
    return depth[node1] + depth[node2] - 2 * depth[lca_node];
}

// Main function to test LCA and shortest distance
int main() {
    // Tree structure
    n = 13;
    tree = {
        {1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12},
        {}, {}, {}, {}, {}, {}, {}, {}, {}
    };

    // Precompute LCA data
    initialize_lca(0);

    // Test LCA queries
    cout << "LCA of 4 and 5: " << find_lca(4, 5) << endl;
    cout << "LCA of 7 and 9: " << find_lca(7, 9) << endl;

    // Test shortest distance
    cout << "Shortest distance between 4 and 7: " << compute_shortest_distance(4, 7) << endl;

    return 0;
}

/*
                                            0
                                            |
                            1               2               3
                        /   |   \       /   |   \       /   |   \
                        4   5   6       7   8   9       10  11  12
*/