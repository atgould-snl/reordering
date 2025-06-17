#ifndef BUCKET_SOLVER_H
#define BUCKET_SOLVER_H

#include "common.h"
#include "blackBoxLop.h"

#include <vector>
#include <set>
#include <map>
#include <utility>
#include <limits>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <queue>

class BucketOrderingSolver; // Forward declaration

class BucketingOption {
public:
    BucketingOption(const BucketOrderingSolver & prob, bool mergeAll=false);

    // Friends
    friend class BucketOrderingSolver;

    // Getters
    [[nodiscard]] auto get_cost() const {return cost;};
    [[nodiscard]] auto get_mergeCost() const {return mergeCost;};
    [[nodiscard]] auto get_nStableBuckets() const {return nStableBuckets;};
    [[nodiscard]] auto get_nMerge() const {return nMerge;};
    [[nodiscard]] auto get_map() const {return physics_to_block_map;};
    [[nodiscard]] auto get_optimisticMergeCostToHitTarget() const {return optimisticMergeCostToHitTarget;};
    bool is_leaf() const {return nStableBuckets==int(buckets.size());}

    std::vector<BucketingOption> makeChildren() const;
    void merge(int a_star, int b);
    void newBucket(int a_star);
    
    // Overload the < operator for BucketingOptions
    bool operator<(const BucketingOption& other) const;

private:
    std::vector<std::set<int>> buckets;
    int nStableBuckets;
    int nMerge=0;
    std::map<int, int> physics_to_block_map;
    std::vector<int> order_of_rows;
    std::vector<int> order;
    double cost;
    double mergeCost;
    double optimisticMergeCostToHitTarget;
    const BucketOrderingSolver* problem;

    void runLOP();
    double bucketSum(const int a, const int b, const bool mergeCosts=false) const;
    std::pair<int, double> get_a_star_L_star() const;
    void updateOrder();
    void updateMap();
    void updateCost();
    void updateMergeCost();
    void updateOptimisticMergeCostToHitTarget();
};

class BucketOrderingSolver {
public:
    BucketOrderingSolver(const BlockNormsViewType & blockNorms, const BlockNormsViewType & blockMergeCosts, const double costTarget, const double tMaxWalltime, bool exhaustive_mode=false);

    // Friends
    friend class BucketingOption;

    [[nodiscard]] auto get_best() const {return best;};
    [[nodiscard]] auto get_leafNodes() const {return leaf_nodes;};
    [[nodiscard]] auto get_internalNodes() const {return internal_nodes;};
    
    void solve();

private:
    const BlockNormsViewType blockNorms;
    const BlockNormsViewType blockMergeCosts;
    const double tMaxWalltime;
    const double costTarget;
    const int N;

    int leaf_nodes=0;
    int internal_nodes=1;

    BucketingOption best;
    BucketingOption base;
    bool exhaustive_mode;
};

#endif // BUCKET_SOLVER_H