#include "branchAndBound.h"
#include <limits>
#include <vector>

void STK_ThrowRequire(bool condition) {
    if (condition) {
        throw std::runtime_error("Condition not met");
    }
}

// Constructor
BucketingOption::BucketingOption(const BucketOrderingSolver* problem)
    : problem(problem), nStableBuckets(0), cost(0), mergeCost(0) {
    // Initialize buckets with 1 physics each
    for (int physics_number = 0; physics_number < problem->N; ++physics_number) {
        buckets.push_back(std::set<int>{physics_number});
    }
    runLOP();
}

std::vector<BucketingOption> BucketingOption::makeChildren() const {
  // Make children by greedy merge selection
  // First find L*, the larges lower diagonal term associated with a non-stable bucket
  auto [a_star, L_star] = get_a_star_L_star();
  // Create a child where a_star is merged to each stable block, tracking the potential of each.
  std::vector<BucketingOption> children = std::vector<BucketingOption>(nStableBuckets+1, *this); // Make copies of self
  for (unsigned long b; b < nStableBuckets; b++){
    children[b].merge(a_star,b); // Also runs LOP
  }
  // Create child where a_star gets its own bucket
  children[nStableBuckets].newBucket(a_star);
  return std::move(children); // Move to not make a copy
}

void BucketingOption::merge(int a_star, int b) {
  STK_ThrowRequire(a_star > b);
  STK_ThrowRequire(a_star > nStableBuckets - 1);
  STK_ThrowRequire(a_star < int(buckets.size()));
  STK_ThrowRequire(1 == int(buckets[a_star].size()));
  buckets[b].insert(buckets[a_star].begin(), buckets[a_star].end());
  buckets.erase(buckets.begin() + a_star);
  runLOP();
}

void BucketingOption::newBucket(int a_star) {
  STK_ThrowRequire(a_star > nStableBuckets - 1);
  STK_ThrowRequire(a_star < int(buckets.size()));
  STK_ThrowRequire(1 == int(buckets[a_star].size()));
  std::set<int> tempBucket = buckets[nStableBuckets];
  buckets[nStableBuckets] = buckets[a_star];
  buckets[a_star] = tempBucket; // Ok if a_star = nStableBuckets
  nStableBuckets++;
  runLOP(); // TODO: Easy optimization, this could be taken away since the order is impacted in a trivial way
}

//bool BucketingOption::operator<(const BucketingOption& other) const {
//  if (mergeCost < other.get_mergeCost()) {
//      return get_nMerge() < other.get_nMerge(); // Merge count is first priority
//  }
//  return cost < other.get_cost(); // Second look to cost
//}

void BucketingOption::runLOP() {
    updateOrder(); // RUNS LOP
    updateCost(); // Based on new order
    updateMergeCost(); // Based on bucketing only
    updateOptimisticMergeCostToHitTarget();
    updateMap(); // Based on new order
}

double BucketingOption::bucketSum(const int a, const int b, const bool mergeCosts) const {
  double T = 0.0; // Initialize T
  STK_ThrowRequire(a!=b);
  for (const auto& i : buckets[a]) { 
      for (const auto& j : buckets[b]) {
          T += mergeCosts ? problem->blockMergeCosts(i,j)+problem->blockMergeCosts(i,j) : std::pow(problem->blockNorms(i, j), 2);
      }
  }
  return T;
}

std::pair<int, double> BucketingOption::get_a_star_L_star() const {
  double L_star = std::numeric_limits<double>::max();
  int a_star = 0;
  for (unsigned long a = nStableBuckets; a < buckets.size(); a++) { // Loop through non-stable only
      for (unsigned long b = 0; b < buckets.size(); b++) {
          if (order[a] <= order[b]) { continue; } // Filter out nonLD blocks
          double L = bucketSum(a, b);
          if (L < L_star) { L_star = L; a_star = a; }
      }
  }
  return std::pair<int, double>(a_star, L_star);
}

void BucketingOption::updateOrder() {
  BlockNormsViewType tournament("normalForm", buckets.size(), buckets.size());
  // Tournament based on the buckets
  for (unsigned long a=0; a < buckets.size(); a++){
    for (unsigned long b=0; b < buckets.size(); b++){
      if (a>b) { continue;} // Block is below diag, filter out
      if (a==b) { 
        tournament(a,a)=0;
        continue;
      }
      // Sum block above diag
      double T=0;
      T += bucketSum(a,b);

      // Take away sum of block below diag
      T -= bucketSum(b, a);
      tournament(a,b) = T;
      tournament(b,a) = -T;
    }
  }
  // Call the LOP solver
  LinearOrderingSolver solver{tournament, problem->tMaxWalltime}; // Is communicator pulled in from the namespace??
  solver.solve();
  order_of_rows = solver.order();
  std::vector<int> order(order_of_rows.size());
  std::iota( order_of_rows.begin(), order_of_rows.end(), 0);
  std::sort(order.begin(), order.end(), [this](int a, int b) { return order_of_rows[a] < order_of_rows[b];});
}

void BucketingOption::updateMap() {
  for (unsigned long bucket_rank=0; bucket_rank<buckets.size(); bucket_rank++){
    auto bucket = buckets[order[bucket_rank]];
    for (const auto& physics_block : bucket) {
      physics_to_block_map[physics_block] = int(bucket_rank);
    }
  }
}

void BucketingOption::updateCost(){ // Find the score based on the order
  cost=0;
  for (unsigned long a=0; a < buckets.size(); a++){
    for (unsigned long b=0; b < buckets.size(); b++){
      cost += order[a]>order[b] ? bucketSum(a,b) : 0;
    }
  }
}

void BucketingOption::updateMergeCost(){ // Find the score based on the order
  mergeCost=0;
  for (unsigned long a=0; a < buckets.size(); a++){
    for (const auto& i : buckets[a]) {
      for (const auto& j : buckets[a]) {
        if (i!=j){
          mergeCost += problem->blockMergeCosts(i,j);
        }
      }
    }
  }
}

void BucketingOption::updateOptimisticMergeCostToHitTarget(){ // Lower bound for bb tree
  optimisticMergeCostToHitTarget = mergeCost;
  if (cost < problem->costTarget) {return;}
  double optimisticCost = cost; // Start at current cost
  struct L_option_struct {
    double cost;
    double mergeCost;
    L_option_struct(double c, double m) : cost(c), mergeCost(m) {}
  };
  std::vector<L_option_struct> L_options;
  // Loop through the buckets and fill out L options
  for (unsigned long a = nStableBuckets; a < buckets.size(); a++) { // Loop through non-stable only
    for (unsigned long b = 0; b < buckets.size(); b++) {
        if (order[a] <= order[b]) { continue; } // Filter out nonLD blocks
        
        double L_cost = bucketSum(a, b);
        double L_mergeCost = bucketSum(a,b,L_mergeCost=true);
        L_options.push_back(L_option_struct(L_cost,L_mergeCost));
    }
  }
  // Sort by ratio
  std::sort(L_options.begin(), L_options.end(),[](const L_option_struct L1, const L_option_struct L2){return L1.cost/L1.mergeCost > L2.cost/L2.mergeCost;}); // Better ratio gets a true and goes first
  // Start adding from the bottom until problem->costTarget would be reached. Then add the last one fractionaly
  int i=0;
  while (i<L_options.size()){
    // Will adding the new cost take us under the target?
    if (optimisticCost-L_options[i].cost < problem->costTarget){
      // Subtract fraction and return
      double fraction = L_options[i].cost / (optimisticCost - problem->costTarget);
      STK_ThrowRequire(fraction<1.); STK_ThrowRequire(fraction>0.);
      optimisticMergeCostToHitTarget += fraction * L_options[i].mergeCost;
      return;
    }
    optimisticMergeCostToHitTarget += L_options[i].mergeCost;
    optimisticCost -= L_options[i].cost;
    i++;
  }
  // This means that no future merges will ever hit the target. No amount of merge cost will be sufficient
  optimisticMergeCostToHitTarget = std::numeric_limits<double>::max()/2;
}






BucketOrderingSolver::BucketOrderingSolver(const BlockNormsViewType & blockNorms, const BlockNormsViewType & blockMergeCosts, double costTarget, double tMaxWalltime)
    : blockNorms(blockNorms), blockMergeCosts(blockMergeCosts), costTarget(costTarget), tMaxWalltime(tMaxWalltime), N(int(blockNorms.extent(0))) {
    BucketingOption base = BucketingOption(this); // Pass this to BucketingOption
    addToOrder(base); // Recursive call
}

void BucketOrderingSolver::addToOrder(const BucketingOption & base) {
  std::vector<BucketingOption> children = base.makeChildren();
  std::vector<int> childrenPriority=std::vector<int>(children.size());
  std::iota(childrenPriority.begin(), childrenPriority.end(), 0);
  // Sort so highest potential goes first
  std::sort(childrenPriority.begin(), childrenPriority.end(),[&children](const int c1, const int c2){return potentials[b1] > potentials[b2];}); // More potential gets a true and goes first
  for (auto i:childrenPriority){
    // Check if child with its potential could be better by comparing to current best
    children
    // If so run LOP on it
    // Check if it still could be better than best
  }
}