#include "branchAndBound.h"
#include "common.h"
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <list>
#include <queue>
#include <vector>

bool almost_equal(double a, double b, double tolerance = 1e-8) {
    return std::fabs(a - b) <= tolerance;
}

void print_matrix_2(const Kokkos::View<double**>& T) {
  const int numRows = T.extent(0);
  const int numCols = T.extent(1);

  std::cout << "\n";
  std::cout << "T = {\n";
  for (int i = 0; i < numRows; i++){
  std::cout << "    {";
  for (int j = 0; j < numCols; j++){ 
      // Print each element, followed by a comma if it's not the last element
      std::cout << std::fixed << std::setprecision(12) << T(i,j);
      if (j < numCols - 1) {
          std::cout << ", ";
      }
  }
  std::cout << "},\n";
  }
  std::cout << "};\n";
}

void print_map_2(const std::map<int, int>& myMap) {
    std::cout << "Map contents:\n";
    for (const auto& pair : myMap) {
        std::cout << "Key: " << pair.first << ", Value: " << pair.second << '\n';
    }
}

void STK_ThrowRequire(bool condition) {
    if (!condition) {
        throw std::runtime_error("Condition not met");
    }
}

// Constructor
BucketingOption::BucketingOption(const BucketOrderingSolver* problem, bool mergeAll)
    : problem(problem), nStableBuckets(0), cost(0), mergeCost(0) {
    if (mergeAll){
      // Merge all into one bucket
      buckets.emplace_back(std::set<int>());
      for (int physics_number = 0; physics_number < problem->N; ++physics_number) {
        buckets[0].insert(physics_number);
        nMerge++;
      }
    }
    else{
      // Initialize buckets with 1 physics each, unless they all start merged
      for (int physics_number = 0; physics_number < problem->N; ++physics_number) {
          buckets.push_back(std::set<int>{physics_number});
      }
    }
    runLOP();
}

std::vector<BucketingOption> BucketingOption::makeChildren() const {
  if (nStableBuckets == int(buckets.size())) return std::vector<BucketingOption>(); // Empty child vector if all buckets are stable already
  // Make children by greedy merge selection
  // First find L*, the largest lower diagonal term associated with a non-stable bucket
  auto [a_star, L_star] = get_a_star_L_star();
  // Create a child where a_star is merged to each stable block, tracking the potential of each.
  std::vector<BucketingOption> children = std::vector<BucketingOption>(nStableBuckets+1, *this); // Make copies of self
  for (int b = 0; b < nStableBuckets; b++){
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
  nMerge++;
  runLOP();
}

void BucketingOption::newBucket(int a_star) {
  STK_ThrowRequire(a_star > nStableBuckets - 1);
  STK_ThrowRequire(a_star < int(buckets.size()));
  STK_ThrowRequire(1 == int(buckets[a_star].size()));

  runLOP();
  std::vector<int> order_old = order;
  double cost_old = cost;
  // Manually update order of buckets since two of them switched indicies
  std::swap(buckets[a_star], buckets[nStableBuckets]);
  std::swap(order[a_star], order[nStableBuckets]);

  nStableBuckets++;
  std::vector<int> order_copy = order;
  updateCost();
  double cost_copy = cost;
  runLOP(); // TODO: Easy optimization, this could be taken away since the order is impacted in a
            // trivial way
  //for (size_t i = 0; i < order_copy.size(); ++i) {
  //  STK_ThrowRequire(order_copy[i] == order[i]);
  //  order_old[i];
  //}
  // STK_ThrowRequire(cost_copy == cost);
  
}

bool BucketingOption::operator<(const BucketingOption& other) const { // Sees if one is better than the other
  return optimisticMergeCostToHitTarget > other.get_optimisticMergeCostToHitTarget(); // > instead of < for priority queue
}

void BucketingOption::runLOP() {
  updateOrder(); // RUNS LOP
  updateCost(); // Based on new order
  updateMergeCost(); // Based on bucketing only
  updateOptimisticMergeCostToHitTarget();
  updateMap(); // Based on new order
}

double BucketingOption::bucketSum(const int a, const int b, const bool mergeCosts) const {
  double T = 0.0; // Initialize T
  //STK_ThrowRequire(a!=b);
  for (const auto& i : buckets[a]) { 
      for (const auto& j : buckets[b]) {
          //STK_ThrowRequire(i!=j);
          // For merge costs, include the cooresponding bucket above the diag too
          T += mergeCosts ? problem->blockMergeCosts(i,j)+problem->blockMergeCosts(i,j) : std::pow(problem->blockNorms(i, j), 2);
      }
  }
  return T;
}

std::pair<int, double> BucketingOption::get_a_star_L_star() const { // Find the largest L
  double L_star = -1;
  int a_star = -1;
  for (int a = 0; a < int(buckets.size()); a++) { // Loop through non-stable only
    for (int b = 0; b < int(buckets.size()); b++) {
      if ((a <= nStableBuckets -1) && (b <= nStableBuckets -1)) continue; // Stable zone, do not search here
      if (order[a] <= order[b]) continue; // Filter out nonLD blocks
      double L = bucketSum(a, b);
      if (L > L_star) { 
        L_star = L; 
        a_star = a>b ? a : b; // Call a_star the greater of a and b
      } 
    }
  }
  STK_ThrowRequire(a_star!=-1);
  return std::pair<int, double>(a_star, L_star);
}

void BucketingOption::updateOrder() {
  BlockNormsViewType tournament("normalForm", buckets.size(), buckets.size());
  // Tournament based on the buckets
  for (int a=0; a < int(buckets.size()); a++){
    for (int b=0; b < int(buckets.size()); b++){
      if (a>b) continue; // Block is below diag, filter out
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
  global_timer_LOP.start();
  LinearOrderingSolver solver{tournament, problem->tMaxWalltime}; // Is communicator pulled in from the namespace??
  solver.solve();
  global_LOP_call_counter++;
  global_timer_LOP.pause();
  order_of_rows = solver.order();
  order = std::vector<int>(order_of_rows.size());
  std::iota( order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [this](int a, int b) { return order_of_rows[a] < order_of_rows[b];});

  /*
  /////// DEBUG LOP ONLY ///////
  updateCost();
  updateMap();
  if (almost_equal(cost, 0.00956199)){
    std::cout << "\nCost: " << cost << std::endl;
    print_map_2(physics_to_block_map);
    std::cout << "TOURNAMENT:\n";
    print_matrix_2(tournament);

    BlockNormsViewType reformed_blockNorms("reformed_blockNorms", buckets.size(), buckets.size());
    // Tournament based on the buckets
    for (int a=0; a < int(buckets.size()); a++){
      for (int b=0; b < int(buckets.size()); b++){
        reformed_blockNorms(a,b) = std::sqrt(bucketSum(a,b));
      }
    }
    std::cout << "BLOCKNORMS:\n";
    print_matrix_2(reformed_blockNorms);

  }
  */
}

void BucketingOption::updateMap() {
  for (int bucket_idx=0; bucket_idx<int(buckets.size()); bucket_idx++){ // Loop through the buckets
    for (const auto& i : buckets[bucket_idx]) { // Loop through that bucket
      physics_to_block_map[i] = int(order[bucket_idx]); // Order says where it should go
    }
  }
}

void BucketingOption::updateCost(){ // Find the score based on the order
  cost=0;
  for (int a=0; a < int(buckets.size()); a++){
    for (int b=0; b < int(buckets.size()); b++){
      cost += order[a]>order[b] ? bucketSum(a,b) : 0;
    }
  }
}

void BucketingOption::updateMergeCost(){ // Find the score based on the order
  mergeCost=0;
  for (int a=0; a < int(buckets.size()); a++){
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
  if (cost <= problem->costTarget) return;
  struct L_option_struct {
    double cost;
    double mergeCost;
    L_option_struct(double c, double m) : cost(c), mergeCost(m) {}
  };
  std::vector<L_option_struct> L_options;
  double stableCost=0;
  // Loop through the buckets and fill out L options, a list of possible new merges
  for (int a = 0; a < int(buckets.size()); a++) { // Loop through non-stable only
    for (int b = 0; b < int(buckets.size()); b++) {
      
      if (order[a] <= order[b]) continue; // Filter out nonLD blocks
      // Stable zone, take sum, these costs are fixed
      if (a <= nStableBuckets -1 && b <= nStableBuckets -1) // Stable zone, do not search here
      {
        stableCost += bucketSum(a, b);
      }
      // Unstable zone, load into L options
      else
      {
        double L_cost = bucketSum(a, b);
        double L_mergeCost = bucketSum(a,b,L_mergeCost=true);
        L_options.push_back(L_option_struct(L_cost,L_mergeCost));
      }
    }
  }
  // Sort by ratio
  std::sort(L_options.begin(), L_options.end(),[](const L_option_struct L1, const L_option_struct L2){
    if (L1.mergeCost == 0) return true; // Filter for zeros. Zero merge costs should go first
    if (L2.mergeCost == 0) return false;
    return L1.cost/L1.mergeCost > L2.cost/L2.mergeCost;
  }); // Better ratio gets a true and goes first
  // Start adding from the bottom until problem->costTarget would be reached. Then add the last one fractionaly
  std::vector<double> optimisticSumCost = std::vector<double>(L_options.size()+1);
  optimisticSumCost[L_options.size()] = stableCost; // Initialize with stable cost
  for (int i=L_options.size()-1; i>=0; i--) // Starting at second from back, start adding each cost as backward cumsum
  { 
    optimisticSumCost[i]=optimisticSumCost[i+1]+L_options[i].cost;
  }
  int i=0;
  while (i<L_options.size() && L_options[i].cost > 0){
    // Will adding the new cost take us under the target?
    // Compute costGap. This is a recomputation of the cost without L_options
    if (optimisticSumCost[i+1] <= problem->costTarget){
      // Subtract fraction and return
      double fraction =  (optimisticSumCost[i] - problem->costTarget) / L_options[i].cost;
      STK_ThrowRequire(fraction <= 1.); STK_ThrowRequire(fraction >= 0.);
      optimisticMergeCostToHitTarget += fraction * L_options[i].mergeCost;
      return;
    }
    optimisticMergeCostToHitTarget += L_options[i].mergeCost;
    i++;
  }
  // This means that no future merges will ever hit the target. No amount of merge cost will be sufficient
  optimisticMergeCostToHitTarget = std::numeric_limits<double>::max()/2;
}

/////// BucketOrderingSolver ////////

BucketOrderingSolver::BucketOrderingSolver(const BlockNormsViewType & blockNorms, const BlockNormsViewType & blockMergeCosts, const double costTarget, const double tMaxWalltime, bool exhaustive_mode)
    : blockNorms(blockNorms), blockMergeCosts(blockMergeCosts), costTarget(costTarget), tMaxWalltime(tMaxWalltime), N(int(blockNorms.extent(0))), best(BucketingOption(this,true)), base(BucketingOption(this, false)), exhaustive_mode(exhaustive_mode) {
    BucketingOption base = BucketingOption(this); // Pass this to BucketingOption
}

void BucketOrderingSolver::solve(){
  std::vector<std::priority_queue<BucketingOption>> optionsQueues = std::vector<std::priority_queue<BucketingOption>>(N+1); // N+1 to account for level 0
  optionsQueues[0].push(base);

  bool all_empty=false;
  while (!all_empty){ // Until all queues are empty or TODO time is reached
    global_best_cost_record[global_best_cost_record.size()-1].push_back(best.mergeCost);
    all_empty=true;
    for (int d = 0; d<N+1 ; d++){
      if (optionsQueues[d].empty()) continue;
      all_empty=false;

      BucketingOption topOption = optionsQueues[d].top(); // Get the option with the greatest criteria
      optionsQueues[d].pop(); // Remove it from the queue

      // If the target is reached, see if this one is better than the current best?
      if (topOption.cost <= costTarget)
      {
        if (topOption.mergeCost < best.mergeCost || best.cost > costTarget) // Best should always be below target, since it defaults to all merged, 
        // but added the || condition in case default changes
        {
          best = topOption; // Make it the new best
        }
      }

      // Is this topOption worth running? Is the best below target and the best is better than the top option?
      // If not, cut the branch
      // We can empty the whole queue
      if ((best.cost <= costTarget) && (topOption.optimisticMergeCostToHitTarget >= best.mergeCost) && (!exhaustive_mode))
      {
        optionsQueues[d]=std::priority_queue<BucketingOption>();
        continue;
      }

      // Otherwise make children from the top option and add them to queue
      auto children = topOption.makeChildren();
      for (const auto & child : children){
        if (child.is_leaf()) leaf_nodes++;
        else internal_nodes++;
        STK_ThrowRequire(d+1!=N+1); // Full depth should not make children
        optionsQueues[child.nMerge].push(child); // Add to queues based on nMerge
        //std::cout << child.nMerge << std::endl;
        //optionsQueues[0].push(child); // Add to queues based on nMerge
      }
    }
  }
  best.runLOP();
}