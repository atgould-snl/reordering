#include "common.h"
#include "branchAndBound.h"
#include "blackBoxLop.cpp"
#include <algorithm>
#include <cassert>
#include <limits>
#include <numeric>
#include <utility>

void STK_ThrowRequire(bool condition){ if (condition){ throw "Condition not met";}}

// BUCKET OPTION CLASS, CONT
class BucketingOption {
    public:
      BucketingOption(const BlockNormsViewType & blockNorms, const BlockNormsViewType & blockMergeCosts, const std::vector<int> & addition_order, const double tMaxWalltime)
      : blockNorms(blockNorms),
      blockMergeCosts(blockMergeCosts), // Merge cost will be added from both i,j and j,i when i and j merge
      tMaxWalltime(tMaxWalltime)
      {
        N = int(addition_order.size());
        // Initialize buckets with 1 physics each
        for (int physics_number = 0; physics_number < N; ++physics_number) {
          buckets.push_back(std::set<int>{physics_number});
        }
        // No buckets are stable
        nStableBuckets=0;
        runLOP();
      }

      // Getters
      [[nodiscard]] double get_cost() const {return cost;} // TODO make sure protection is correct
      [[nodiscard]] double get_nMerge() const {return N-int(buckets.size());}
      [[nodiscard]] int get_nStableBuckets() const {return nStableBuckets;}
      [[nodiscard]] auto get_map() const {STK_ThrowRequire(cost_is_current); return physics_to_block_map;}
      //STK_ThrowRequire(cost_is_current); 

      std::vector<BucketingOption> makeChildren() const{ // List of children with potentials, last child is no merge
        STK_ThrowRequire(cost_is_current);
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
      void merge(int a_star, int b){
        STK_ThrowRequire( a_star > b);
        STK_ThrowRequire( a_star > nStableBuckets - 1);
        STK_ThrowRequire( a_star < int(buckets.size()));
        STK_ThrowRequire( 1 == int(buckets[a_star].size()));
        potential = order[a_star]>order[b] ? bucketSum(a_star, b) : bucketSum(b, a_star);
        buckets[b].insert(buckets[a_star].begin(), buckets[a_star].end());
        buckets.erase(buckets.begin()+a_star);
        cost_is_current = false;
      }
      void newBucket(int a_star){ // but a_star into the stable buckets
        STK_ThrowRequire( a_star > nStableBuckets - 1);
        STK_ThrowRequire( a_star < int(buckets.size()));
        STK_ThrowRequire( 1 == int(buckets[a_star].size()));
        std::set<int> tempBucket = buckets[nStableBuckets];
        buckets[nStableBuckets]=buckets[a_star];
        buckets[a_star]=tempBucket; // Ok if a_star = nStableBuckets
        nStableBuckets++;
        // Assign potential of the highest next availible merge by finding a new a_star
        auto [a_star_next, L_star] = get_a_star_L_star();
        potential = L_star;
        cost_is_current = false;
      }
      // Overload the < opp for BucketingOptions
      bool operator<(const BucketingOption& other) const {
        if (get_nMerge() != other.get_nMerge()){ return get_nMerge() < other.get_nMerge();} // Merge count is first priority
        return cost < other.get_cost(); // Second look to cost
      }
      void hasSufficientLB(BucketingOption best){
        // Make a copy of self
        // On the copy, play out optomistic merging until merge counts are the same
        // Then compare

      }
    
    private:
      BlockNormsViewType blockNorms;
      BlockNormsViewType blockMergeCosts; // All 1s for now
      double tMaxWalltime;
      std::vector<std::set<int>> buckets;
      int nStableBuckets;
      std::map<int,int> physics_to_block_map;
      std::vector<int> order_of_rows;
      std::vector<int> order;
      double cost;
      double mergeCost;
      double potential;
      int N;
      bool cost_is_current = false;

      void runLOP(){
        updateOrder(); // RUNS LOP
        updateCost(); // Based on new order
        updateMergeCost(); // Based on bucketing only
        updateMap(); // Based on new order
        cost_is_current=true;
      }
      double bucketSum(int a, int b) const {
        double T;
        for (const auto& i : buckets[a]) { 
            for (const auto& j : buckets[b]){
              T += std::pow(blockNorms(i,j),2);
            }
        }
        return T;
      }
      std::pair<int,double> get_a_star_L_star() const{
        double L_star=std::numeric_limits<double>::max();
        int a_star=0;
        for (unsigned long a=nStableBuckets; a < buckets.size(); a++){ // Loop through non-stable only
            for (unsigned long b; b < buckets.size(); b++){
                if (order[a]<=order[b]) {continue;} // Filter out nonLD blocks. (for <>, concider if order=0,1,2,3... a<b is an UD BLOCK)
                double L = bucketSum(a,b);
                if (L < L_star){ L_star = L; a_star = a;}
            }
        }
        return std::pair<int,double>(a_star,L_star);
      }
      void updateOrder(){
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
        LinearOrderingSolver solver{tournament, tMaxWalltime}; // Is communicator pulled in from the namespace??
        solver.solve();
        order_of_rows = solver.order();
        std::vector<int> order(order_of_rows.size());
        std::iota( order_of_rows.begin(), order_of_rows.end(), 0);
        std::sort(order.begin(), order.end(), [this](int a, int b) { return order_of_rows[a] < order_of_rows[b];});
      }
  
      void updateMap(){
        for (unsigned long bucket_rank=0; bucket_rank<buckets.size(); bucket_rank++){
          auto bucket = buckets[order[bucket_rank]];
          for (const auto& physics_block : bucket) {
            physics_to_block_map[physics_block] = int(bucket_rank);
          }
        }
      }
  
      void updateCost(){ // Find the score based on the order
        cost=0;
        for (unsigned long a=0; a < buckets.size(); a++){
          for (unsigned long b=0; b < buckets.size(); b++){
            cost += order[a]>order[b] ? bucketSum(a,b) : 0;
          }
        }
      }

      void updateMergeCost(){ // Find the score based on the order
        mergeCost=0;
        for (unsigned long a=0; a < buckets.size(); a++){
          for (const auto& i : buckets[a]) {
            for (const auto& j : buckets[a]) {
              if (i!=j){
                mergeCost += blockMergeCosts(i,j);
              }
            }
          }
        }
      }
  };
  
  class BucketOrderingSolver {
    public:
      BucketOrderingSolver(const BlockNormsViewType & blockNorms, double tMaxWalltime)
      : blockNorms(blockNorms),
      tMaxWalltime(tMaxWalltime)
      {
        N = int(blockNorms.extent(0));
      }
      void solve(){
        std::vector<int> memberPresort=std::vector<int>(N);
        std::iota(memberPresort.begin(), memberPresort.end(), 0); // Just make a vector from 0 to N-1 for now
  
        BucketingOption base = BucketingOption(blockNorms,memberPresort, tMaxWalltime);
        addToOrder(base);
      }
    private:
      BlockNormsViewType blockNorms;
      double tMaxWalltime;
      int N;
      void addToOrder(const BucketingOption & base){
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
      
  };