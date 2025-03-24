#include "common.h"
#include "branchAndBound.h"
#include "blackBoxLop.cpp"
#include <algorithm>
#include <cassert>
#include <limits>
#include <utility>

void STK_ThrowRequire(bool condition){ if (condition){ throw "Condition not met";}}

// BUCKET OPTION CLASS, CONT
class BucketingOption {
    public:
      BucketingOption(const BlockNormsViewType & blockNorms, const std::vector<int> & addition_order, const double tMaxWalltime)
      : blockNorms(blockNorms),
      tMaxWalltime(tMaxWalltime)
      {
        N = int(addition_order.size());
        // Initialize buckets with 1 physics each
        for (int physics_number = 0; physics_number < N; ++physics_number) {
          buckets.push_back(std::set<int>{physics_number});
        }
        // No buckets are stable
        nStableBuckets=0;
      }

      // Getters
      [[nodiscard]] double get_cost() const {return cost;} // TODO make sure protection is correct
      [[nodiscard]] double get_nMerge() const {return N-int(buckets.size());}
      [[nodiscard]] int get_nStableBuckets() const {return nStableBuckets;}
      [[nodiscard]] auto get_map() const {STK_ThrowRequire(cost_is_current); return physics_to_block_map;}
      //STK_ThrowRequire(cost_is_current); 
      
      void runLOP(){
        updateOrder(); // RUNS LOP
        updateCost(); // Based on new order
        updateMap(); // Based on new order
        cost_is_current=true;
      }
      std::vector<BucketingOption*> makeChildren(){ // Second returns a list of children sorted by potential
        STK_ThrowRequire(cost_is_current);
        // Make children by greedy merge selection
        // First find L*, the larges lower diagonal term associated with a non-stable bucket
        double L_star=std::numeric_limits<double>::max();
        int a_star;
        for (unsigned long a=nStableBuckets; a < buckets.size(); a++){ // Loop through non-stable only
            for (unsigned long b; b < buckets.size(); b++){
                if (order[a]<=order[b]) {continue;} // Filter out nonLD blocks. (for <>, concider if order=0,1,2,3... a<b is an UD BLOCK)
                double L = bucketSum(a,b);
                if (L < L_star){ L_star = L; a_star = a;}
            }
        }
        // Create a child where a_star is merged to each stable block, tracking the potential of each.
        children = std::vector<BucketingOption>(nStableBuckets+1, *this); // Make copies of self
        std::vector<double> potentials = std::vector<double>(nStableBuckets+1);
        for (unsigned long b; b < nStableBuckets; b++){
            potentials[b] = bucketSum(a_star,b);
            children[b].merge(a_star,b);
        }
        // Create child where a_star gets its own bucket
        children[nStableBuckets].newBucket(a_star);

        // Sort the list by potential
      }
      void merge(int a_star, int b){
        STK_ThrowRequire( a_star > b);
        STK_ThrowRequire( a_star > nStableBuckets - 1);
        STK_ThrowRequire( a_star < int(buckets.size()));
        STK_ThrowRequire( 1 == int(buckets[a_star].size()));
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
        cost_is_current = false;
      }

      // Overload the < opp for BucketingOptions
      bool operator<(const BucketingOption& other) const {
        if (get_nMerge() != other.get_nMerge()){ return get_nMerge() < other.get_nMerge();} // Merge count is first priority
        return cost < other.get_cost(); // Second look to cost
      }
    
    private:
      BlockNormsViewType blockNorms;
      double tMaxWalltime;
      std::vector<std::set<int>> buckets;
      int nStableBuckets;
      std::map<int,int> physics_to_block_map;
      std::vector<int> order;
      double cost;
      int N;
      std::vector<BucketingOption> children;
      bool cost_is_current = false;
      double bucketSum(int a,int b){
        double T;
        for (const auto& i : buckets[a]) { 
            for (const auto& j : buckets[b]){
              T += std::pow(blockNorms(i,j),2);
            }
        }
        return T;
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
        order = solver.order();
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
        std::vector<BucketingOption> branchOptions = std::vector<BucketingOption>(base.get_nStableBuckets(), base); // Fill with copies. Preallocated for speed
        std::vector<double> potentials;
        for (int intoBucket=0; intoBucket <= base.get_nStableBuckets(); intoBucket ++){
          potentials.push_back(branchOptions[intoBucket].potentialOfAddingNextPhysicsToBucket(intoBucket));
        }
        // Sort so highest potential goes first
        std::vector<int> branchPriority=std::vector<int>(N);
        std::iota(branchPriority.begin(), branchPriority.end(), 0); // Just make a vector from 0 to N-1 for now
        std::sort(branchPriority.begin(), branchPriority.end(),[&potentials](const int b1, const int b2){return potentials[b1] > potentials[b2];}); // More potential gets a true and goes first
        // Start with the highest priority branch, if its possible that it 
      }
      
  };