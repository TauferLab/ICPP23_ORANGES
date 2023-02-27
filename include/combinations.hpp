#ifndef COMBINATIONS_HPP
#define COMBINATIONS_HPP

#include <Kokkos_Core.hpp>
using namespace std;

//class Combinations
//{
//    public:
//        // The main function that prints all 
//        // combinations of size r in arr[] 
//        // of size n. This function mainly 
//        // uses combinationUtil() 
//        // output -> 
//        void getCombination(int arr[], int n, int r, vector<vector<int>>* output)
//        {
//            // A temporary array to store 
//            // all combination one by one 
//            int* data = new int[r];
//
//            // Print all combination using 
//            // temprary array 'data[]' 
//            combinationUtil(arr, n, r, 0, data, 0, output);
//
//            delete[] data;
//            
//        }
//     private:
//        /* arr[] ---> Input Array
//        n ---> Size of input array
//        r ---> Size of a combination to be printed
//        index ---> Current index in data[]
//        data[] ---> Temporary array to store current combination
//        i ---> index of current element in arr[] */
//        void combinationUtil(int arr[], int n, int r,
//            int index, int data[], int i, vector<vector<int>>* combilist)
//        {
//            // Current cobination is ready, print it 
//            if (index == r)
//            {
//                vector<int> vec;
//                for (int j = 0; j < r; j++)
//                {
//                    //cout << data[j] << " ";
//                    vec.push_back(data[j]);
//                }
//                combilist->push_back(vec);
//                return ;
//            }
//
//            // When no more elements are there to put in data[] 
//            if (i >= n)
//                return;
//            // current is included, put next at next location 
//            data[index] = arr[i];
//            combinationUtil(arr, n, r, index + 1, data, i + 1, combilist);
//            // current is excluded, replace it with next (Note that 
//            // i+1 is passed, but index is not changed) 
//            combinationUtil(arr, n, r, index, data, i + 1, combilist);
//        }
//
//
//};

// Class for creating combinations of k indices
class CombinationGenerator
{
  public:
    Kokkos::View<int*> k_combination;
    bool done;
    int n, k;
    int num_comb;

    KOKKOS_INLINE_FUNCTION CombinationGenerator() {}

//    KOKKOS_INLINE_FUNCTION CombinationGenerator(int num_elements, int combo_size)
//    {
//      n = num_elements;
//      k = combo_size;
//      num_comb = 1;
//      done =  n<1 || k>n;
//      k_combination = Kokkos::View<int*>("Combination indices", k);
//      for(int i=0; i<k; i++) 
//      {
//        k_combination(i) = i;
//      }
//    }

    template<class IndexView>
//    KOKKOS_INLINE_FUNCTION CombinationGenerator(int num_elements, int combo_size, Kokkos::View<int*>& indices)
    KOKKOS_INLINE_FUNCTION CombinationGenerator(int num_elements, int combo_size, IndexView& indices)
    {
      n = num_elements;
      k = combo_size;
      done =  n<1 || k>n;
      for(int i=0; i<k; i++) 
      {
        indices(i) = i;
      }
    }

    KOKKOS_INLINE_FUNCTION void set_n(int num_elements)
    {
      n = num_elements;
    }

    KOKKOS_INLINE_FUNCTION void set_k(int num_selection)
    {
      k = num_selection;
    }

    template<class IndexView, class SetView, class CombinationView>
    KOKKOS_INLINE_FUNCTION void kokkos_get_combination(IndexView& indices, int set_size, const SetView& set, CombinationView& combo)  const
    {
      for(int i=0; i<k; i++) {
        if(i < set_size)
          combo(i) = set(indices(i));
      }
      return;
    }

    template<class IndexView>
    KOKKOS_INLINE_FUNCTION bool kokkos_next(IndexView& indices)
    {
      done = true;
      num_comb += 1;
      for(int i=k-1; i>=0; i--) {
        if(indices(i) < n-k+i) {
          indices(i) += 1;
          for(int j=i+1; j<k; j++) {
            indices(j) = indices(j-1)+1;
          }
          done = false;
          break;
        } 
      }
      return done;
    }

    KOKKOS_INLINE_FUNCTION int get_num_comb() {
      return num_comb;
    }
};

// Class for creating combinations of k indices
class CombinationGeneratorVector
{
  public:
    vector<int> combination; // Indices for combination
    bool done;
    int n, k;

    CombinationGeneratorVector(int num_elements, int combo_size)
    {
      n = num_elements;
      k = combo_size;
      done =  n<1 || k>n;
      for(int i=0; i<k; i++) 
      {
        combination.push_back(i);
      }
    }

    void get_combination(const vector<int>& set, vector<int>& combo) 
    {
      combo.clear();
      for(int i=0; i<k; i++) {
        if(static_cast<uint32_t>(combination[i]) < set.size())
          combo.push_back(set[combination[i]]);
      }
      return;
    }

    bool next()
    {
      done = true;
      for(int i=k-1; i>=0; --i)
      {
        if(combination[i] < n - k + i + 1)
        {
          int j = combination[i] + 1;
          while(i <= k-1)
            combination[i++] = j++;
          done = false;
          break;
        }
      }
      return done;
    }
};

KOKKOS_INLINE_FUNCTION
int get_node_from_combinations(const matrix_type& graph, Kokkos::View<uint64_t*>::HostMirror& num_combinations, uint64_t k_interval, uint64_t i) {
  uint64_t threshold = i*k_interval;
  uint64_t counter = 0;
  for(uint32_t j=0; j<num_combinations.extent(0); j++) {
    counter += num_combinations(j);
    if(counter > threshold) {
      return j;
    }
  }
  return graph.numRows();
}

KOKKOS_INLINE_FUNCTION
int get_node_from_combinations(const matrix_type& graph, const Kokkos::View<uint64_t*>& num_combinations, const uint64_t k_interval, const uint64_t i) {
  uint64_t threshold = i*k_interval;
  uint64_t counter = 0;
  for(uint32_t j=0; j<num_combinations.extent(0); j++) {
    counter += num_combinations(j);
    if(counter > threshold) {
      return j;
    }
  }
  return graph.numRows();
}

#endif
