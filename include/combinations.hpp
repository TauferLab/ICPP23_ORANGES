#ifndef COMBINATIONS_HPP
#define COMBINATIONS_HPP

#include <vector>
#include <stdexcept>
//adapted from Nigel's Kokkos branch

// #include "structure_defs.hpp"
// #include "input_to_network.hpp"
// #include "printout_others.hpp"
// #include "printout_network.hpp"
// #include "ADJ/find_Xneighbors.hpp"
// template<typename T>

// gives all combinations from 1..k
class AllCombinationGenerator
{
    public:
    int n;
    int max_k;
    int k;
    size_t combo_cnt;
    bool done;
    std::vector<int> indices;  // combination k indices of n ie if n=5,k=3 => [0,2,4]

    AllCombinationGenerator(int num_elements, int max_size)
    :n(num_elements),max_k(max_size),combo_cnt(1)
    {
        if (n < 1 || max_k > n)
        {
            throw std::invalid_argument("Invalid parameters.");
        }
        k = max_k;
        done = false;
        for(int i=0; i<k; i++) 
        {
            indices.push_back(i);
        }
    }

    bool next()
    {
        done = true;
        for(int i=k-1; i>=0; --i)
        {
            if(indices[i] < n - k + i)
            {
                int j = indices[i] + 1;
                while(i <= k-1)
                    indices[i++] = j++;
                combo_cnt++;
                done = false;
                break;
            }
        }
        if (done && k > 1)
        {
            k--;
            indices.clear();
            for (int i=0; i<k; i++)
            {
                indices.push_back(i);
            }
            combo_cnt++;
            done = false;
        }
        return done;
    }
}; // end class AllCombinationGenerator

#endif
// class Combinations
// {
//     public:
//         // The main function that prints all 
//         // combinations of size r in arr[] 
//         // of size n. This function mainly 
//         // uses combinationUtil() 
//         // output -> 
//         void getCombination(int arr[], int n, int r, std::vector<std::vector<int>>* output)
//         {
//             // A temporary array to store 
//             // all combination one by one 
//             int* data = new int[r];

//             // Print all combination using 
//             // temprary array 'data[]' 
//             combinationUtil(arr, n, r, 0, data, 0, output);

//             delete[] data;
            
//         }
//      private:
//         /* arr[] ---> Input Array
//         n ---> Size of input array
//         r ---> Size of a combination to be printed
//         index ---> Current index in data[]
//         data[] ---> Temporary array to store current combination
//         i ---> index of current element in arr[] */
//         void combinationUtil(int arr[], int n, int r,
//             int index, int data[], int i, std::vector<std::vector<int>>* combilist)
//         {
//             // Current cobination is ready, print it 
//             if (index == r)
//             {
//                 std::vector<int> vec;
//                 for (int j = 0; j < r; j++)
//                 {
//                     //cout << data[j] << " ";
//                     vec.push_back(data[j]);
//                 }
//                 combilist->push_back(vec);
//                 return ;
//             }

//             // When no more elements are there to put in data[] 
//             if (i >= n)
//                 return;
//             // current is included, put next at next location 
//             data[index] = arr[i];
//             combinationUtil(arr, n, r, index + 1, data, i + 1, combilist);
//             // current is excluded, replace it with next (Note that 
//             // i+1 is passed, but index is not changed) 
//             combinationUtil(arr, n, r, index, data, i + 1, combilist);
//         }


// };


// // Class for creating combinations of k indices
// class CombinationGeneratorVector
// {
//   public:
//     std::vector<int> combination; // Indices for combination
//     bool done;
//     int n, k;

//     CombinationGeneratorVector(int num_elements, int combo_size)
//     {
//       n = num_elements;
//       k = combo_size;
//       done =  n<1 || k>n;
//       for(int i=0; i<k; i++) 
//       {
//         combination.push_back(i);
//       }
//     }

//     void get_combination(const std::vector<int>& set, std::vector<int>& combo) 
//     {
//       combo.clear();
//       for(int i=0; i<k; i++) {
//         if(combination[i] < set.size())
//           combo.push_back(set[combination[i]]);
//       }
//       return;
//     }

//     bool next()
//     {
//       done = true;
//       for(int i=k-1; i>=0; --i)
//       {
//         if(combination[i] < n - k + i + 1)
//         {
//           int j = combination[i] + 1;
//           while(i <= k-1)
//             combination[i++] = j++;
//           done = false;
//           break;
//         }
//       }
//       return done;
//     }
// };

// #endif