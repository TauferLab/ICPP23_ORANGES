#ifndef COMBINATIONS_HPP
#define COMBINATIONS_HPP

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
using namespace std;

class Combinations
{
    public:
        // The main function that prints all 
        // combinations of size r in arr[] 
        // of size n. This function mainly 
        // uses combinationUtil() 
        // output -> 
        void getCombination(int arr[], int n, int r, vector<vector<int>>* output)
        {
            // A temporary array to store 
            // all combination one by one 
            int* data = new int[r];

            // Print all combination using 
            // temprary array 'data[]' 
            combinationUtil(arr, n, r, 0, data, 0, output);

            delete[] data;
            
        }
     private:
        /* arr[] ---> Input Array
        n ---> Size of input array
        r ---> Size of a combination to be printed
        index ---> Current index in data[]
        data[] ---> Temporary array to store current combination
        i ---> index of current element in arr[] */
        void combinationUtil(int arr[], int n, int r,
            int index, int data[], int i, vector<vector<int>>* combilist)
        {
            // Current cobination is ready, print it 
            if (index == r)
            {
                vector<int> vec;
                for (int j = 0; j < r; j++)
                {
                    //cout << data[j] << " ";
                    vec.push_back(data[j]);
                }
                combilist->push_back(vec);
                return ;
            }

            // When no more elements are there to put in data[] 
            if (i >= n)
                return;
            // current is included, put next at next location 
            data[index] = arr[i];
            combinationUtil(arr, n, r, index + 1, data, i + 1, combilist);
            // current is excluded, replace it with next (Note that 
            // i+1 is passed, but index is not changed) 
            combinationUtil(arr, n, r, index, data, i + 1, combilist);
        }


};

// Class for creating combinations of k indices
class CombinationGenerator
{
  public:
    vector<int> combination; // Indices for combination
    bool done;
    int n, k;

    CombinationGenerator(int num_elements, int combo_size)
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
        if(combination[i] < set.size())
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

#endif
