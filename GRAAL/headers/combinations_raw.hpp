#ifndef COMBINATIONS_RAW_HPP
#define COMBINATIONS_RAW_HPP

#include "structure_defs.hpp"
#include "input_to_network.hpp"
#include "printout_others.hpp"
#include "printout_network.hpp"
#include "ADJ/find_Xneighbors.hpp"
using namespace std;

class Combinations_raw
{
    public:
        // The main function that prints all 
        // combinations of size r in arr[] 
        // of size n. This function mainly 
        // uses combinationUtil() 
        // output -> 
        void getCombination(int arr[], int n, int r, intvecvec& output)
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
            int index, int data[], int i, intvecvec& combilist)
        {
            // Current cobination is ready, print it 
            if (index == r)
            {
                intvec vec = new_intvec(r);
                for (int j = 0; j < r; j++)
                {
                    //cout << data[j] << " ";
                    //vec.push_back(data[j])
		    pushback_intvec(vec, data[j]);
                }
                //combilist->push_back(vec);
		pushback_intvecvec(combilist, vec);
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

#endif
