#include <iostream>
#include "../include/combinations.hpp"

using namespace std;

int main()
{
    AllCombinationGenerator combo_gen(5, 3);
    // AllCombinationGenerator combo_gen(2,2);
    while(!combo_gen.done)
    {
        cout << "(" << combo_gen.combo_cnt << ") ";
        for (auto idx: combo_gen.indices)
        {
            cout << idx << " ";
        }
        cout << endl;
        combo_gen.next();
    }

}