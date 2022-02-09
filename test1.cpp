#include<cstdlib>
#include<iostream>
#include <algorithm>
using namespace std;

class graph_combination
{
   public:
   int parent;
   int child;
};

class subgraph
{
    public:
    int root;
    int no_of_nodes; 
    vector<graph_combination> list;
    vector<int> list_of_nodes;
    vector<int> end_nodes;
    bool is_subgraph;
};

void trace_path_till_root(int root, vector<graph_combination> &list, int iterator)
{
    int current_node = list[iterator].parent;
    while(current_node! = root)
    {
        
    }
}

int main()
{
    int no_of_nodes;
    vector<graph_combination> list;
    cin >> no_of_nodes;
    int root;
    cin>> root;
    graph_combination parent_child_comb;
    for(int i = 0; i < no_of_nodes; i++)
    {
        cin>>parent_child_comb.parent;
        cin>>parent_child_comb.child;
        list.push_back(parent_child_comb);
    }
    if(no_of_nodes >= 3)
    {
        int iterator = 2;
        while(iterator < no_of_nodes)
        {
            trace_path_till_root(root,list,iterator);
            iterator++;
        }    
    }

    return 0;
}