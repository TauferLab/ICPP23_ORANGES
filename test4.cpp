#include<iostream>
#include<cstdlib>
#include<vector>
using namespace std;

class graph_combination
{
   public:
   int parent;
   int child;
};
void update(graph_combination a,vector<graph_combination> &b) 
{
    a.parent = 100;
    a.child = 99;
    b.push_back(a);
}

int main()
{
    graph_combination a,c;
    //vector<graph_combination> b;
    cin>>a.parent;
    cin>>a.child;
    cout<<"a parent is"<<a.parent<<endl;
    cout<<"a child is"<<a.child<<endl;
    c=a;
    c.parent = 99;
    c.child = 100;
    //update(a,b);
    cout<<"a parent is"<<a.parent<<endl;
    cout<<"a child is"<<a.child<<endl;
    cout<<"****************************"<<endl;
    cout<<"c parent is"<<c.parent<<endl;
    cout<<"c child is"<<c.child<<endl;
return 0;    
}