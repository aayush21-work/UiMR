#include <iostream>
#include <future>


int add(int a, int b)
{
  return a+b;
}


int diff(int a, int b)
{
  return a-b;
}



int main()
{
  int a=10;
  int b=5;

  auto f1 = std::async(std::launch::async, add, a, b);
  auto f2 = std::async(std::launch::async, diff, a, b);

  std::cout << "Sum: " << f1.get() << "\n";
  std::cout << "Diff: " << f2.get() << "\n";



}
