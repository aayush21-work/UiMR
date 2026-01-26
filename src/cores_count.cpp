#include <iostream>
#include <thread>

int main() {
    std::cout << "Logical cores: " << std::thread::hardware_concurrency() << std::endl;
    return 0;
}   
