#include <iostream>
#include <future>
#include <chrono>

// Artificially expensive add
long long add(long long a, long long b) {
    long long sum = 0;
    for (long long i = 0; i < 500000000; ++i) {
        sum += (a + b) % 7;
    }
    return sum;
}

// Artificially expensive subtract
long long subtract(long long a, long long b) {
    long long diff = 0;
    for (long long i = 0; i < 500000000; ++i) {
        diff += (a - b) % 7;
    }
    return diff;
}

int main() {
    long long a = 10, b = 5;

    // ---------------- SEQUENTIAL ----------------
    auto start_seq = std::chrono::high_resolution_clock::now();

    auto s1 = add(a, b);
    auto s2 = subtract(a, b);

    auto end_seq = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> time_seq = end_seq - start_seq;

    std::cout << "Sequential results:\n";
    std::cout << "Sum  = " << s1 << "\n";
    std::cout << "Diff = " << s2 << "\n";
    std::cout << "Sequential time = " << time_seq.count() << " seconds\n\n";

    // ---------------- PARALLEL (ASYNC) ----------------
    auto start_par = std::chrono::high_resolution_clock::now();

    auto f1 = std::async(std::launch::async, add, a, b);
    auto f2 = std::async(std::launch::async, subtract, a, b);

    auto p1 = f1.get();
    auto p2 = f2.get();

    auto end_par = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> time_par = end_par - start_par;

    std::cout << "Parallel (async) results:\n";
    std::cout << "Sum  = " << p1 << "\n";
    std::cout << "Diff = " << p2 << "\n";
    std::cout << "Parallel time = " << time_par.count() << " seconds\n\n";

    // ---------------- SPEEDUP ----------------
    std::cout << "Speedup ≈ " << time_seq.count() / time_par.count() << "x\n";
}

